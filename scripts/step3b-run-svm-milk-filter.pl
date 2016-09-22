#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef forkExecWait);
use gcconfig;

&initRef($ref);

my @chrs = @ARGV;

if ( $#chrs < 0 ) {
    die "Usage: [command] [list of chromosomes separated by space]\n";
}

my $milkDir = "$out/aux/milk";
my $svmDir = "$out/svm";
my $pasteDir = "$out/paste";
my @posVcfs = ($omnivcf,$hapmapvcf);
my @ignores = qw(AC AN AF GC GN HWEAF HWDGF);
my @includes = ();
my $checkNA = 1;
my $cutoff = -0.5;
my $bfIntercept = 2;
my $bfSlope = 2;

unless ( -e "$svmDir" ) {
    mkdir("$svmDir") || die "Cannot create directory $pasteDir\n";
}

my @names = ();
my $ncols = 0;
my %hIgnores = map { $_ => 1 } @ignores;
my %hIncludes = ();
my @includeKeys = ();
my @includeDefaultValues = ();
for(my $i=0; $i < @includes; ++$i) {
    my ( $key, $defaultVal ) = split(/,/,$includes[$i]);
    $hIncludes{$key} = $i;
    push(@names,$key);
    push(@includeKeys,$key);
    push(@includeDefaultValues,defined($defaultVal) ? $defaultVal : 0);
}
my $nIncludes = $#includes + 1;


foreach my $chr (@chrs) {
    my $szchr = $hszchrs{$chr}->[3];        
    my $out = "$svmDir/$chr\_1\_$szchr\_milk_svm";

    my %hpos = ();

    print STDERR "Loading positive labels..\n";
    ## LABEL positive samples
    foreach my $posVcf (@posVcfs) {
	open(IN,"$tabix $posVcf $chr:1-$szchr| grep -w PASS|") || die "Cannot open file\n";
	while(<IN>) {
	    next if ( /^#/ );
	    my ($chrom,$pos,$id,$ref,$alt) = split;
	    #next if ( $pos > $szchr/2 ); ## temporary for CV
	    #$hpos{"$chrom:$pos:$ref:$alt"} = 1;
	    $hpos{"$chrom:$pos"} = "$ref:$alt";
	}
    }

    my %hneg = ();
    my $milkVcf = "$milkDir/$chr\_1\_$szchr\_milk.sites.vcf.gz";
    print STDERR "Identifying negative labels from $milkVcf using intercept $bfIntercept, log10-scale slope $bfSlope..\n";        

    open(IN,"zcat $milkVcf | cut -f 1-8| grep -v ^#|") || die "Cannot open file\n";
    while(<IN>) {
	my @F = split;
	next if ( $F[4] =~ /,/ );
	next if ( ( ( length($F[3]) > 1 ) || ( length($F[4]) > 1 ) ) && ( $F[6] ne "." ) );

	my $ac = $1 if ( $F[7] =~ /AC=(\d+)/ );
	my $bf = $1 if ( $F[7] =~ /MILK_BF=([^;]+)/ );
	my $abe = $1 if ( $F[7] =~ /ABE=([^;]+)/ );
	my $stz = $1 if ( $F[7] =~ /STZ=([^;]+)/ );
	my $fic = $1 if ( $F[7] =~ /IBC=([^;]+)/ );
	my $ior = $1 if ( $F[7] =~ /IOR=([^;]+)/ );
	my $cyz = $1 if ( $F[7] =~ /CYZ=([^;]+)/ );
	my $qual = $F[5];
	my @tcnts = split(/,/,$1) if ( $F[7] =~ /TRIO_CONC_THRES=([^;]+)/ );
	my @dcnts = split(/,/,$1) if ( $F[7] =~ /DUP_CONC_THRES=([^;]+)/ );	
	my @discords = (0,1,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,1,0);
	my @dupdisc = (0,1,1,1,0,1,1,1,0);	
	my ($ndisc,$nconc,$ddisc,$dconc) = (0,0,0,0);
	for(my $i=1; $i < 27; ++$i) {
	    if ( $tcnts[$i] > 0 ) {
		if ( $discords[$i] == 0 ) { $nconc += $tcnts[$i] }
		else { $ndisc += $tcnts[$i]; }
	    }
	}

	for(my $i=1; $i < 9; ++$i) {
	    if ( $dcnts[$i] > 0 ) {
		if ( $dupdisc[$i] == 0 ) { $dconc += $dcnts[$i] }
		else { $ddisc += $dcnts[$i]; }		
	    }
	}		

	my $nfail = 0;
	++$nfail if ( $abe > 0.7 );
	++$nfail if ( abs($stz) > 5 );
	++$nfail if ( $fic < -0.1 );
	++$nfail if ( $ior > 2 );
	++$nfail if ( $cyz < -5 );
	++$nfail if ( $qual < 5 );

	next if ( $ac < 2 );  ## do not include singletons
	#if ( ( $bf < 0 - $bfIntercept - log($ac+1)/log(10)*$bfSlope ) || ( $nfail > 2 ) ) {
	if ( ( $bf < 0 - $bfIntercept - log($ac+1)/log(10)*$bfSlope ) || ( ( $ndisc > 2 ) && ( $ndisc > 0.2 * $nconc ) ) || ( ( ( $ddisc > 1 ) && ( $ddisc > 0.2 * $dconc ) ) ) ) {
	    $hneg{"$F[0]:$F[1]:$F[3]:$F[4]"} = 1;
	}	
    }
    close IN;

    print STDERR "Writing the feature information for libsvm..\n";
    open(RAW,">$out.raw") || die "Cannot open $out.raw for writing\n";
    open(SITE,">$out.site") || die "Cannot open $out.raw for writing\n";
    open(IN,"zcat $pasteDir/$chr\_1\_$szchr\_paste.sites.vcf.gz|") || die "Cannot open file\n";
    my @hdrs = ();
    while(<IN>) {
	if ( /^#/ ) {
	    push(@hdrs,$_);
	    next;
	}
	my ($chrom,$pos,$id,$ref,$alt,$qual,$filt,$info) = split(/[ \t\r\n]+/);
	$info .= ";QUAL=$qual" unless ( $qual eq "." );
	my @infos = split(/;/,$info);
	my @values = ();
	my $k = 0;
	for(my $j=0; $j < @includeKeys; ++$j) {
	    push(@values,$includeDefaultValues[$j]);
	    ++$k;
	}	    
	
	for(my $j=0; $j < @infos; ++$j) {
	    my ($key,$val) = split(/=/,$infos[$j]);
	    #next if ( defined($hIgnores{$key})  || (!defined($hIncludes{$key}))); ## skip if ignored, or not in includes key
	    next if ( defined($hIgnores{$key}) );
	    next unless defined($val); ## skip keys without any values
	    if ( defined($hIncludes{$key}) ) {
		if ( !($checkNA) || ($val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) ) {
		    $values[$hIncludes{$key}] = $val; # set value if given
		}
	    }
	    else {
		if ( !($checkNA) || ( $val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) ) {
		    push(@values,$val);
		}
		else {
		    push(@values,0);
		}
		
		if ( $ncols == 0 ) {
		    push(@names,$key);
		}
		else {
		    die "Cannot recognize $key in $infos[$j], supposed to be $names[$j] at $j, $chrom:$pos:$ref:$alt $info\n" unless ($names[$k] eq $key );
		}
		++$k;
	    }
	}
	if ( $ncols == 0 ) {
	    $ncols = $#names+1;
	    print STDERR "Recording following $ncols features : @names\nThe info field was $info\n@infos\n";
	    die if ( $ncols == 0 );
	}
	elsif ( $ncols != $#values+1 ) {
	    die "Number of columns are not identical at $chrom:$pos:$ref:$alt\n";
	}
	print SITE join("\t",$chrom,$pos,$id,$ref,$alt,$qual,$filt,$info)."\n";	    
	print RAW join("\t",@values)."\n";
    }
    close IN;
    close RAW;
    close SITE;

    print STDERR "Performing quantile normalization of features..\n";    
    my $cmd = "$invNorm --in $out.raw --out $out.norm";
    &forkExecWait($cmd);

    open(SITE,"$out.site") || die "Cannot open $out.sites\n";
    open(NORM,"$out.norm") || die "Cannot open $out.norm\n";
    open(FTR,">$out.feature") || die "Cannot open $out.feature\n";
    open(LBL,">$out.label") || die "Cannot open $out.label\n";

    my ($npos,$nneg,$noth) = (0,0,0);
    while(<SITE>) {
	my ($chrom,$pos,$id,$ref,$alt) = split;
	my @z = split(/[ \t\r\n]+/,<NORM>);
	my $ln = "";
	for(my $i=0; $i < @names; ++$i) {
	    $ln .= " ".($i+1).":$z[$i]";
	}
	#if ( defined($hpos{"$chrom:$pos:$ref:$alt"}) ) {
	if ( defined($hpos{"$chrom:$pos"}) ) {	    
	    if ( ( $hpos{"$chrom:$pos"} ne "$ref:$alt" ) || defined($hneg{"$chrom:$pos:$ref:$alt"}) ) {
		print FTR "0 $ln\n";
		++$noth;
	    }
	    else {
		print FTR "1 $ln\n";
		print LBL "1 $ln\n";
		++$npos;
	    }
	}
	elsif ( defined($hneg{"$chrom:$pos:$ref:$alt"}) ) {
	    print FTR "-1 $ln\n";
	    print LBL "-1 $ln\n";
	    ++$nneg;
	}
	else {
	    print FTR "0 $ln\n";
	    ++$noth;
	}
    }
    close FTR;
    close LBL;
    close SITE;

    print STDERR "Training SVM classification model using $npos positive samples, $nneg negative samples with $noth additional samples\n";

    $cmd = "$svmlearn -s 0 -t 2 $out.label $out.svm.model";
    &forkExecWait($cmd);

    print STDERR "Applying trained SVM model on $noth additional samples\n";

    $cmd = "$svmclassify $out.feature $out.svm.model $out.svm.pred";
    &forkExecWait($cmd);

    print STDERR "Writing filtered site VCF files with SVM scores..\n";
    open(SITE,"$out.site") || die "Cannot open $out.sites\n";    
    open(PRED,"$out.svm.pred") || die "Cannot open $out.svm.pred file\n";
    open(OUT,"| $bgzip -c >$out.sites.vcf.gz") || die "Cannot open $out.sites.vcf.gz\n";
    print OUT join("",@hdrs);

    my %hbest = ();
    
    while(<SITE>) {
	my @F = split;
	my $pred = <PRED>;
	chomp($pred);

	my @filts = ( ($F[6] eq "PASS") || ( $F[6] eq "." ) ) ? () : split(/;/,$F[6]);
	if ( $pred < $cutoff ) {
	    push(@filts,"SVM");
	}
	$F[6] = ($#filts < 0) ? "PASS" : join(';',@filts);
	$F[7] .= ";SVM=$pred";
	
	print OUT join("\t",@F)."\n";

	#next if ( length($F[3]) + length($F[4]) > 2 );
	if ( defined($hbest{$F[1]}) ) {  ## keep best SNPs among possible ones
	    if ( $hbest{$F[1]}->[2] < $pred ) {
		$hbest{$F[1]} = [$F[3],$F[4],$pred];		
	    }
	}
	else {
	    $hbest{$F[1]} = [$F[3],$F[4],$pred];
	}	
    }
    close OUT;
    close PRED;
    close SITE;

    #my $ns = 10597;
    #my $acbreaks = join(",",1,2,3,int($ns*0.002),int($ns*0.02),int($ns*0.2));

    print STDERR "Producing VCF summary for SNPs..\n";
    open(IN,"zcat $out.sites.vcf.gz|") || die "Cannot open file\n";
    #open(OUT1, "| perl $vcfsummary --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr --info AC --info-breaks $acbreaks > $out.sites.summary_ac") || die "Cannot open file\n";
    open(OUT2, "| perl $vcfsummary --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr > $out.sites.summary") || die "Cannot open file\n";
    open(OUT3, "| perl $vcfsummary2 --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr > $out.sites.summary_v2") || die "Cannot open file\n";    
    while(<IN>) {
	next if ( /^#/ );
	my @F = split;
	#next if ( length($F[3]) + length($F[4]) > 2 );
	next unless ( ( defined($hbest{$F[1]}) ) && ( $hbest{$F[1]}->[0] eq $F[3] ) && ( $hbest{$F[1]}->[1] eq $F[4] ) );
	#print OUT1 $_;
	print OUT2 $_;
	print OUT3 $_;
    }
    close IN;
    #close OUT1;
    close OUT2;
    close OUT3;
}
