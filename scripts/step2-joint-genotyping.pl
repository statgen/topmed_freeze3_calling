#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef);
use gcconfig;

&initRef($ref);

my @chrs = @ARGV;

if ( $#chrs < 0 ) {
    die "Usage: [command] [list of chromosomes separated by space]\n";
}

my $outDir = "$out/paste";
my $unionDir = "$out/aux/union";

mkpath($outDir) unless ( -e $outDir ); # || die "Cannot! create directory $outDir";

my $outf = "$outDir/chr".join("_",@chrs);
open(MAK,">$outf.Makefile") || die "Cannot open file $outf.Makefile for writing\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all:";
foreach my $chr (@chrs) {
    print MAK " $outDir/$chr.OK";
}
print MAK "\n\n";

my $unit = $genotypeUnit;

foreach my $chr (@chrs) {
    unless ( -e "$outDir/$chr" ) {
	mkdir("$outDir/$chr") || die "Cannot create directory $outDir/$chr\n";
    }

    my $szchr = $hszchrs{$chr}->[3];
    my @cmds = ();
    my @tgts = ();
    my @deps = ();    
    my @outvcfs = ();

    for(my $iD=1; $iD < $szchr; $iD += $discoverUnit) {
	my $begD = $iD;
	my $endD = ( $iD + $discoverUnit > $szchr ) ? $szchr : ($iD + $discoverUnit - 1);
	my $sitevcf = "$unionDir/$chr\_$begD\_$endD.sites.bcf";
	
	for(my $i=$begD; $i < $endD; $i += $unit) {
	    my $beg = $i;
	    my $end = ( $i + $unit > $szchr ) ? $szchr : $i + $unit - 1;
	    my $outvcf = "$outDir/$chr/$chr\_$beg\_$end\_paste.bcf";
	    my $cmd = "REF_PATH=$md5 $vt joint_genotype_sequential -r $ref -L $index -i $chr:$beg-$end -o $outvcf $sitevcf 2> $outvcf.log && $vt index $outvcf && touch $outvcf.OK"; 
	    push(@cmds,$cmd);
	    push(@outvcfs,$outvcf);
	    push(@tgts,"$outvcf.OK");
	    push(@deps,"$sitevcf");
	}
    }
    my $svcf = "$outDir/$chr\_1\_$szchr\_paste.sites.vcf.gz";

    print MAK "$outDir/$chr.OK : @tgts\n";
    print MAK "\t($bcftools view -G $outvcfs[0]";
    for(my $i=1; $i < @outvcfs; ++$i) {
	print MAK "; $bcftools view -G -H $outvcfs[$i]";
    }
    print MAK ") | $bgzip -c > $svcf\n";
    print MAK "\t$tabix -pvcf $svcf\n";    
    print MAK "\ttouch $outDir/$chr.OK\n\n";
    
    for(my $j=0; $j < @tgts; ++$j) {
	print MAK "$tgts[$j]: $deps[$j]\n";
	print MAK "\t$cmds[$j]\n\n";
    }
}
close MAK;

print "Run make -f $outf.Makefile -j [numjobs] to complete this step\n";
