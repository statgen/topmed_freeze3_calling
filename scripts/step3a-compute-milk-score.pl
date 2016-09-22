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

my $milkDir = "$out/aux/milk";
my $pasteDir = "$out/paste";
my $unit = $genotypeUnit;

unless ( -e "$milkDir" ) {
    mkpath("$milkDir") || die "Cannot create directory $pasteDir\n";
}

my $outf = "$milkDir/milk.$chrs[0].$chrs[$#chrs]";

open(MAK,">$outf.Makefile") || die "Cannot open file $outf.Makefile for writing\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all:";
foreach my $chr (@chrs) {
    my $szchr = $hszchrs{$chr}->[3];    
    my $catvcf = "$milkDir/$chr\_1\_$szchr\_milk";
    print MAK " $catvcf.sites.vcf.gz.tbi"
}
print MAK "\n\n";

foreach my $chr (@chrs) {
    unless ( -e "$milkDir/$chr" ) {
	mkdir("$milkDir/$chr") || die "Cannot create directory $pasteDir\n";
    }
    
    my $szchr = $hszchrs{$chr}->[3];
    my @cmds = ();
    my @tgts = ();
    my @milkVcfs = ();

    for(my $i=1; $i < $szchr; $i += $unit) {
	my $beg = $i;
	my $end = ( $i + $unit > $szchr ) ? $szchr : $i + $unit - 1;
	my $pasteBcf = "$pasteDir/$chr/$chr\_$beg\_$end\_paste.bcf";
	my $milkVcf = "$milkDir/$chr/$chr\_$beg\_$end\_milk.vcf.gz";
	my $cmd = "$vt milk_filter -f $pedf -b $pasteBcf -o $milkVcf && $tabix -pvcf $milkVcf";
	push(@cmds,$cmd);
	push(@milkVcfs,$milkVcf);
	push(@tgts,"$milkVcf.tbi");
    }

    my $catvcf = "$milkDir/$chr\_1\_$szchr\_milk";

    print MAK "$catvcf.sites.vcf.gz.tbi: @tgts\n";
    print MAK "\t($tabix -h $milkVcfs[0] NA:0; zcat @milkVcfs | grep -v ^#;) | $bgzip -c > $catvcf.full.vcf.gz\n";
    print MAK "\t$tabix -pvcf $catvcf.full.vcf.gz\n";    
    print MAK "\tzcat $catvcf.full.vcf.gz | cut -f 1-8 | $bgzip -c > $catvcf.sites.vcf.gz\n";    
    print MAK "\t$tabix -pvcf $catvcf.sites.vcf.gz\n\n";
    
    for(my $j=0; $j < @tgts; ++$j) {
	print MAK "$tgts[$j]:\n";
	print MAK "\t$cmds[$j]\n\n";
    }    
}
close MAK;

print "Run make -f $outf.Makefile -j [numjobs] to complete this step\n";
