package gcconfig;

use base qw/Exporter/;
use Cwd qw(realpath);
use File::Basename qw(dirname);
use POSIX qw(pow sqrt);
use FindBin;

## Variables and methods shared across the package
our @EXPORT = qw($ref $md5 $bgzip $samtools $bcftools $bamUtil $tabix $index $pedf $out $vt $discoverUnit $genotypeUnit $omnivcf $hapmapvcf $dbsnp $invNorm $svmlearn $svmclassify $vcfsummary $vcfsummary2);

############################################################
### MODIFY THESE VARIABLES TO YOUR COMPUTING ENVIRONMENT
our $index = "data/trio_data.index";
our $pedf = "data/trio_data.ped";
our $out = "out";
our $discoverUnit = 20000000;
our $genotypeUnit = 1000000;
############################################################
### MODIFY THESE VARIABLES TO IF REFERENCE IS LOCATED ELSEWHERE
our $refDir = "$FindBin::Bin/../gotcloud.ref";
our $md5 = "$refDir/md5/%2s/%2s/%s";
our $ref = "$refDir/hs37d5.fa";
our $dbsnp = "$refDir/dbsnp_142.b37.vcf.gz";
our $hapmapvcf = "$refDir/hapmap_3.3.b37.sites.vcf.gz";
our $omnivcf = "$refDir/1000G_omni2.5.b37.sites.PASS.vcf.gz";
############################################################
### MODIFY THESE VARIABLES TO IF EXTERNAL BINARIES ARE USED
our $bgzip = "$FindBin::Bin/../htslib/bgzip";
our $tabix = "$FindBin::Bin/../htslib/tabix";
our $vt = "$FindBin::Bin/../vt/vt";
our $samtools = "$FindBin::Bin/../samtools/samtools";
our $bcftools = "$FindBin::Bin/../bcftools/bcftools";
our $bamUtil = "$FindBin::Bin/../gotcloud/src/bin/bamUtil";
our $svmlearn = "$FindBin::Bin/../gotcloud/src/bin/svm-train";
our $svmclassify = "$FindBin::Bin/../gotcloud/src/bin/svm-predict";
our $invNorm = "$FindBin::Bin/../gotcloud/src/bin/invNorm";
our $vcfsummary = "$FindBin::Bin/vcf-summary";
our $vcfsummary2 = "$FindBin::Bin/vcf-summary-v2";

1;
