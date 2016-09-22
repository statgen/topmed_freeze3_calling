package gcconfig;

use base qw/Exporter/;
use Cwd qw(realpath);
use File::Basename qw(dirname);
use POSIX qw(pow sqrt);
use FindBin;

## Variables and methods shared across the package
our @EXPORT = qw($ref $md5 $bgzip $tabix);

############################################################
### MODIFY THESE VARIABLES TO YOUR COMPUTING ENVIRONMENT
our $md5 = "/data/local/ref/gotcloud.ref/md5/%2s/%s/%s";
our $ref = "/data/local/ref/gotcloud.ref/hs37d5.fa";
our $index = "data/trio_data.index";
############################################################

our $bgzip = "$FindBin::Bin/../cramore/lib/htslib/bgzip";
our $tabix = "$FindBin::Bin/../cramore/lib/htslib/tabix";
our

1;
