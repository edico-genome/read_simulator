#!/usr/bin/perl

use strict; use warnings;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname abs_path $0) . '/AltAware';
use AltAware::Coords qw(alt_to_pri_pos);

my ($pri_name, $pri_start) = alt_to_pri_pos("/opt/bwakit-0.7.12-0/resource-GRCh38/hs38DH.fa.alt", "chr6_GL000251v2_alt", 10);

print "$pri_name:$pri_start\n";
