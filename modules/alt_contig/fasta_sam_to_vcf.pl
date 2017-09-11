#!/usr/bin/perl
#
# Usage: fasta_sam_to_vcf.pl FASTA ALT_SAM SEQ_RANGE1 SEQ_RANGE2 > OUT.VCF
#
# This tool requires a single FASTA, an ALT contig SAM, and 2 overlapping sequence 
# ranges (i.e. haplotypes). SEQ_RANGE1 must be on an ALT contig. On the command line,
# the two ranges are specified like this:
#        chr1_GL383518v1_alt:17254-17263 chr1:153700531-153700540
# The ALT SAM file provides alignments of ALT contigs to the primary assembly contigs.
#
# This tool uses `samtools faidx $input_fasta $chr_range` to fetch sequences. The tool
# requires that the two input ranges align to exactly the same region on the primary
# assembly contig. The first range arg must be fully specified; i.e. CONTIG:START-END.
# If the second input range arg consists of just a contig name (lacks base pos coords), 
# then the range on this contig corresponding to the first range arg will be determined
# automatically. For instance, if user specifies:
#     chr1_GL383518v1_alt:17254-17263 chr1
# for the two range args then the tool will find the range on chr1 corresponding to
# chr1_GL383518v1_alt:17254-17263.
# 
# The FASTA must include each primary assembly contig and ALT contig specified on the
# command line. The FASTA must also include the primary assembly contigs given in the
# RNAME fields in the ALT SAM file. The ALT SAM file must contain primary assembly 
# alignment(s) for each ALT contig specified on the command line.
#
# This tool may fail, or produce ill-defined results, if an input ALT contig range does
# not begin in a 'M','=', or 'X' operation in the ALT's CIGAR. In future, the tool may
# find the nearest range that satisfies this requirement. 
#
# The tool's output is VCF format w/ one record per variant starting POS. Variant alleles
# are in the ALT contig(s) w.r.t. the corresponding primary assembly region. The variants 
# are phased and zygosity is provided by the VCF genotype value.

use strict; use warnings;

use v5.10.1;

use File::Basename;
use List::MoreUtils qw(uniq);

my $dname = dirname($0);
my $CONVERT_ALT_TO_PRI = "$dname/convert_alt_to_pri_coords.pl";
my $CONVERT_PRI_TO_ALT = "$dname/convert_pri_to_alt_coords.pl";

my $MAX_STD_INDEL = 50;

# The following works for hg38.
my @pri_contig_list = 1..22;
@pri_contig_list = map { "chr" . $_ } @pri_contig_list;
push @pri_contig_list, ("chrM", "chrX", "chrY");

sub usage {

   print "Usage:\n  fasta_seqs_to_vcf.pl your.fasta alt_contig.sam seq_range1 seq_range2 > output.vcf\n";
   return;
}

if (@ARGV != 4) {
   usage();
   exit(1);
}

my ($fasta, $alt_sam, $range1, $range2) = @ARGV;


# This subroutine requires that $alt_sam already be set.
sub prim_contig_range {

   my $seq_name = $_[0];
   my $seq_start = $_[1];
   my $seq_end = $_[2];

   my $seq_pri_contig_name;
   my $seq_pri_contig_start;
   my $seq_pri_contig_end;

   if ( grep( /^$seq_name$/, @pri_contig_list ) ) {
       $seq_pri_contig_name = $seq_name;
       $seq_pri_contig_start = $seq_start;
       $seq_pri_contig_end = $seq_end;
   } elsif ($seq_name =~ /_alt$/ or $seq_name =~ /^HLA/) {

       # This is an ALT contig. Find primary assembly region overlapping this region.
       if ($seq_start < 1 or  $seq_end < 1) {
           die "Start and end positions must be specified for an ALT contig."
       }
       my $ret_str = `$CONVERT_ALT_TO_PRI $alt_sam $seq_name $seq_start $seq_end`;
       chomp $ret_str;
       my @field = split "\t", $ret_str;
       $seq_pri_contig_name = $field[0];
       $seq_pri_contig_start = $field[1];
       $seq_pri_contig_end = $field[2];
   } else {
       die "$seq_name is not a valid contig name for this tool.";
   }

   return ($seq_pri_contig_name, $seq_pri_contig_start, $seq_pri_contig_end);
}

sub cigar_array {

    my @result = ();
    my ($cigar) = @_;
    my @op_arr = $cigar =~ /[\d]+[MDISHN=X]/g;

    for my $op (@op_arr) {
        my @len = ($op =~ /^[\d]+/g);
        my @symbol = ($op =~ /[MDISHN=X]$/g);
        push (@result,["@symbol","@len"]);
    }
    return @result;
}

sub reverse_complement {
  my $sequence = shift;
  $sequence =~ tr/CATG/GTAC/;
  $sequence = reverse($sequence);
  return $sequence;
}

# Usage example:  @seq_vars = find_variants($alt_range, $pri_range, $alt_sam, $fasta);
sub find_variants {

   my $alt_range = $_[0];
   my $pri_range = $_[1];
   my $alt_sam = $_[2];
   my $fasta = $_[3];

   # The first (and last) bases in $alt_seq and $pri_seq should align w/ each other.
   # Use the convert_alt_to_pri_coords.pl tool to get corresponding base positions.
   #print STDERR "samtools faidx $fasta $alt_range | tail -n +2 | tr -d '\n'";
   #print STDERR "samtools faidx $fasta $pri_range | tail -n +2 | tr -d '\n'";
   my $alt_seq = `samtools faidx $fasta $alt_range | tail -n +2 | tr -d '\n'`;
   my $pri_seq = `samtools faidx $fasta $pri_range | tail -n +2 | tr -d '\n'`;

   $alt_seq = uc($alt_seq);
   $pri_seq = uc($pri_seq);

   #my $foo = substr $pri_seq, 0, 10;
   #print "$foo\n";
   #exit(1);

   my ($alt_contig_name, $alt_start, $alt_end) = split /[:-]/, $alt_range;
   my ($pri_contig_name, $pri_start, $pri_end) = split /[:-]/, $pri_range;

   my @field = split "\t", $alt_sam;
   my $sam_flag = $field[1];
   my $sam_pri_contig = $field[2];
   my $priSamPOS = $field[3];
   my $cigar = $field[5];
   my @cigar_array = cigar_array($cigar);

   die "Primary sequence name mismatch in find_variants()." if ($pri_contig_name ne $sam_pri_contig);
   die "pri_start < priSamPOS detected in find_variants()." if ($pri_start < $priSamPOS);

   # If alt_seq is RC-oriented w.r.t. prim contig, then reverse complement the sequence so our
   # CIGAR is in-sync w/ alt_seq.
   if ($sam_flag & 0x10) {
      print STDERR "Got here.\n";
      $alt_seq = reverse_complement($alt_seq);
   }

   my $a1 = substr($alt_seq,0,15);  my $p1 = substr($pri_seq,0,15);
   print STDERR "alt_seq=$a1\tpri_seq=$p1\n";

   # Init our base position indeces
   my $altCurrPos = 0;               # position of end of current cigar op (1-based index into whole ALT contig)
   my $priCurrPos = $priSamPOS - 1;  # position of end of current cigar op (1-based index into whole primary contig)
   my $altPrevPos;  # position of end of previous cigar op (1-based index into whole ALT contig)
   my $priPrevPos;  # position of end of previous cigar op (1-based index into whole primary contig)

   # Init our variant hash array
   my %var_arr = ();

   # We will start adding variants to @var_arr when 
   #   [$altCurrPos >= $alt_start]
   # We stop adding variants when 
   #   [$altCurrPos >= $alt_end]


   foreach my $cig_op (@cigar_array) {

      my ($op, $len) = @$cig_op;

      given ($op) {
          when (($_ eq "M") or ($_ eq "X") or ($_ eq "=")) {
              $priCurrPos += $len;
              $altCurrPos += $len;
          }
          when ($_ eq "I") {
              $altCurrPos += $len;
          }
          when ($_ eq "D") {
              $priCurrPos += $len;
          }
          when ($_ eq "S") {
              $altCurrPos += $len;
          }
      }

      if (not defined $altPrevPos and $altCurrPos >= $alt_start) {
          # This is the 1st CIGAR op that moves us at least to $alt_start.
          # We may start finding variants in ALT w.r.t. primary contig now.
          # If we're not in an 'M' op then the input seq range for this ALT
          # contig was specified incorrectly.
          die "ALT start position must be located within M CIGAR op." if ($op ne "M");

          # Sync up $altPrevPos and $priPrevPos
          $altPrevPos = $alt_start-1;   # will be 1-based as we walk through CIGAR
          $priPrevPos = $pri_start-1;   # will be 1-based
      }

      if (defined $altPrevPos) {
          # We're in the VC region. If we're in a M|=|X op, we'll detect base mismatches.
          # If we're in a D or I op, then we'll simply append the deletion or insertion
          # to our variant list. If we're in an M op and the last base in the op is a 
          # mismatch, it can lead to a problem when the next op is D or I!

          # First check if we've progressed to or beyond $alt_end. If so then also
          # confirm that we're in an 'M' op.
          if ($altCurrPos >= $alt_end and $op ne "M") {
              die "ALT range end position must be located within M CIGAR op.";
          }

          my $altAllele;
          my $refAllele;

          #my $pri_seq_len = length($pri_seq);
          #my $alt_seq_len = length($alt_seq);
          #print "prim_len=$pri_seq_len\tpriPrevPos=$priPrevPos\tpri_start=$pri_start\n";
          #print "alt_len=$alt_seq_len\taltPrevPos=$altPrevPos\talt_start=$alt_start\n";

          # Remember, when $op eq "=", there are no variants.
          given ($op) {
              when (($_ eq "M") or ($_ eq "X")) {

                 # We walk through the op, detecting base mismatches as we go.
                 my $nbases = ($alt_end < $altCurrPos) ? $alt_end : $altCurrPos;
                 $nbases -= $altPrevPos;
                 # print "nbases=$nbases\n";
                 for (my $jj=1; $jj<=$nbases; $jj++) {

                     $refAllele = substr $pri_seq, $priPrevPos+$jj-$pri_start, 1;
                     $altAllele = substr $alt_seq, $altPrevPos+$jj-$alt_start, 1;
                     if ($altAllele ne $refAllele) {
                         # we have a base mismatch so add the SNV to our var hash
                         my $key = $priPrevPos+$jj;

                         # we should never hit the following, but just in case...
                         die "Hash key $key already exists." if ( exists $var_arr{$key} );

                         $var_arr{$key} = [$refAllele, $altAllele];
                     }
                 }
              }
              when ($_ eq "I") {
                  $refAllele = substr($pri_seq, $priPrevPos-$pri_start, 1);
                  $altAllele = $refAllele . substr($alt_seq, $altPrevPos-$alt_start+1, $len);
                  my $key = $priPrevPos;
                  if ( exists $var_arr{$key} ) {
                      my $len_ref = length($var_arr{$key}[0]);
                      my $len_alt = length($var_arr{$key}[1]);
                      my $prev_len = ($len_ref > $len_alt) ? $len_ref : $len_alt;
                      $len_ref = length($refAllele);
                      $len_alt = length($altAllele);
                      my $all_len = ($len_ref > $len_alt) ? $len_ref : $len_alt;
                      $all_len = ($all_len > $prev_len) ? $all_len : $prev_len;
                      $refAllele = substr($refAllele, 0, 1);
                      $var_arr{$key} = [$refAllele, "<*>", $all_len]; 
                  } elsif ($len > $MAX_STD_INDEL) {
                      $refAllele = substr($refAllele, 0, 1);
                      $var_arr{$key} = [$refAllele, "<INS>", $len];
                  } else {
                      $var_arr{$key} = [$refAllele, $altAllele];
                  }
              }
              when ($_ eq "D") {
                  $refAllele = substr($pri_seq, $priPrevPos-$pri_start, 1+$len);
                  $altAllele = substr($refAllele,0,1);
                  my $key = $priPrevPos;
                  if ( exists $var_arr{$key} ) {
                      my $len_ref = length($var_arr{$key}[0]);
                      my $len_alt = length($var_arr{$key}[1]);
                      my $prev_len = ($len_ref > $len_alt) ? $len_ref : $len_alt;
                      $len_ref = length($refAllele);
                      $len_alt = length($altAllele);
                      my $all_len = ($len_ref > $len_alt) ? $len_ref : $len_alt;
                      $all_len = ($all_len > $prev_len) ? $all_len : $prev_len;
                      $refAllele = substr($refAllele, 0, 1);
                      $var_arr{$key} = [$refAllele, "<*>", $all_len];
                  } elsif ($len > $MAX_STD_INDEL) {
                      $refAllele = substr($refAllele, 0, 1);
                      $var_arr{$key} = [$refAllele, "<DEL>", -$len];
                  } else {
                      $var_arr{$key} = [$refAllele, $altAllele];
                  }
              }
              when ($_ eq "S" or $_ eq "H") {
                  die "Not supposed to ever reach clipping CIGAR op when w/i our VC region.";
              }
          }
          $altPrevPos = $altCurrPos;
          $priPrevPos = $priCurrPos;
      }

      #print "altCurrPos=$altCurrPos\talt_end=$alt_end\n";
      last if ($altCurrPos >= $alt_end);
   }

   return \%var_arr;
}

my $seq1_contig_name; my $seq1_start; my $seq1_end;
my $seq2_contig_name; my $seq2_start; my $seq2_end;
($seq1_contig_name, $seq1_start, $seq1_end) = split /[:-]/, $range1;
($seq2_contig_name, $seq2_start, $seq2_end) = split /[:-]/, $range2;

if (not defined $seq1_contig_name or not defined $seq1_start or not defined $seq1_end) {
   die "SEQ1_RANGE1 arg must be in CONIG:INT1-INT2 format.";
}

die "SEQ2_RANGE contig must be specified." if (not defined $seq2_contig_name);

$seq2_start = 0 if (not defined $seq2_start);
$seq2_end = 0 if (not defined $seq2_end);

# For SEQ1_RANGE, find corresponding range on primary assembly contig
my ($seq1_pri_contig_name, $seq1_pri_contig_start, $seq1_pri_contig_end) = 
   prim_contig_range($seq1_contig_name, $seq1_start, $seq1_end);


# SEQ2_RANGE can be a fully specified range (CONTIG:INT1-INT2 format), or
# just a contig name w/o position vals. If the latter, then SEQ2_RANGE may
# be determined from the primary assembly alignment of SEQ1_RANGE. We deal
# wth SEQ2_RANGE in one of several ways when the base coords are missing.

my ($seq2_pri_contig_name, $seq2_pri_contig_start, $seq2_pri_contig_end);

if ($seq2_start < 1 or $seq2_end < 1) { 

   # range2 base coords are unspecified, e.g. = "chr3"

   if ( $seq2_contig_name eq $seq1_contig_name ) {

       # This is homozygous haplotypes case
       $seq2_start = $seq1_start;
       $seq2_end = $seq1_end;
       ($seq2_pri_contig_name, $seq2_pri_contig_start, $seq2_pri_contig_end) =
          ($seq1_pri_contig_name, $seq1_pri_contig_start, $seq1_pri_contig_end);

       $range2 = $range1;
   } elsif ( $seq2_contig_name eq $seq1_pri_contig_name ) {

       # This is the SEQ2_RANGE = "chr7" case. Given the specified seq2 contig
       # agrees with $seq1_pri_contig_name, the following is guaranteed to work.

       ($seq2_pri_contig_name, $seq2_pri_contig_start, $seq2_pri_contig_end) =
          ($seq1_pri_contig_name, $seq1_pri_contig_start, $seq1_pri_contig_end);

       $seq2_start = $seq1_pri_contig_start;
       $seq2_end = $seq1_pri_contig_end;
       $range2 = "$seq2_contig_name:$seq2_start-$seq2_end";
   } elsif ( $seq2_contig_name =~ /_alt$/ or $seq2_contig_name =~ /^HLA/ ) {

       # This is the heterozygous case when both ranges are on ALT contigs.
       # Must get $range2 start/end positions via $seq1_pri range.
       my $ret_str = `$CONVERT_PRI_TO_ALT $alt_sam $seq2_contig_name $seq1_pri_contig_start $seq1_pri_contig_end`;
       chomp $ret_str;
       my @field = split "\t", $ret_str;
       if ($field[0] ne $seq1_pri_contig_name) {
          die "SEQ1_RANGE and SEQ2_RANGE do not overlap.";
       }
       $seq2_start = $field[4];
       $seq2_end = $field[5];
       $range2 = "$seq2_contig_name:$seq2_start-$seq2_end";

         # We could just copy seq1_pri_contig range to seq2_pri_contig, but
         # let's do it the hard way to catch some unexpected errors.
       ($seq2_pri_contig_name, $seq2_pri_contig_start, $seq2_pri_contig_end) = 
          prim_contig_range($seq2_contig_name, $seq2_start, $seq2_end);

       #print STDERR "PRI_RANGE2 = $seq2_pri_contig_name:$seq2_pri_contig_start-$seq2_pri_contig_end\tRANGE2 = $range2\n";
   } else {
       die "SEQ1_RANGE and SEQ2_RANGE do not overlap."
   }

} else {
   # This is the fully-specified CONTIG:START-END case for range2.
   # Later we'll check if SEQ2_RANGE pri alignment agrees with SEQ1_RANGE's
 
   ($seq2_pri_contig_name, $seq2_pri_contig_start, $seq2_pri_contig_end) =
      prim_contig_range($seq2_contig_name, $seq2_start, $seq2_end);
}

# More error checking...
die "Specified seq ranges do not map to same primary contig." if ( $seq1_pri_contig_name ne $seq2_pri_contig_name);

if ($seq1_pri_contig_start != $seq2_pri_contig_start or $seq1_pri_contig_end != $seq2_pri_contig_end) {
   die "Input seq ranges not aligned to same region on $seq1_pri_contig_name."
}

# Now we can walk through our 2 seq's and find variants w.r.t. primary assembly.

my $pri_name = $seq1_pri_contig_name;
my $pri_start = $seq1_pri_contig_start;
my $pri_end = $seq1_pri_contig_end;
my $primary_range = "$pri_name:$pri_start-$pri_end";
#print STDERR "$range1 and $range2 align to primary assembly region $primary_range\n";

my $seq1_alt_sam = "";
if ($seq1_contig_name =~ /_alt$/ or $seq1_contig_name =~ /^HLA/) {
   open my $altFH, "$alt_sam" or die $!;
   while (my $line=<$altFH>) {
      chomp $line;
      next if ($line =~ /^@/);
      my @t = split "\t", $line;
      if ($t[0] eq $seq1_contig_name and $t[1] < 200) {
          $seq1_alt_sam = $line;
          last;
      }
   }
   close $altFH;
}

my $seq2_alt_sam = "";
if ($seq2_contig_name =~ /_alt$/ or $seq2_contig_name =~ /^HLA/) {
   open my $altFH, "$alt_sam" or die $!;
   while (my $line=<$altFH>) {
      chomp $line;
      next if ($line =~ /^@/);
      my @t = split "\t", $line;

      # for now we are only handling the primary alignments of ALT's
      if ($t[0] eq $seq2_contig_name and $t[1] < 200) {
          $seq2_alt_sam = $line;
          last;
      }
   }
   close $altFH;
}


my $seq1_vars;  # This will be a reference to a hash.
if ($seq1_alt_sam) {

   # Each element of seq1_vars is a list, e.g. [43322356, "A", "G"]
   # The first item is the primary contig position. The second
   # item is the REF allele, and the third item is the ALT allele.
   # All variants are on $seq1_pri_contig_name.

   $seq1_vars = find_variants($range1, $primary_range, $seq1_alt_sam, $fasta);
}

my @var1_keys = (sort { $a <=> $b } keys(%$seq1_vars) );
my $foo = scalar @var1_keys;
#print "length seq1_vars = $foo\n";
#for my $key (@var1_keys) {
#   print "$key  $seq1_vars->{$key}[0] $seq1_vars->{$key}[1]\n";
#}

my $seq2_vars;
if ($seq2_alt_sam) {
   $seq2_vars = find_variants($range2, $primary_range, $seq2_alt_sam, $fasta);
}

my @var2_keys = (sort { $a <=> $b } keys(%$seq2_vars) );
$foo = scalar @var2_keys;
#print "length seq2_vars = $foo\n";
#for my $key (@var2_keys) {
#   print "$key  $seq2_vars->{$key}[0] $seq2_vars->{$key}[1]\n";
#}

my @sorted_keys = uniq (sort { $a <=> $b } (keys(%$seq1_vars), keys(%$seq2_vars) ) );


# Now write out our VCF with the variants in the two haplotypes.

my $VARID = ".";
my $QUAL = ".";
my $FILTER=".";
my $INFO = ".";
my $FORMAT = "GT";

for my $poskey (@sorted_keys) {

   my $hap1_ref = "";
   my $hap2_ref = "";
   my $hap1_alt = "";
   my $hap2_alt = "";
   my $ref = "";
   my $alt = "";
   my $gtype = "";
   my $hap1_len;
   my $hap2_len;
   my $info = $INFO;
   if ( exists $seq1_vars->{$poskey} ) {
       my $sv = $seq1_vars->{$poskey}; # ref to the anon array
       ($hap1_ref, $hap1_alt, $hap1_len) = @$sv;  # this may not touch $hap1_len
   }
   if ( exists $seq2_vars->{$poskey} ) {
       my $sv = $seq2_vars->{$poskey};
       ($hap2_ref, $hap2_alt, $hap2_len) = @$sv;  # this may not touch $hap2_len
   }
   if ( $hap1_ref ne "" and $hap2_ref eq "" ) {
       $ref = $hap1_ref;
       $alt = $hap1_alt;
       $gtype = "1|0";
   } elsif ( $hap1_ref eq "" and $hap2_ref ne "" ) {
       $ref = $hap2_ref;
       $alt = $hap2_alt;
       $gtype = "0|1";
   } elsif ( $hap1_ref ne "" and $hap2_ref ne "" ) {

       if (substr($hap1_ref,0,1) ne substr($hap2_ref,0,1)) {
          print STDERR "Mismatched REF alleles detected at $pri_name:$poskey\n";
          print STDERR "    hap1_ref=$hap1_ref\thap2_ref=$hap2_ref\t\n";
          die;
       }

       if ( length($hap1_ref) != length($hap2_ref) ) {
          # The following will force consistent ref alleles
          if (length($hap1_ref) > 1) {
             $hap1_alt = "<DEL>";
             $hap1_len = -(length($hap1_ref)-1);
             $hap1_ref = substr($hap1_ref,0,1);
          }
          if (length($hap2_ref) > 1) {
             $hap2_alt = "<DEL>";
             $hap2_len = -(length($hap2_ref)-1);
             $hap2_ref = substr($hap2_ref,0,1);
          }
       }

       die "REF alleles not equal at $pri_name:$poskey" if ($hap1_ref ne $hap2_ref);

       $ref = $hap1_ref;
       $alt = $hap1_alt;
       $gtype = "1|1";
       if ($hap1_alt eq "<DEL>" and $hap2_alt eq "<DEL>") {
           $hap2_alt = "<DEL2>" if ($hap1_len != $hap2_len);
       }
       if ($hap1_alt eq "<INS>" and $hap2_alt eq "<INS>") {
           $hap2_alt = "<INS2>" if ($hap1_len != $hap2_len);
       }
       if ($hap1_alt eq "<*>" and $hap2_alt eq "<*>") {
           $hap2_alt = "<*2>" if ($hap1_len != $hap2_len);
       }
       if ( $hap1_alt ne $hap2_alt ) {
           $alt = "$hap1_alt,$hap2_alt";
           $gtype = "1|2";
       }
   }

   $info = "SVLEN=$hap1_len" if ( defined $hap1_len );
   if ( defined $hap2_len ) {
      if (not defined $hap1_len) {
         $info = "SVLEN=$hap2_len";
      } elsif ($hap1_len != $hap2_len) {
         $info .= ";SVLEN2=$hap2_len";
      }
   }
   print "$pri_name\t$poskey\t$VARID\t$ref\t$alt\t$QUAL\t$FILTER\t$info\t$FORMAT\t$gtype\n"
}
