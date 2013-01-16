#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# Gbreak: Convert a genbank file into a fasta file
# Written by Charles, Tue Jun 16 11:29:50 PDT 2009
# Modified Oct 23 by Simon Eng, to add option to extract sequences only from
# CDSes without /product qualifiers

my $INPUT;
my $result;
my %region;
my $OUTPUT;
my $AA;
my $DNA;
my $REALDNA;
my $LENGTH=60;
my $NOPRODUCT;
my $NOBLANK;

# If there are command line arguments:
if (defined $ARGV[0]) {

  # Process command line arguments
  $result=GetOptions(
    "output=s"=>\$OUTPUT,  # One output file
    "aa"=>\$AA,  # Amino acids from CDS
    "dna"=>\$DNA,  # The DNA section
    "realdna"=>\$REALDNA,  # The DNA section
    "length:i"=>\$LENGTH,  # Wrap length
    "productless+" => \$NOPRODUCT,    # Process only productless CDSes
    "noblank" => \$NOBLANK, # No empty entries will be printed
  );

  # The input files without flags, so you can use wildcards.
  $INPUT=$ARGV[0];
}

# Usage, or bad arguments:
if (!$result or !defined $OUTPUT or !defined $INPUT) {
  die("Usage: $0 [-r] [-n] [-d] [-a] [-l length] -o output.fas input.gb [input2.gb...]

  This program will convert the input genbank file(s) into a fasta file.

  -a : save the translation of each CDS (the amino acid sequences).
  -d : save the nucleotides of each locus (genes and intergenic sequences).
  -r : save the nucleotides of each CDS (the open reading frames).

  -n : no blank entries will be output

  -l length : wrap sequences at this many columns (default 60).  Use 0 for no
              wrapping.
");
}

my $sum=($AA?1:0)+($DNA?1:0)+($REALDNA?1:0);
if ($sum!=1) {
  die("Please specify one of -a, -d, or -r: what do you want to see?\n");
}

# Process all files:
open(OUT,">$OUTPUT") or die("$OUTPUT: $!");
foreach my $i (@ARGV) {dofile($i);}
close(OUT);

# Error checking:
if (-z $OUTPUT) {
  unlink($OUTPUT);
  die("No output was written to $OUTPUT, deleting empty file\n");
}

exit 0;

# Deal with each file:
sub dofile {
  my $INPUT=$_[0];
  if(!-f $INPUT || !-r _ ) { die("$INPUT isn't a readable file!");}
  open(IN,"<$INPUT") or die("$INPUT: $!");
  local $/="\n//\n";
  while(my $locus=<IN>) { dolocus($locus,$INPUT);}
  close(IN);
}

# Deal with each locus in the file:
sub dolocus {
  my ($locus,$filename)=@_;
  if ($locus =~ m/^\s*$/s) {return;}
  if ($locus !~ m/^\s*LOCUS +(\S+)\s(.*)\nFEATURES[^\n]*\n(.*)\nORIGIN[^\n]*(\n?.*)\n\/\/\n\s*$/s) { die("$filename: not a genbank file?\nlocus='$locus'"); }
  my ($locusname,$features,$seq)=($1,$3,$4);
  my @feature=split(/(^     source *|\n     gene *|\n     CDS   *)/,$features);
  my %cds;
  $seq=~s/\s+\d+\s+//gs;$seq=~s/\s+//gs;

  # User asked for un-edited DNA sequence:
  if (defined $DNA) {
    if (length($seq) or !defined $NOBLANK) { print OUT ">$locusname\n".wrapit($seq,$LENGTH)."\n"; }
    return;
  }

  # User asked for CDS sequences, translated or not:
  my $ltnum="0001";
  for (my $x=1;$x<=$#feature;$x++) { # Start at one, not zero, because 1 is 'source'.
    #if ($feature[$x] =~m/^\s*CDS/s) {%cds=(%cds,docds($feature[$x+1],$seq,$locusname,$filename,$ltnum++));next;}
    if ($feature[$x] =~m/^\s*CDS/s) {my ($aa,$bb)=docds($feature[$x+1],$seq,$locusname,$filename,$ltnum++);$cds{$aa}=$bb;next;}
  }
  print OUT map {">$_\n".wrapit($cds{$_},$LENGTH)."\n"} grep {length($cds{$_}) or !defined $NOBLANK} sort keys %cds;
}

# Deal with each CDS as requested:
sub docds {
  my ($cds,$seq,$locusname,$filename,$ltnum)=@_;

  # Unwrap lines:
  while ($cds=~s/([^"])\n +([^ \/])/$1 $2/sg) {;}

  # Generate or get the locustag:
  my $ltag="${locusname}_$ltnum";
  if ($cds=~m/\n *\/locus_tag="([^"]*)"\n/s) {$ltag=$1;}

  # if Aminoacid requested, return translation:
  if (defined $AA) {
    if ($cds!~m/\n *\/translation="([^"]*)"/s) { die("no translation found for $filename locus $locusname: $cds");}
    my $trans=$1;
    $trans=~s/\s+//g; # Unwrap code left spaces between fragments
    return ($ltag,$trans);
  }

  if (!defined $REALDNA) {die("Shouldn't get here");}

  # Find the region:
  #              (1              (2        (3         (4     (5(6           (7       
  if ($cds!~m/^ *(complement[(])?(join[(])?(<?\d+)\.\.(\d+>?)(,(<?\d+)[.][.](\d+>?))?/s) {
    die("cds doesn't have a recognizeable region: ==>$cds<==");
  }

  # Read the numbers:
  my ($comp,$joined,$start,$end,$start2,$end2)=($1,$2,$3,$4,$6,$7);

  # Find the length:
  my $slen=$end-$start+1;
  if ($slen<0) {die("cds $ltag in $filename has start after end, not supposed to happen.\n");}

  # Get the subsequence
  my $subseq=substr($seq,$start-1,$slen); # Base 1 
  if (defined $joined) {
    my $slen2=$end-$start+1;
    $subseq.=substr($seq,$start2-1,$slen2);
  }
  if (defined $comp) {$subseq=revcomplement($subseq);}

  return ($ltag,$subseq);
}

# Reverse-complement a sequence.  Doesn't do IUPAC ambiguity symbols.
sub revcomplement {
  my $in=$_[0];
  $in=reverse($in);
  $in=~y/ATCGatcg/TAGCtagc/;
  return $in;
}

sub wrapit {
  my ($in,$len)=@_;
  if (!$len) {return $in;}

  my $out="";
  while (defined $in and $in ne "") {
    $out.=substr($in,0,$len)."\n";
    #if ($len>length($in)) {$in="";} else {$in=substr($in,$len);}
    $in=substr($in,(length($in)<$len)?length($in):$len);
  }
  $out=~s/\n$//;
  return $out;
}
