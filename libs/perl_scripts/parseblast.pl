#!/usr/bin/perl -wT
# TODO: blast -m8?
# TODO: This needs to take a BLAST score ratio file
use strict;

# Written by Charles Howes Fri Jun  5 10:28:34 PDT 2009
# Updated by Simon Eng Thu Oct 15
# Updated again by Charles Mon Mar  1 09:08:52 PST 2010
# Wed Mar  9 12:16:28 PST 2011 Charles updated for blast 2.2.24+

use Getopt::Long;
use Cwd;

use Data::Dumper;

# Where am I?
my $HOME=$0;
$HOME=~s/[^\/]*$//;
if ($HOME eq "") {$HOME=cwd();}

my $result;

# Global configuration variables:
my $OUTPUT;
my $BLAST;
my $SCORE;

my $MIN_HIT_COVERAGE;
my $MIN_QUERY_COVERAGE;
my $MAX_EXPECT = 1e-6;
my $MIN_LENGTH = 20;
my $MAX_LENGTH;
my $MIN_IDENTITY;
my $MAX_IDENTITY;
my $MIN_POSITIVES;
my $MAX_GAPS;
my $BSR_FILE;
my $BLAST_SCORE_RATIO;
my $HEADER;
my $LIMIT=5;

# If there are command line arguments:
if (defined $ARGV[0]) {

  # Process command line arguments
  $result=GetOptions(
    "output=s"=>\$OUTPUT,     # Output filename
    "min-score:f"=>\$SCORE,     # Bit score cutoff
    "min-hit-coverage:f" => \$MIN_HIT_COVERAGE,     # Minimum hit coverage
    "min-query-coverage:f" => \$MIN_QUERY_COVERAGE, # Minimum query coverage
    "max-expect:f" => \$MAX_EXPECT,         # Maximum expect
    "min-length:i" => \$MIN_LENGTH,         # Minimum length
    "max-length:i" => \$MAX_LENGTH,         # Maximum length
    "min-identity:f" => \$MIN_IDENTITY,     # Minimum identity
    "max-identity:f" => \$MAX_IDENTITY,     # Maximum identity
    "min-positives:f" => \$MIN_POSITIVES,   # Minimum positives
    "max-gaps:f" => \$MAX_GAPS,             # Maximum gaps
    "limit:i" => \$LIMIT,	# How many hits to keep
    "input=s"=>\$BLAST,       # COG BLAST file
    "header"=>\$HEADER,       # Show the header?
    "bsr-file=s" => \$BSR_FILE,       # BLAST score ratio file
    "min-bsr-ratio:f" => \$BLAST_SCORE_RATIO    # Minimum BLAST score ratio threshold
  );
}

# Check for limit
if ($LIMIT < 1) {die("The -l limit needs to be greater than zero, not $LIMIT.\n");}

# Usage, or bad arguments:
if (!$result or !defined $OUTPUT or !defined $BLAST) {
  die("Usage: $0 [options] -i blast_file -o output_file 

This program will read a BLASTed file and write a file showing the summary of
all hits.

Required flags:

  -i blast_file, --input=blast_file
                        the BLASTed file
  -o output_file, --output=output_file
                        where to write the CSV-formatted output file

Options:

  --min-score=SCORE     bitscore threshold (an integer or decimal)
  --max-expect=EXPECT   maximum expect-value threshold (a decimal; default 1e-6)

  --min-length=LENGTH, --max-length=LENGTH
                        hit length thresholds (an integer; default 20 to
                        infinity)
  --min-hit-coverage    minimum hit coverage threshold (a proportion)
  --min-query-coverage  minimum query sequence coverage threshold (a
                        proportion)

  --min-identity=IDENTITY, --max-identity=IDENTITY
                        identity thresholds (proportions)
  --min-positives=POSITIVES
                        minimum positives threshold (a proportion)
  --max-gaps=GAPS       maximum gaps threshold (a proportion)

  --bsr-file=FILE       uses the reference bit score value file FILE to
                        generate BLAST score ratios (should end in .csv)
  --min-bsr-ratio=RATIO minimum BLAST score ratio threshold (a proportion)

  -l limit              limit the number of hits to 'limit' (default $LIMIT)

  -h                    show header in output file
");
}
if ($BLAST!~m/^(.+)$/i) {die("The BLAST file had bad characters in the filename\n");}
$BLAST=$1;
if (! -r $BLAST) {die("The BLAST file wasn't readable\n");}

# Check the output filename:
if ($OUTPUT!~m/^([\/a-z0-9_.-]+)$/i) {die("The output file had bad characters in the filename\n");}
$OUTPUT=$1;

# Check the score:
if (defined $SCORE) {
  if ($SCORE!~m/^([0-9]{1,3})$/) {
    die("The score had bad characters; expecting a number between 0 and 999\n");
  }
  $SCORE=($1)+0;
}

# Check other parameters

my $int_regexp = qr/^\d+$/;
my $float_regexp = qr/^\d+(\.\d+)?|\de\-?\d+$/;

if (defined $MIN_HIT_COVERAGE) {
    if ($MIN_HIT_COVERAGE !~ $float_regexp) {
        die('The hit_coverage must be a non-negative decimal number');
    }
}

if (defined $MAX_EXPECT) {
    if ($MAX_EXPECT !~ $float_regexp) {
        die('The expect must be a non-negative decimal number');
    }
}

if (defined $MIN_LENGTH) {
    if ($MIN_LENGTH !~ $int_regexp) {
        die('The minimum length must be a non-negative integer');
    }
}

if (defined $MAX_LENGTH) {
    if ($MAX_LENGTH !~ $int_regexp) {
        die('The maximum length must be a non-negative integer');
    }
}

if (defined $MIN_IDENTITY) {
    if ($MIN_IDENTITY !~ $float_regexp) {
        die('The minimum length must be a non-negative decimal number');
    }
}

if (defined $MAX_IDENTITY) {
    if ($MAX_IDENTITY !~ $float_regexp) {
        die('The maximum length must be a non-negative decimal number');
    }
}

if (defined $MIN_POSITIVES) {
    if ($MIN_POSITIVES !~ $float_regexp) {
        die('The minimum positives must be a non-negative decimal number');
    }
}

if (defined $MAX_GAPS) {
    if ($MAX_GAPS !~ $float_regexp) {
        die('The minimum gaps must be a non-negative decimal number');
    }
}

# Check that the BSR file exists; print a warning if nonexistent

my %bsrs;

if (defined $BSR_FILE) {
    open BSR, $BSR_FILE or die "Cannot open $BSR_FILE";
    while (<BSR>) {
        chomp;
        #s/ //g;
        if ($_ =~ /^(.+)\s+(\d+(\.\d+)?(e[\+\-]\d{2,})?)$/) {
            my $locus = $1;
            my $max_score = $2;
            $locus =~ s/\s+$//g;
            $locus =~ s/\s//g; #XXX
            $bsrs{$locus} = $max_score;
        }
    }
    close BSR;
}
else {
    print STDERR "Warning: no blast score ratio file specified: -b bsrfile.csv\n";
}

#############################################################################
# Reading the BLAST file:
open(IN,"<$BLAST") or die("$BLAST: $!");

# Try to recognize the file:
my $firstline=<IN>;
if ($firstline !~ m/^T?BLAST[PXN] /) { die("$0: $BLAST is not a BLAST result file\n"); }

my %letters;  # Length of each protein - always

# Main loop: Break input into queries:
open(OUT,">$OUTPUT") or die("$OUTPUT: $!");
my $query="";
while (<IN>) {   # For each line in the file:
  if (/^(T?BLAST[PXN] |Query=|  Database:)/) {processquery($query);$query="";}
  $query.=$_;
}
processquery($query);  # Catch the final sequence
close(OUT);

sub processquery {
  my $block=$_[0];

  my $dbg=0;

  # Find, cut off and dissect the query string:

  if ($block !~m/^Query/) {return;}

  if ($block=~m/^Query=\s*((\S+).*)\n +\(([,\d]+) letters\)\n/s) {return processoldblast($block);}
  if ($block=~m/^Query=\s*((\S+).*)\nLength=([,\d]+)\n/s) {return processnewblast($block);}

  die($block);
}

sub processoldblast {
  my $block=$_[0];
  if ($block!~s/^Query=\s*((\S+).*)\n +\(([\d,]+) letters\)\n//s) {die("processoldblast: failed regex:\n$block");}
  my $qb=$1;  # Query block
  my $q=$2;   # Just the first word of the query
  if (defined $letters{$q}) {warn("$q appears more than once in $BLAST\n");}
  $letters{$q}=$3;  # How many letters in this sequence
  $letters{$q}=~s/,//g;

  $qb=~s/\n */ /g;
  $qb =~ s/\s+$//;

  if ($block=~m/\* No hits found \*/s) {return;} # No-hit block

  if ($block!~s/.*Sequences producing significant alignments:[^\n]*\n\n//s) {
      print $block."\n";
      #die("No summary header?\n$block");
      return;
  }

  if ($block!~s/^((\S[^\n]+\n)+)\n//s) {die("No summary?\n$block");}
  my $sum=$1;  # The summary is not actually used

  # Split hits into an array:
  $block=~s/^>//s;
  my @hits=split(/\n>/,$block);

  # Pull down information for each hit
  my @entries;
  for my $hit (grep {defined} @hits[0..$LIMIT-1]) {

  if ($hit !~ m/
    ^(\S+)(.*?)                                     # 1=hit name, 2=product
    \n
    [ ]{10}Length[ ]=[ ]([\d,]+)                       # 3=length
    \n
    \n
    [ ]*Score[ ]*=[ ]*([\d.e+]+)[ ]*bits[ ]*\((\d+)\),    # 4=bitscore, 5=raw bitscore (unused)
    [ ]*Expect(?:[()\d]*)[ ]*=[ ]*([\d.e+-]+)                      # 6=evalue
    (?:,[ ]*Method:\s*(?:[^\n]+))?                      # optional method section
    \n
    [ ]*Identities[ ]*=[ ]*(\d+)\/(\d+)[ ]*\((\d+)%\)  # 7=identities, 8=matchlength, 9=idpercent
    (?:,[ ]*Positives[ ]*=[ ]*(\d+)\/(\d+)[ ]*\((\d+)%\))?    # 10=positives, 11=matchlength, 12=pospercent
    ([ ]*Gaps[ ]*=[ ]*(\d+)\/(\d+)[ ]*\((\d+)%\))?      # 12=gaps, 13=matchlength, 14=gappercent
    /sx
) {
    return;  
   # die("Hit doesn't match the pattern:\n===Hit:\n$hit fail");

}

   
    my %entry = (
      name => $1,
      product => $1.$2,
      query_length => $letters{$q},
      length => $3,
      bitscore => $4,
      expect => $6,
      identities => $7,
      identities_pct => $9,
      hit_length => $8,
      positives => $10||0,
      positives_pct => $12||0,
      gaps => 0,
      gaps_pct => 0.0,
      max_bitscore => undef,
      bsr => undef,
      hit_coverage => undef,
      query_coverage => undef,
      ec_numbers => [],
      kegg_ids => []
    );
  


    $entry{product} =~ s/\s+/ /g; # Turn whitespaces into space
    $entry{product} =~ s/>[^#]*#/#/; # delete multiple names prior to any hash tags
    $entry{product} =~ s/>[^#]*$//; # delete multiple names if no hash tags
    $entry{product} =~ s/^ //; # delete leading space
    $entry{product} =~ s/ $//; # delete trailing space
    $entry{product} =~ s/"//g; # delete quotes

    $entry{expect} =~ s/^e/1e/; # Fix malformed e-values (e-107 => 1e-107)

    
    if (defined $BSR_FILE) {
      if (!defined $bsrs{$q}) { 
         print("$q does not have a bsr value\n"); 
         next;
      }
      $entry{max_bitscore} = $bsrs{$q};
      $entry{bsr} = $entry{bitscore} / $entry{max_bitscore};
    }

    push @entries, \%entry;
  }

  # Filter by bit scores:
  if (defined $SCORE) {
    if (grep {!defined $_->{bitscore}} @entries) { die("Some hits don't have bit scores in query $q"); }
    @entries = grep {$_->{bitscore} >= $SCORE} @entries;
  }

  # BSRs
  if (defined $BLAST_SCORE_RATIO) {
    if (my @bad=map {$_->{name}} grep {!defined $_->{bsr}} @entries) { die("Some hits don't have bsr values in query $q (@bad)"); }
    @entries = grep {$_->{bsr} >= $BLAST_SCORE_RATIO} @entries;
  }

  # Expect
  if (defined $MAX_EXPECT) {
    if (grep {!defined $_->{expect}} @entries) { die("Some hits don't have expect values in query $q"); }
    @entries = grep {$_->{expect} <= $MAX_EXPECT} @entries;
  }


  # Length
  if (defined $MIN_LENGTH) {
    if (grep {!defined $_->{hit_length}} @entries) { die("Some hits don't have lengths for query $q"); }
    @entries = grep {$_->{hit_length} >= $MIN_LENGTH} @entries;
  }

  if (defined $MAX_LENGTH) {
    if (grep {!defined $_->{hit_length}} @entries) { die("Some hits don't have lengths for query $q"); }
    @entries = grep {$_->{hit_length} <= $MAX_LENGTH} @entries;
  }

  # Positives
  if (defined $MIN_POSITIVES) {
    if (grep {!defined $_->{positives_pct}} @entries) { die("Some hits don't have positives for query $q"); }
    @entries = grep {$_->{positives_pct} >= $MIN_POSITIVES} @entries;
  }

  # Identities
  if (defined $MIN_IDENTITY) {
    if (grep {!defined $_->{identities_pct}} @entries) { die("Some hits don't have identities for query $q"); }
    @entries = grep {$_->{identities_pct} >= $MIN_IDENTITY} @entries;
  }

  if (defined $MAX_IDENTITY) {
    if (grep {!defined $_->{identities_pct}} @entries) { die("Some hits don't have identities for query $q"); }
    @entries = grep {$_->{identities_pct} <= $MAX_IDENTITY} @entries;
  }

  # Gaps
  if (defined $MAX_GAPS) {
    if (grep {!defined $_->{gaps_pct}} @entries) { die("Some hits don't have gaps for query $q"); }
    @entries = grep {$_->{gaps_pct} >= $MAX_GAPS} @entries;
  }

  # For each entry, calculate percentage hit and query coverage
  for (@entries) {
    $_->{hit_coverage} = $_->{length} ? ($_->{hit_length} - $_->{gaps}) / $_->{length} : undef;
    $_->{query_coverage} = $_->{length} ? ($_->{hit_length} - $_->{gaps}) / $_->{query_length} : undef;
  }
  if (defined $MIN_HIT_COVERAGE) {
    @entries = grep {defined $_->{hit_coverage} and $_->{hit_coverage} >= $MIN_HIT_COVERAGE} @entries;
  }
  if (defined $MIN_QUERY_COVERAGE) {
    @entries = grep {defined $_->{query_coverage} and $_->{query_coverage} >= $MIN_QUERY_COVERAGE} @entries;
  }

  if (!@entries) { return; } # No hits survived

  # Extract some other things
  for my $entry (@entries) {
    # Deal with EC numbers
    for ($entry->{product} =~ /\[EC:(.+?)\]/g) {
        my $ec_string = $1;
        for ($ec_string =~ /(((\d+|\-)\.){3}(\d+|\-))/g) {
            push @{$entry->{ec_numbers}}, $1;
        }
    }
    my %seen;
    my @ec_numbers = grep {not $seen{$_}++} @{$entry->{ec_numbers}};
    $entry->{ec_numbers} = \@ec_numbers;

    # Deal with KEGG identifiers
    for ($entry->{product} =~ /([a-z]{3}:\S+)/g) {
        push @{$entry->{kegg_ids}}, $1;
    }

    my %kegg_seen = map {($_||"")=>1} @{$entry->{kegg_ids}};;
    my @kegg_ids = keys %kegg_seen;
    $entry->{kegg_ids} = \@kegg_ids;

    # Death to whitespace:
    $entry->{product} =~ s/\s+/ /g;
    $entry->{product} =~ s/^\s*//;
    $entry->{product} =~ s/\s*$//;
    $entry->{product} =~ s/\s+;/;/g;

    # Separate products with semicolons
    $entry->{product} =~ s/\s([a-z]+\|[A-Za-z0-9_]+\.\d+\|)/; $1/g;
  }

  my $sep = "\t";
  if (defined $HEADER) {
    print OUT "# query";
    print OUT $sep;
    print OUT "q_length";
    print OUT $sep;
    print OUT "product";
    print OUT $sep;
    print OUT "bitscore";
    print OUT $sep;
    print OUT "max_bitscore";
    print OUT $sep;
    print OUT "bsr";
    print OUT $sep;
    print OUT "expect";
    print OUT $sep;
    print OUT "hit_length";
    print OUT $sep;
    print OUT "length";
    print OUT $sep;
    print OUT "identities";
    print OUT $sep;
    print OUT "identities_pct";
    print OUT $sep;
    print OUT "positives";
    print OUT $sep;
    print OUT "positives_pct";
    print OUT $sep;
    print OUT "gaps";
    print OUT $sep;
    print OUT "gaps_pct";
    print OUT $sep;
    print OUT "hit_coverage";
    print OUT $sep;
    print OUT "query_coverage";
    print OUT $sep;
    print OUT "ec_numbers";
    print OUT $sep;
    print OUT "kegg_ids";
    print OUT "\n";
  }

  for (@entries) {
    print OUT "\"$q\"";
    print OUT $sep;
    print OUT $_->{query_length};
    print OUT $sep;
    print OUT "\"$_->{product}\"";
    print OUT $sep;
    print OUT $_->{bitscore};
    print OUT $sep;
    print OUT (defined $_->{max_bitscore} ? $_->{max_bitscore} : 'null');
    print OUT $sep;
    print OUT (defined $_->{bsr} ? $_->{bsr} : 'null');
    print OUT $sep;
    print OUT (defined $_->{expect} ? $_->{expect} : 'null');
    print OUT $sep;
    print OUT $_->{hit_length};
    print OUT $sep;
    print OUT $_->{length};
    print OUT $sep;
    print OUT $_->{identities};
    print OUT $sep;
    print OUT $_->{identities_pct};
    print OUT $sep;
    print OUT $_->{positives};
    print OUT $sep;
    print OUT $_->{positives_pct};
    print OUT $sep;
    print OUT (defined $_->{gaps} ? $_->{gaps} : 'null');
    print OUT $sep;
    print OUT (defined $_->{gaps_pct} ? $_->{gaps_pct} : 'null');
    print OUT $sep;
    print OUT $_->{hit_coverage};
    print OUT $sep;
    print OUT $_->{query_coverage};
    print OUT $sep;
    print OUT "\"" . (join ', ', @{$_->{ec_numbers}}) . "\"";
    print OUT $sep;
    print OUT "\"" . (join ', ', @{$_->{kegg_ids}}) . "\"";
    print OUT "\n";
  }
}

sub processnewblast {
  my $block=$_[0];
  if ($block!~s/Query=\s*(\S+)(.*?)\nLength=([\d,]+)\n//s) {die("processnewblast: failed regex:\n$block");}
  my $qb=$1.$2;  # Query block
  my $q=$1;   # Just the first word of the query
  if (defined $letters{$q}) {warn("$q appears more than once in $BLAST\n");}
  $letters{$q}=$3;  # How many letters in this sequence
  $letters{$q}=~s/,//g;

  $qb =~ s/\s+/ /gs;
  $qb =~ s/\s$//;
  $qb =~ s/^\s//;

  if ($block=~m/\* No hits found \*/s) {return;} # No-hit block


  if ($block!~s/^ {70}Score\s+E\n//s) {die("processnewblast: failed regex2:\n$block");}
  if ($block!~s/^Sequences producing significant alignments:\s+\(Bits\)\s+Value\n\n//s) {die("processnewblast: failed regex3:\n$block");}
  if ($block!~s/^( *\S[^\n]*\n)*\n\n//s) {die("processnewblast: failed regex4:\n$block");}

  # Split hits into an array:
  $block=~s/^>//s;
  my @hits=split(/\n\n>/,$block);

  # Pull down information for each hit, up to the limit:
  my @entries;
  for my $hit (grep {defined} @hits[0..$LIMIT-1]) {

    if ($hit !~ s/
      ^\s*(\S+)(.*?)\n    # 1=hit name, 2=product
      Length=(\d+)\ *\n      # 3=length
      \ *\n
      //sx) {die("Pattern fail 1:\n$hit\n ");}

    my %entry = (
      name => $1,
      product => $1.$2,
      query_length => $letters{$q},
      length => $3,

      max_bitscore => undef,
      bsr => undef,
      hit_coverage => undef,
      query_coverage => undef,
      ec_numbers => [],
      kegg_ids => [],
    );

    # This could be a while loop:
    if ($hit !~ m/
      [ ]Score[ ]=[ ]*([\d.e+-]+)[ ]bits[ ]\((\d+)\),    # 1=bitscore, 2=raw bitscore (unused)
      [ ][ ]Expect(?:\(\d+\))?[ ]=[ ]*([\d.e+-]+),?                  # 3=evalue
      (?:[ ]Method:\s*([^\n]+))?\n                          # 4=method

      [ ]Identities[ ]=[ ](\d+)\/(\d+)[ ]\((\d+)%\), # 5=identities, 6=matchlength, 7=idpercent
      (?:[ ]Positives[ ]=[ ](\d+)\/(\d+)[ ]\((\d+)%\),)?  # 8=positives, 9=matchlength, 10=pospercent
      [ ]Gaps[ ]=[ ](\d+)\/(\d+)[ ]\((\d+)%\)\n      # 11=gaps, 12=matchlength, 13=gappercent
      /sx) {die("Pattern fail 2:\n$hit\n ");}

# Score =  536 bits (594),  Expect = 5e-149
# Identities = 466/572 (81%), Gaps = 6/572 (1%)
# Strand=Plus/Plus


    %entry=(%entry,
      bitscore => $1,
      expect => $3,
      identities => $5,
      identities_pct => $7,
      hit_length => $6,
      positives => $8,
      positives_pct => $10,
      gaps => $11,
      gaps_pct => $13,
    );

    $entry{product} =~ s/\s+/ /g;
    $entry{product} =~ s/^ //;
    $entry{product} =~ s/ $//;
    $entry{product} =~ s/"//g;

    if (defined $BSR_FILE && defined $bsrs{$q}) {
      $entry{max_bitscore} = $bsrs{$q};
      $entry{bsr} = $entry{bitscore} / $entry{max_bitscore};
    }

    push @entries, \%entry;
  }


  # Filter by bit scores:
  if (defined $SCORE) {
    if (grep {!defined $_->{bitscore}} @entries) { die("Some hits don't have bit scores in query $q"); }
    @entries = grep {$_->{bitscore} >= $SCORE} @entries;
  }

  # BSRs
  if (defined $BLAST_SCORE_RATIO) {
    if (grep {!defined $_->{bsr}} @entries) { die("Some hits don't have bsr values in query $q"); }
    @entries = grep {$_->{bsr} >= $BLAST_SCORE_RATIO} @entries;
  }

  # Expect
  if (defined $MAX_EXPECT) {
    if (grep {!defined $_->{expect}} @entries) { die("Some hits don't have expect values in query $q"); }
    @entries = grep {$_->{expect} <= $MAX_EXPECT} @entries;
  }

  # Length
  if (defined $MIN_LENGTH) {
    if (grep {!defined $_->{hit_length}} @entries) { die("Some hits don't have lengths for query $q"); }
    @entries = grep {$_->{hit_length} >= $MIN_LENGTH} @entries;
  }

  if (defined $MAX_LENGTH) {
    if (grep {!defined $_->{hit_length}} @entries) { die("Some hits don't have lengths for query $q"); }
    @entries = grep {$_->{hit_length} <= $MAX_LENGTH} @entries;
  }

  # Positives
  if (defined $MIN_POSITIVES) {
    if (grep {!defined $_->{positives_pct}} @entries) { die("Some hits don't have positives for query $q"); }
    @entries = grep {$_->{positives_pct} >= $MIN_POSITIVES} @entries;
  }

  # Identities
  if (defined $MIN_IDENTITY) {
    if (grep {!defined $_->{identities_pct}} @entries) { die("Some hits don't have identities for query $q"); }
    @entries = grep {$_->{identities_pct} >= $MIN_IDENTITY} @entries;
  }

  if (defined $MAX_IDENTITY) {
    if (grep {!defined $_->{identities_pct}} @entries) { die("Some hits don't have identities for query $q"); }
    @entries = grep {$_->{identities_pct} <= $MAX_IDENTITY} @entries;
  }

  # Gaps
  if (defined $MAX_GAPS) {
    if (grep {!defined $_->{gaps_pct}} @entries) { die("Some hits don't have gaps for query $q"); }
    @entries = grep {$_->{gaps_pct} >= $MAX_GAPS} @entries;
  }

  # For each entry, calculate percentage hit and query coverage
  for (@entries) {
    $_->{hit_coverage} = $_->{length} ? ($_->{hit_length} - $_->{gaps}) / $_->{length} : undef;
    $_->{query_coverage} = $_->{length} ? ($_->{hit_length} - $_->{gaps}) / $_->{query_length} : undef;
  }
  if (defined $MIN_HIT_COVERAGE) {
    @entries = grep {defined $_->{hit_coverage} and $_->{hit_coverage} >= $MIN_HIT_COVERAGE} @entries;
  }
  if (defined $MIN_QUERY_COVERAGE) {
    @entries = grep {defined $_->{query_coverage} and $_->{query_coverage} >= $MIN_QUERY_COVERAGE} @entries;
  }

  if (!@entries) { return; } # No hits survived

  # Extract some other things
  for my $entry (@entries) {
    # Deal with EC numbers
    for ($entry->{product} =~ /\[EC:(.+?)\]/g) {
        my $ec_string = $1;
        for ($ec_string =~ /(((\d+|\-)\.){3}(\d+|\-))/g) {
            push @{$entry->{ec_numbers}}, $1;
        }
    }
    my %seen;
    my @ec_numbers = grep {not $seen{$_}++} @{$entry->{ec_numbers}};
    $entry->{ec_numbers} = \@ec_numbers;

    # Deal with KEGG identifiers
    for ($entry->{product} =~ /([a-z]{3}:\S+)/g) {
        push @{$entry->{kegg_ids}}, $1;
    }

    my %kegg_seen = map {($_||"")=>1} @{$entry->{kegg_ids}};;
    my @kegg_ids = keys %kegg_seen;
    $entry->{kegg_ids} = \@kegg_ids;

    # Death to whitespace:
    $entry->{product} =~ s/\s+/ /g;
    $entry->{product} =~ s/^\s*//;
    $entry->{product} =~ s/\s*$//;
    $entry->{product} =~ s/\s+;/;/g;

    # Separate products with semicolons
    $entry->{product} =~ s/\s([a-z]+\|[A-Za-z0-9_]+\.\d+\|)/; $1/g;
  }

  my $sep = "\t";
  if (defined $HEADER) {
    print OUT "# query";
    print OUT $sep;
    print OUT "q_length";
    print OUT $sep;
    print OUT "product";
    print OUT $sep;
    print OUT "bitscore";
    print OUT $sep;
    print OUT "max_bitscore";
    print OUT $sep;
    print OUT "bsr";
    print OUT $sep;
    print OUT "expect";
    print OUT $sep;
    print OUT "hit_length";
    print OUT $sep;
    print OUT "length";
    print OUT $sep;
    print OUT "identities";
    print OUT $sep;
    print OUT "identities_pct";
    print OUT $sep;
    print OUT "positives";
    print OUT $sep;
    print OUT "positives_pct";
    print OUT $sep;
    print OUT "gaps";
    print OUT $sep;
    print OUT "gaps_pct";
    print OUT $sep;
    print OUT "hit_coverage";
    print OUT $sep;
    print OUT "query_coverage";
    print OUT $sep;
    print OUT "ec_numbers";
    print OUT $sep;
    print OUT "kegg_ids";
    print OUT "\n";
  }

  for (@entries) {
    print OUT "\"$q\"";
    print OUT $sep;
    print OUT $_->{query_length};
    print OUT $sep;
    print OUT "\"$_->{product}\"";
    print OUT $sep;
    print OUT $_->{bitscore};
    print OUT $sep;
    print OUT (defined $_->{max_bitscore} ? $_->{max_bitscore} : 'null');
    print OUT $sep;
    print OUT (defined $_->{bsr} ? $_->{bsr} : 'null');
    print OUT $sep;
    print OUT (defined $_->{expect} ? $_->{expect} : 'null');
    print OUT $sep;
    print OUT $_->{hit_length};
    print OUT $sep;
    print OUT $_->{length};
    print OUT $sep;
    print OUT $_->{identities};
    print OUT $sep;
    print OUT $_->{identities_pct};
    print OUT $sep;
    print OUT $_->{positives}||'null';
    print OUT $sep;
    print OUT $_->{positives_pct}||'null';
    print OUT $sep;
    print OUT (defined $_->{gaps} ? $_->{gaps} : 'null');
    print OUT $sep;
    print OUT (defined $_->{gaps_pct} ? $_->{gaps_pct} : 'null');
    print OUT $sep;
    print OUT $_->{hit_coverage};
    print OUT $sep;
    print OUT $_->{query_coverage};
    print OUT $sep;
    print OUT "\"" . (join ', ', @{$_->{ec_numbers}}) . "\"";
    print OUT $sep;
    print OUT "\"" . (join ', ', @{$_->{kegg_ids}}) . "\"";
    print OUT "\n";


  }

}
