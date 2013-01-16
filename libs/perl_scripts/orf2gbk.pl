#!/usr/bin/perl -w
#Yaojie 2010-OCT-20;
# previous version: cghtogbk.pl
# modified features: can read mutiple bamgOUT and fasta files at the same time. 
# 		    using a genbank.yaml file as an input
#Yaojie 2010-NOV-02;
# previous version: orf_to_gbk.pl
# modified features: 1) will also generate ganbank files for non-ORF sequences.
#                    2) change name to "orf2gbk.pl" which is ...kind of simple to type...P
#
#Yaojie 2010-NOV-18;
# previous version: orf2gbk.pl
# modified features: 1) modified to accept converted format from prodigal;
#                    2) organism part to accept temparory organism information
#
###########################################################################################
use strict;
use warnings;
use lib './internal-tools/perllib';
use YAML qw(LoadFile);
use POSIX qw(strftime);
use Text::Wrap;
use Getopt::Long;
use Data::Dumper;

sub get_parameters;
#######
my @fasta_in;
my @bamgout_in;
my $bamgout_file;
my $file_id = 0;
my $locus_tag;
my @gene_tags;
my $num_genes = 0;
my $len_seq;
my $errout;
my $command_mark;
my $fasta_file;
my $gcode; #this parameter can be got from parameter file.
my $para_file;
my @output_files;
my $sample_name;

########implement for orphelia files
my @orphelia_in;

########implement for metagene files
my @metagene_in;

my @months = ('January','February','March','April', 'May', 'June', 'July', 'August', 'September', 'Octorber', 'November', 'December');


if (defined $ARGV[0]) {
 	$command_mark = GetOptions(
		"errors=s"   => \$errout,	# Error file
		"sample-name=s" => \$sample_name,	# sample name
		"fasta=s" => \@fasta_in,	# fasta files in read
		"bamg=s"  => \@bamgout_in,	# bamgOUT files in read
		"orphelia=s" => \@orphelia_in, 
		"metagene=s" => \@metagene_in,
		"out=s"   => \@output_files,   # If there is no output parameter defined the gene names would be used as the file name;	
	);
}

if ( @bamgout_in and scalar @fasta_in != scalar @bamgout_in) { die ("ERROR: Number of fasta files different with Number of Bamg files\n\tNotice: bamg files and fasta files had better be in corresponding order.\n");}

if ( @metagene_in and scalar @fasta_in != scalar @metagene_in) { die ("ERROR: Number of fasta files different with Number of MetaGene files\n\tNotice: Metagene files and fasta files had better be in corresponding order.\n");}

if ( @orphelia_in and scalar @fasta_in != scalar @orphelia_in) { die ("ERROR: Number of fasta files different with Number of Orphelia files\n\tNotice: Orphelia files and fasta files should be in corresponding order.\n");}


if (!@fasta_in or  !@fasta_in  or (!@metagene_in and !@bamgout_in and !@orphelia_in)) {
	die ("
Usage: perl $0 --bamg/metagene/orphelia file1.predict.orf
    file2.predict.orf [...] --fasta file1.fasta file2.fasta [...]
    --parameter-file parameter.file --out(optional) dna.1.gbk dna.2.gbk
    dna.3.gbk [...]

Notice:
1)  the fasta files and the software predicted files are better in
    corresponding order to make sure the number of the files are the same.
2)  this script can accept three software outputs: FgenesB, MetaGene and
    Orphelia. You can just use --bamg / --orphelia / --metagene to specify
    which program you are using.
3)  the number of output files is better to be equal to the number of genes.

"
	);
}

###parameter file checking part####

my $para_info = get_parameters($sample_name);
#print Dumper($para_info);

#extract usefule information in that file:
#	organism
#	version
#	Translation table

my @default_paras = ('Organism','Strain','Substrain','Keywords','Definition','Taxonomy','Accession','Version','Chromosome','Database link','GI number','Translation table','Codon start','References','Comments');

my @missing_parameters;
push (@missing_parameters, 'Organism'  ) unless (defined $para_info->{'Organism'});
push (@missing_parameters, 'Version'   ) unless (defined $para_info->{'Version'});
push (@missing_parameters, 'Taxonomy'  ) unless (defined $para_info->{'Taxonomy'});
push (@missing_parameters, 'Definition') unless (defined $para_info->{'Definition'});
push (@missing_parameters, 'Accession' ) unless (defined $para_info->{'Accession'});

if (@missing_parameters) {
	die ("ERROR: required parameters @missing_parameters are missing in parameter file $para_file.\n");
}

my $definition = $para_info->{'Definition'};
($definition = '.') unless ($definition);
my $organism = $para_info->{'Organism'};
my @organism_part = split (/\s+/,$organism);

if ($#organism_part >= 2) {

	$organism = $organism_part[0]." ".$organism_part[1];

} else {

	$organism = $organism;

}

#test
#foreach (@organism_part)
#{
#
#	print "==$_==\n";
#}


my $strain;
my $substrain;
if (!defined $para_info->{'Substrain'} and !defined $para_info->{'Strain'} and $organism_part[3] =~ /substr\./ ) {
	$substrain = $organism_part[4]; 
	$strain    = $organism_part[2];
} else {
	$substrain = (defined $para_info->{'Substrain'})?($para_info->{'Substrain'}):('');
	$strain    = (defined $para_info->{'Strain'})?($para_info->{'Strain'}):('');
}
die ("Error in parameter file $para_file: parameter 'Strain' required, or at least contained in parameter 'Organism' using mark 'substr.'\n") unless ($strain);
#print "Warning: No information 'Substain' contained in the parameter file." unless ($substrain);
my $locus_tag_base = uc(substr($organism_part[0], 0, 1)).lc(substr($organism_part[1], 0, 3)).$strain;

   #$locus_tag_base =~ s/\W//;

$locus_tag_base = $sample_name;

$gcode = $para_info->{'Translation table'} if (defined $para_info->{'Translation table'} and $para_info->{'Translation table'}=~/^\d+$/);
my $version;
if (defined $para_info->{'GI number'}) {
	$version = $para_info->{'Version'}."    GI:".$para_info->{'GI number'};
}else {
	$version = $para_info->{'Version'};
}
my $dblink    = (defined $para_info->{'Database link'})?$para_info->{'Database link'}:'.';
my $keywords  = (defined $para_info->{'Keywords'})?$para_info->{'Keywords'}:'.';
my $accession = $para_info->{'Accession'};

###other parameter checking part####
my @unreadable_files = grep {!-r} @fasta_in;
push @unreadable_files, grep {!-r} @bamgout_in;
die "ERROR: Files unreadable:\n\t@unreadable_files.\n" if (@unreadable_files);
####Load up the codon table, copied from cghtogbk.pl###

my $table=join("",grep {!m/^--/} <DATA>);
$table=~s/.*::= {\n {\n//s;
$table=~s/\n }\n}$//;
my @gcode=grep {/ id $gcode ,/} split(/ },\n {\n/,$table);
if (!defined $gcode[0]) {die("Genetic code $gcode not found in NCBI table\n");}
my ($nuc,$stn,$b1,$b2,$b3)=$gcode[0]=~m/ncbieaa  "([^"]{64})",.*sncbieaa "([^"]{64})".*-- Base1  ([ACTG]{64}).*-- Base2  ([ACTG]{64}).*-- Base3  ([ACTG]{64})/s;
my %codons=map {substr($b1,$_,1).substr($b2,$_,1).substr($b3,$_,1)=>substr($nuc,$_,1)} 0..63;   # Codon table

#print Dumper(\%codons);

#my %stcodon=map {substr($b1,$_,1).substr($b2,$_,1).substr($b3,$_,1)=>substr($stn,$_,1)} 0..63;  # Start codons (not needed in this program)
# Free up space:
$table=undef;
@gcode=();
####codon table loaded, copied from cghtogbk.pl###

my %dna;
my %fasta_file_dna;

foreach $fasta_file (@fasta_in) {

	open FASTA, "< $fasta_file" or die ("ERROR: Failed to open file $fasta_file : $!\n");
	#print "FASTA file inread: $fasta_file.\n";
	my $origin = $/;
	$/ = '>';
	my @fasta_seq = <FASTA>;
	$/ = $origin;
	shift @fasta_seq;
	chop @fasta_seq;
	close FASTA;
	
	foreach (@fasta_seq) {
		if ($_ =~ /^([^\n]+?)\n(.*)$/si) {
			my $dna_name;
			my $header = $1; my $seq = $2;
			if ($1 =~ /(\S+)\s*\S*/) {$dna_name = $1;}
			$dna{$dna_name} = [$header,$seq];
			$dna{$dna_name}->[1] =~ s/\s+//sig;
			push @{$fasta_file_dna{$fasta_file}}, $dna_name; 

		} else {
			die "ERROR: $_ (File: $fasta_file) cannot map to standard format.\n";
		}
	}
}

#test
#foreach (sort keys %dna) {
#	print $_,"-->\n\t",length $dna{$_}->[1],"\n";
#
#}
#
#foreach (sort keys %fasta_file_dna) {
#
#	print $_,"->",@{$fasta_file_dna{$_}},"\n";
#}


#/test

my %gene_list;
my %aa;
my %extra_info;
my $tag = 0;

if ( @metagene_in or @orphelia_in){
	my @files = (@metagene_in)?@metagene_in:@orphelia_in;
	foreach my $file (@files){
		print STDERR "FILE $file processing..\n";
		$file_id ++;
		my @contents;
		open ORP, "< $file" or die ("ERROR: Failed to open file $file : $!\n");
		@contents = <ORP>;	
		chop @contents;
		close ORP;
	
		my $block;
		my $gene_name = '';
		my $start_end;
		my $line_content;
		my $dna_name;
		my $line_mark = 0;	
		my $start;
		my $end;
		my $nt; 
		my $count_gene;
        my $gene_number = 0;
        my %seen_dna_name;

		foreach $line_content (@contents){
			$line_mark ++;
			if ($line_content =~ /^>(\d+),(\d+)_(\d+)_(\d+)_([+-])_([123\?])_([CI])_(\S+)/si) { #Orphelia contents match
				#print "File format matched to Orphelia.\n";
				next if ($7 =~ /I/);
				$tag ++;
				$dna_name = $8;
                if( defined($seen_dna_name{$dna_name}) ) {
                   $gene_number++;
                }
                else{
                   $gene_number=0;
                   $seen_dna_name{$dna_name} = 1;
                }
				$extra_info{$gene_name} = '';
				die ("ERROR: Cannot find DNA sequence $dna_name## ($file) in corresponding fasta file.\n") unless (defined $dna{$dna_name});
				$gene_name = $8."_".$2;
                #print $gene_name."\n";
				push @{$gene_list{$dna_name}}, $gene_name;
				$start = ($5 eq '+')?$3:$4;
				$end = ($5 eq '+')?$4:$3;
				$start_end = ($5 eq '+')?"$3..$4":"$4..$3";
			} elsif ($line_content =~ /#\s+gc = .*$/) {
				$count_gene = 0;
				if ($contents[$line_mark-2] =~ /^#\s+(.+?)\s*$/){
					#print "DNA is $1\n";
					$dna_name = $1;
					die ("ERROR: Cannot find DNA sequence $dna_name## ($file) in corresponding fasta file.\n") unless (defined $dna{$dna_name});
				}else {
					die ("ERROR: $file contents format error in $line_mark..\n");
				}
				next;

			} elsif ($line_content =~ /(\d+)\s+(\d+)\s+([+-])\s+[012]\s+[\d\.]+\s+(.*)$/){	#MetaGene Contents Match
				#print " match here:::$line_content \n";
				$count_gene ++;
				next if $4 =~ /partial/;
				$tag ++;
				$extra_info{$gene_name} = '';
				$gene_name = $dna_name."_".$count_gene;
				push @{$gene_list{$dna_name}}, $gene_name;
				$start = ($3 eq '+')?$1:$2;
				$end = ($3 eq '+')?$2:$1;
				$start_end = $start."..".$end;
			} else {
				next;
			}



			if ($start < $end) { 
				#print "$dna_name + $start_end with start $start end $end\n";
				$nt = substr($dna{$dna_name}->[1], $start-1, $end-$start+1);
				#print "$nt\n";
			} else {
				#print "$dna_name - $start_end with start $start end $end\n";
				$nt = revcomp (substr ($dna{$dna_name}->[1], $end-1, $start-$end+1));
				#print "$nt\n";
			}
			my $translated_aa = translate($nt);
			$translated_aa =~ s/[*]$//;
			$aa{$gene_name} = [$translated_aa,$start_end,$tag,$line_mark,$file,$dna_name, $gene_number];
			$extra_info{$gene_name} = '';	
			
		}
	}
}

if (@bamgout_in){
    foreach $bamgout_file (@bamgout_in){

	my @bamg_contents;
	#print "bamgout file in read: ",$bamgout_file,"\n";
	open BAMG, "< $bamgout_file" or die ("ERROR: Failed to open file $bamgout_file : $!\n");
	@bamg_contents = <BAMG>;	
	chomp @bamg_contents;
	close BAMG;
	
	my $block;
	my $gene_name = '';
	my $start_end;
	my $seq_header;
	my $line_content;
	my $dna_name;
	my $line_mark = 0;	
    my %seen_seq_header;
    my $orf_Number = 0;
	foreach $line_content (@bamg_contents){
		$line_mark ++;
		if ($line_content =~ /^\s+Seq name:\s+(\S+)\s+[^\n]*/) {
			$dna_name = $1;
			die ("ERROR: Cannot find DNA sequence $dna_name ($bamgout_file) in corresponding fasta file.\n") unless (defined $dna{$dna_name});
		} 
        elsif ($line_content =~ m/^>/) {

			$block = '';
			$line_content =~ s/^>(.+?)\s+GENE\s+/>$1_/s;			

			if ($line_content =~ m/^>((.+?)_(\d+)\s+(\d+)\s+-\s+(\d+)\s.*chain\s+([+-]))/si) {

				die ("DNA name $2 different with Seq name $dna_name in file $bamgout_file line $line_mark.\n") unless ($2 eq $dna_name);
				$start_end = ($6 eq '+')?"$4..$5":"$5..$4";
				$gene_name = $2."_".$3;
				$seq_header = $1;
                if( defined($seen_seq_header{$seq_header}) ) {
                   $orf_Number++;
                }
                else{
                   $orf_Number=0;
                   $seen_seq_header{$seq_header} = 1;
                }
				$tag ++;
				$extra_info{$gene_name} = '';

				if ($seq_header =~ m/ (COG\d{4} )/) { $extra_info{$gene_name} .= $1;}
				if ($seq_header =~ m/ Function: ([^#\n]+)\s/) {	$extra_info{$gene_name} .= $1;}

				push @{$gene_list{$dna_name}}, $gene_name;

			} else {
				die "ERROR: Substituted seq $line_content (line $line_mark) cannot map to standard format!\n";
			}

		} 
        elsif ($line_content =~ m/^(\S+.*)/ and $line_content !~ m/Predicted protein\(s\)/i) {
			$block .= $1;
			$aa{$gene_name} = [$block, $start_end, $tag, $seq_header, $bamgout_file, $dna_name, $orf_Number];
		}

	}
	$file_id ++;
    }
}
#test
#foreach (sort keys %aa) {
#
#	print $_,"--->\n\t",$aa{$_}->[0],"\n\t",$aa{$_}->[2],"\n";
#}
#foreach (sort keys %gene_list) {
#
#	print $_,"\n";foreach (@{$gene_list{$_}}) {print "\t",$_,"\n"};
#
#}


if (@bamgout_in and defined $errout){ 

	open ERR, ">$errout" or die ("ERROR: Cannot open error file $errout.\n");
	my $nt;
	my $local_predic_aa;
	my $gene_name;
	my $dna_name;
	foreach $gene_name (sort keys %aa) {
	
		$dna_name = $aa{$gene_name}->[5];
		my ($start, $end) = ($aa{$gene_name}->[1] =~ /(\d+)\D+(\d+)/);
		if ($start < $end) {
			$nt = substr($dna{$dna_name}->[1], $start-1, $end-$start+1);
		} else {
			$nt = revcomp (substr ($dna{$dna_name}->[1], $end-1, $start-$end+1));
		}

		$local_predic_aa = translate($nt);
		$local_predic_aa =~ s/[*]$//;

		next if ($local_predic_aa eq $aa{$gene_name}->[0] or $local_predic_aa eq $aa{$gene_name}->[0]."*");
		print ERR "GENE NAME: $gene_name\nSOURCE: ",$aa{$gene_name}->[5],"\nBamg TRANSLATED SEQ: ",$aa{$gene_name}->[0],"\t",$aa{$gene_name}->[1],"\nLOCAL TRANSLATED SEQ: ",$local_predic_aa,"\n";
		print ERR "Nucleotide from the predicted file--> \n$nt\n";
	}

	close ERR;
}

## output

my $date = uc(strftime("%d-%b-%Y", localtime ()));
my $count_dna_id = 0;
# each DNA has a corresponding gbk file except those have no predicted ORFs.



die "ERROR: Genbank files specified unequal to the number of fasta files.\n\tSolve: specify genbankfiles with a same number to fasta files OR leave --out parameter free (the program would name genbank files according to the names of the fasta files).\n" if (@output_files and (scalar @output_files != scalar @fasta_in)) ;
my $count_file_id = 0;

foreach $fasta_file (@fasta_in) {
    my $output;
    if (@output_files and scalar @output_files == scalar @fasta_in){
		$output = $output_files[$count_file_id];
    } else {
	$output = $fasta_file;
	$output =~ s/\.[^.]+$/\.gbk/;
    }
    $count_file_id ++;
    print STDERR "\n****************************************************************\nCurrently processing $fasta_file...Corresponding Genbank file $output\n";
    if (defined $fasta_file_dna{$fasta_file}) {
    	open OUT, "> $output" or die ("ERROR: Cannot open output file $output.\n");
	} else {
#		print STDERR "No ORFs are predicted in fasta file $fasta_file.\n";
		print STDERR "No sequences contained in fasta file $fasta_file.\n";
		next;
	}

    my $total_dna_count=@{$fasta_file_dna{$fasta_file}};
    my $dna_count=0;
    my $start=time();
    my $last=$start-1;
    foreach my $dna_name (@{$fasta_file_dna{$fasta_file}}){
        $dna_count++;

	my @genes = ();
	if (exists $gene_list{$dna_name}) { @genes = @{$gene_list{$dna_name}};} # generating genbank files for non-orf sequences too.
	$count_dna_id ++;

	if ($last!=time()) {
	  $last=time();
	  my $elapsed=(time()-$start);
	  my $percent=$dna_count/$total_dna_count;
	  my $eta=$start+$elapsed/$percent;
	  my $etatext=strftime("%Y-%m-%d %H:%M:%S",localtime($eta));
	  #print STDERR "\r    Currently processing DNA $dna_name... ";
	  print STDERR sprintf("\r    Currently processing DNA %s... %.2f%% done,  ETA: %s    ",$dna_name,100*$percent,$etatext);
	}
	printf OUT "LOCUS       %-18s  %4d bp   DNA           BCT      %-11s\n",$dna_name,  length($dna{$dna_name}->[1]), $date;
	print OUT mywrap("DEFINITION  $definition\n", 12);
	print OUT mywrap("ACCESSION   $accession\n", 12);
	print OUT mywrap("VERSION     $version\n", 12);
	print OUT "DBLINK      $dblink\n";
	print OUT mywrap("KEYWORDS    $keywords\n", 12);
	print OUT "SOURCE";
	print OUT mywrap("      $organism\n", 12);
	print OUT mywrap("  ORGANISM  $organism.\n", 12);
	my $source_organism = join("; ",@{$para_info->{'Taxonomy'}});
	print OUT mywrap("            $source_organism\n", 12);
	my $ref_id = 1;

	foreach my $ref (@{$para_info->{'References'}}){
		print OUT "REFERENCE   $ref_id";
		if (defined $ref->{'Bases'}) {
			print OUT "  (bases ", $ref->{'Bases'}, ")\n";
		} 
		if (defined $ref->{'Authors'}){ #### Need to add more restrictions here for Name processing.
			foreach (@{$ref->{'Authors'}}) {
				my $previous_name = $_;
#				print "$_ in previous format\t";
				if ($_ =~ /^(\w)+,\s*(\w\.\s*)+$/) {				
					$_ =~ s/\s+//;
#					print "$_ after processed\n";
				}elsif ($_ =~ /^(\w+\s+)+\w+$/) { 
					my @name_parts = split (/\s+/, $_);
					$_ = $name_parts[-1].',';
					foreach my $name_part (@name_parts) {
						$_ .= uc substr ($name_part, 0, 1).'.';
					}
					chop $_;
					chop $_;
#					print "$_ after processed.\n";
					print STDERR "\n\t------------------------------------\n\tWARNINGS: Author's name $previous_name (reference $ref_id) is processed with western name conventions ($_). If not, please modify the names directly in short formats.\n";
				}elsif ($_ =~ /^(\w+,)\s*((\w+\s*)+)$/){
					my @name_parts = split (/\s+/, $2);
					$_ = $1;
					foreach my $name_part (@name_parts) {
						$_ .= uc substr ($name_part, 0, 1).'.';
					}
#					print "$_ after processed\n";
				}elsif ($_ =~ /^(\w+,)\s*((\w+\s+)+)((\w\.)+)$/){
					my @name_parts = split (/\s+/, $2);
					$_ = $1;
					foreach my $name_part (@name_parts) {
						$_ .= uc substr ($name_part, 0, 1).'.';
					}
					$_ .= $4;
#					print "$_ after processed\n";
				}else {
					print STDERR "\n\t------------------------------------\n\ERROR: Author's name $_ (reference $ref_id) is unrecognizable. Please modify it in short name format.\n";
					die;
				}
			}
			my $authors = join(", ", @{$ref->{'Authors'}});
			print OUT mywrap("  AUTHORS   $authors\n", 12);
		}
		if (defined $ref->{'Consortium'}){
			print OUT mywrap("  CONSRTM   $ref->{'Consortium'}\n",12);
		}

		if (defined $ref->{'Title'}){
			print OUT mywrap("  TITLE     $ref->{'Title'}\n", 12);
		}

		if (defined $ref->{'Publications'}){
			die ("ERROR: Parameter file reference $ref_id has publication information without journal name.\n") unless (defined $ref->{'Publications'}->{'Journal'});
			die ("ERROR: Parameter file reference $ref_id has publication information without pubmed id.   \n") unless (defined $ref->{'Publications'}->{'PubMed ID'});		
			print OUT mywrap("  JOURNAL   $ref->{'Publications'}->{'Journal'}\n", 12);
			print OUT mywrap("   PUBMED   $ref->{'Publications'}->{'PubMed ID'}\n", 12);
		
		} elsif (defined $ref->{'Direct submission'}) {
			die ("ERROR: Parameter file reference $ref_id has direct submission information without submission date.\n") unless (defined $ref->{'Direct submission'}->{'Date'});
			die ("ERROR: Parameter file reference $ref_id has direct submission information without submission address.\n") unless (defined $ref->{'Direct submission'}->{'Address'});	
			$ref->{'Direct submission'}->{'Date'} = uc ($ref->{'Direct submission'}->{'Date'});
			if ($ref->{'Direct submission'}->{'Date'} =~ /\d{2}-(\w{3})-\d{4}/) {
				die "ERROR: Unrecognizable month in DATE parameter. Month should be from JAN to DEC (not $1).\n" unless ($1 eq 'FEB' || $1 eq 'JAN' || $1 eq 'MAR' || $1 eq 'APR' || $1 eq 'MAY' ||$1 eq 'JUN' ||$1 eq 'JUL' ||$1 eq 'AUG' ||$1 eq 'SEP' ||$1 eq 'ORG' || $1 eq 'NOV' ||$1 eq 'DEC');
			} elsif ($ref->{'Direct submission'}->{'Date'} =~ /(\w+)\s+(\d{1,2}),\s?(\d{4})/i) {
				my $day =($2 =~ /^\d$/)? ('0'.$2):$2;
				my $year = $3;
#				print "day is $day, year is $year\n";
				foreach my $month (@months) {
					if ($1 =~ /$month/i){
						$ref->{'Direct submission'}->{'Date'} = $day.'-'.uc substr ($month, 0, 3).'-'.$year;
						next;
					}
				}
			} else{
				print STDERR "\n\t------------------------------------\n\tWARNINGS: The date format in parameter file is $ref->{'Direct submission'}->{'Date'},\n";
				print STDERR "\tWARNINGS: The format of DATE parameter in parameter file should be transformed into format:  \"DD-MON-YEAR\" OR \"MONTH DD, YYYY\".\n\t------------------------------------\n";
				die "\t------------------------------------\n\tERROR: date format ERROR: Unrecognizable date in parameter file \'$ref->{'Direct submission'}->{'Date'}\',\n";
			}
			print OUT mywrap("  TITLE     Direct Submission\n",12);
			print OUT mywrap("  JOURNAL   Submitted ($ref->{'Direct submission'}->{'Date'}) $ref->{'Direct submission'}->{'Address'}\n", 12);		
		
		}
	    if (defined $ref->{'Remark'}){
			print OUT mywrap("  REMARK    $ref->{'Remark'}\n", 12);
		}


		$ref_id ++;
	}

	if (defined $para_info->{'Comments'}) {

		my $comments = join("\n", @{$para_info->{'Comments'}});
		$comments = "COMMENT     ".$comments."\n";
		print OUT mywrap($comments, 12);
	
	}

#	my @extra_info = grep {!@default_paras} keys %{$para_info};
#	print @extra_info;


	print OUT "FEATURES             Location/Qualifiers\n";
	my $length = length($dna{$dna_name}->[1]);
	print OUT "     source          1..$length\n";
	print OUT mywrap("                     /organism=\"$organism\"\n", 21);
	print OUT mywrap("                     /strain=\"$strain\"\n", 21);
	print OUT mywrap("                     /chromosome=\"$para_info->{'Chromosome'}\"\n", 21) if (defined $para_info->{'Chromosome'});



    if (@genes) {

	foreach my $gene_name (@{$gene_list{$dna_name}}) {
	
        #print Dumper($aa{$gene_name});
		my $locus_tag = sprintf '%s_%08d_%05d',$locus_tag_base, $aa{$gene_name}->[5], $aa{$gene_name}->[6];
        #print $locus_tag."\n";
		print OUT mywrap("     gene           ", 21);

		$aa{$gene_name}->[1] =~ m/(\d+)\.\.(\d+)/;		
		if ($1<$2) {
			print OUT "$1..$2\n";
		} else {
			print OUT "complement($2..$1)\n";
		}
		print OUT mywrap("                     /locus_tag=\"$locus_tag\"\n",21);
		print OUT mywrap("     CDS             ",21);
		if ($1<$2) {
			print OUT "$1..$2\n";
		} else {
			print OUT "complement($2..$1)\n";
		}

		print OUT mywrap("                     /locus_tag=\"$locus_tag\"\n", 21);
		print OUT mywrap("                     /codon_start=$para_info->{'Codon start'}\n", 21) if (defined $para_info->{'Codon start'});
		print OUT mywrap("                     /transl_table=$gcode\n",21);		
		print OUT mywrap("                     /translation=\"$aa{$gene_name}->[0]\"\n",21);
	
	}


    }




	print OUT "ORIGIN\n";
	print OUT dnawrap($dna{$dna_name}->[1]);
	print OUT "//\n";
	
    }
    close OUT;

}
print STDERR "\n";


#Test
#	foreach (sort keys %aa){
#
#		print "$_===>\n","\t",$aa{$_}->[0],"\n\t",,$aa{$_}->[1],"\n\t",,$aa{$_}->[2],"\n\t",,$aa{$_}->[3],"\n\n\n\n";
#
#	}
#
#
#	foreach (sort keys %gene_list) {
#
#		print "DNA_name: $_====>\n";
#		foreach (@{$gene_list{$_}}){
#
#			print "\t $_\n";
#		}
#
#
#	}

exit 0;
sub revcomp {
    print "-----\n";
    print $_[0]."\n";
	my $out = reverse($_[0]);
	$out =~ y/atcgATCG/tagcTAGC/;
    print $out."\n";
	return $out;
}


sub translate {  # Input: dna and codon table, Output: protein
	my ($in) = @_;
	$in = uc($in);
	my $l = length($in);
	my @c;
	while(defined $in and $in ne "") {
		push @c,substr($in,0,3);
		if (length($in)<3) {last;}
		$in=substr($in,3);
	}
	my @aa = grep {defined} @codons{@c};  # Slice!
	return join("",@aa);
}

sub mywrap {
  my ($in,$ks)=@_;
  my $spc="";
  if ($in=~s/^(\s+)//) {$spc=$1;} # Preserve initial leading space
  my $out=wrap($spc," "x$ks,$in);  # Hardcoded 21
  $out=~s/\t/        /g;
  return $out;
}

sub dnawrap {
	my $in=$_[0];
	$in = lc ($in);
	my $out="";
	my $p=1;
	while (defined $in and $in ne "") {
		my $l=substr($in,0,60);
		if (length($in)<60) {$in="";} else {$in=substr($in,60);}
		$l=~s/(\S{10})/$1 /g;
		$l=~s/ $//;
		$out.=sprintf("%9d %s\n",$p,$l);
		$p+=60;  
	}
	return $out;
}


sub get_parameters() {
      my $sample_name = shift;
      
      my %para_info;
      $para_info{'Organism'}=$sample_name;
      $para_info{'Version'}=1.0;
      @{$para_info{'Taxonomy'}}=("Metagenome");
      $para_info{'TaxonId'}='12908';
      $para_info{'Definition'}=$sample_name;
      $para_info{'Chromosome'}='1';
      $para_info{'Accession'}='.';
      $para_info{'Strain'}=1;
      $para_info{'GI number'}=".";
      $para_info{'Database link'}= $sample_name;
      $para_info{'Version'}=".";
      $para_info{'Keywords'}=".";
      $para_info{'Translation table'}=11;
      $para_info{'Codon start'}=1;
      @{$para_info{'References'}} = ( 
                                   {
                                     'Bases' => '1 to XXXXX',
                                     'Authors' => ['YYYYYY, X.'],
                                     'Consortium' => 'XXXXX',
                                     'Title'=>'XXXXX',
                                     'Publications'=> {'Journal'=>'XXXXX', 'PubMed ID'=>'XXXXX'},
                                     'Remark'=>'XXXXX'
                                    }
                                  );

      @{$para_info{'Comments'}} =  (
                                     'PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI  review \n
                                     COMPLETENESS: XXXXX'
                                   ) ;


     return \%para_info;

}
exit 0;



__DATA__
--**************************************************************************
--  This is the NCBI genetic code table
--  Initial base data set from Andrzej Elzanowski while at PIR International
--  Addition of Eubacterial and Alternative Yeast by J.Ostell at NCBI
--  Base 1-3 of each codon have been added as comments to facilitate
--    readability at the suggestion of Peter Rice, EMBL
--  Later additions by Taxonomy Group staff at NCBI
--
--  Version 3.9
--     Code 14 differs from code 9 only by translating UAA to Tyr rather than
--     STOP.  A recent study (Telford et al, 2000) has found no evidence that
--     the codon UAA codes for Tyr in the flatworms, but other opinions exist.
--     There are very few GenBank records that are translated with code 14,
--     but a test translation shows that retranslating these records with code
--     9 can cause premature terminations.  Therefore, GenBank will maintain
--     code 14 until further information becomes available.
--
--  Version 3.8
--     Added GTG start to Echinoderm mitochondrial code, code 9
--
--  Version 3.7
--     Added code 23 Thraustochytrium mitochondrial code
--        formerly OGMP code 93
--        submitted by Gertraude Berger, Ph.D.
--
--  Version 3.6
--     Added code 22 TAG-Leu, TCA-stop
--        found in mitochondrial DNA of Scenedesmus obliquus
--        submitted by Gertraude Berger, Ph.D.
--        Organelle Genome Megasequencing Program, Univ Montreal
--
--  Version 3.5
--     Added code 21, Trematode Mitochondrial
--       (as deduced from: Garey & Wolstenholme,1989; Ohama et al, 1990)
--     Added code 16, Chlorophycean Mitochondrial
--       (TAG can translated to Leucine instaed to STOP in chlorophyceans
--        and fungi)
--
--  Version 3.4
--     Added CTG,TTG as allowed alternate start codons in Standard code.
--        Prats et al. 1989, Hann et al. 1992
--
--  Version 3.3 - 10/13/95
--     Added alternate intiation codon ATC to code 5
--        based on complete mitochondrial genome of honeybee
--        Crozier and Crozier (1993)
--
--  Version 3.2 - 6/24/95
--  Code       Comments
--   10        Alternative Ciliate Macronuclear renamed to Euplotid Macro...
--   15        Bleharisma Macro.. code added
--    5        Invertebrate Mito.. GTG allowed as alternate initiator
--   11        Eubacterial renamed to Bacterial as most alternate starts
--               have been found in Achea
--
--
--  Version 3.1 - 1995
--  Updated as per Andrzej Elzanowski at NCBI
--     Complete documentation in NCBI toolkit documentation
--  Note: 2 genetic codes have been deleted
--
--   Old id   Use id     - Notes
--
--   id 7      id 4      - Kinetoplast code now merged in code id 4
--   id 8      id 1      - all plant chloroplast differences due to RNA edit
--
--*************************************************************************

Genetic-code-table ::= {
 {
  name "Standard" ,
  name "SGC0" ,
  id 1 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M---------------M---------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Vertebrate Mitochondrial" ,
  name "SGC1" ,
  id 2 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
  sncbieaa "--------------------------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Yeast Mitochondrial" ,
  name "SGC2" ,
  id 3 ,
  ncbieaa  "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "----------------------------------MM----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
    name "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate
 Mitochondrial; Mycoplasma; Spiroplasma" ,
  name "SGC3" ,
  id 4 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--MM---------------M------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Invertebrate Mitochondrial" ,
  name "SGC4" ,
  id 5 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "---M----------------------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear" ,
  name "SGC5" ,
  id 6 ,
  ncbieaa  "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Echinoderm Mitochondrial; Flatworm Mitochondrial" ,
  name "SGC8" ,
  id 9 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Euplotid Nuclear" ,
  name "SGC9" ,
  id 10 ,
  ncbieaa  "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Bacterial, Archaeal and Plant Plastid" ,
  id 11 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "---M---------------M------------MMMM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Alternative Yeast Nuclear" ,
  id 12 ,
  ncbieaa  "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-------------------M---------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Ascidian Mitochondrial" ,
  id 13 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
  sncbieaa "---M------------------------------MM---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 },
 {
  name "Alternative Flatworm Mitochondrial" ,
  id 14 ,
  ncbieaa  "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Blepharisma Macronuclear" ,
  id 15 ,
  ncbieaa  "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Chlorophycean Mitochondrial" ,
  id 16 ,
  ncbieaa  "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Trematode Mitochondrial" ,
  id 21 ,
  ncbieaa  "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Scenedesmus obliquus Mitochondrial" ,
  id 22 ,
  ncbieaa  "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "-----------------------------------M----------------------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 } ,
 {
  name "Thraustochytrium Mitochondrial" ,
  id 23 ,
  ncbieaa  "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  sncbieaa "--------------------------------M--M---------------M------------"
  -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 }
}
