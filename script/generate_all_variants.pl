#!/user/bin/perl

use strict;
use warnings;

our %rate=();
our %gene_rate=();
our $cutoff=0.025;
our @nucleotides=("A","T","C","G");

our %trichanges=();
our %amino_tables=();
my $path="src/";
my $file_rate="$path/3mer_table.txt";
## read amino acid change table
my $file_aminoacid="$path/amino_acid_codon.txt";


my %complement=(
"A" =>"T",
"T" =>"A",
"C" => "G",
"G" => "C",
"a" =>"t",
"t"=>"a",
"c"=>"g",
"g"=>"c"	
);

sub main {
		
		my $gene="";
		my $range="";
		my $ftranscript=$_[0];
		my $fout=$_[1];
		print "$ftranscript\t$fout\n";

		my %gene_rate=();
	   open my $OUT, ">$fout";
		print $OUT "#CHROM\tPOS\tREF\tALT\tCODON\tCODONChange\tAA\tAAChange\tTranscriptID\tExon\tExonicFunc\tMutationRate\tTranscriptID\tStrand\n";
	  my $cmd="cat";
	  if($ftranscript =~/.gz$/){$cmd="zcat";}
	       	open IN, "$cmd $ftranscript|" or die $ftranscript." cannot find!\n";
	   my @lastgenes=();
		my @lastpos=();
		my $line=<IN>;
	    my $strseq="";
		my @exon_starts=();
        my @exon_ends=();
		my $strand="+";
		do{
			if(!$line){last;}
			$line=~ s/range=chr//g;
			$line=~ s/strand=//g;
			my @sets=split(/\s+/,$line);
			if(@sets<6){die "uncomplete info in $line ";}
			my @pos_info=split(/:|-/,$sets[1]); ## $pos_info[1] start, $pos_info[2] end
			my @gene_info=split(/_|\./,$sets[0]);
			if(@gene_info==5){$gene_info[3]="";}
			my $chr=$pos_info[0]; ##chr
			## out put #chr pos ref alt  codon alt_codon mutation rate
			$strand=$sets[4];   	## depend on - or +
			## read next line until start with :\>:
			$line=<IN>;
			my $exonseq="";
			while($line&& !($line =~ /^>/)){
				chomp($line);
				$exonseq=$exonseq."".$line;
 				#print length($line);
				$line=<IN>;
			}
			## exonseq
			my @seqs=split("",$exonseq);
			if($strand eq "-"){
			## from end to start and reverse

			for(my $i=1;$i<@seqs-1;$i++){
			    my $j=$pos_info[2]-$i;
			    foreach my $nc(@nucleotides){
			    	if($nc eq uc($complement{$seqs[$i]})){next;}
			    	my $codon=$complement{$seqs[$i+1]}.$complement{$seqs[$i]}.$complement{$seqs[$i-1]};
			    	my $alt_codon=$complement{$seqs[$i+1]}.$nc.$complement{$seqs[$i-1]};
			    	my $rate=0;
			    	if(exists($trichanges{uc($codon).">".uc($alt_codon)})){$rate=$trichanges{uc($codon).">".uc($alt_codon)};}
					print $OUT $pos_info[0]."\t".$j."\t".uc($complement{$seqs[$i]})."\t".$nc."\t".$codon."\t".$alt_codon."\t".$rate."\t".$gene_info[2]."\t".$strand."\n";
				}
			}
			}else{
			## from start to end and reverse
			
			for(my $i=1;$i<@seqs-1;$i++){
			    my $j=$pos_info[1]+$i;
			    foreach my $nc(@nucleotides){
			    	if($nc eq uc($seqs[$i])|| $nc eq $seqs[$i]){next;}
			    	my $codon=$seqs[$i-1].$seqs[$i].$seqs[$i+1];
			    	my $alt_codon=$seqs[$i-1].$nc.$seqs[$i+1];
			    	my $rate=0;
			    	if(exists($trichanges{uc($codon).">".uc($alt_codon)})){$rate=$trichanges{uc($codon).">".uc($alt_codon)};}
					print $OUT $pos_info[0]."\t".$j."\t".uc($seqs[$i])."\t".$nc."\t".$codon."\t".$alt_codon."\t".$rate."\t".$gene_info[2]."\t".$strand."\n";
				}
			}
			
			}
		} while($line);
		#		my $rs=callinfo($strseq,$lastgenes[2],\@lastpos, \@exon_starts,\@exon_ends,$strand,$OUT);
		close IN;
	   close $OUT;
		
	}
	
	
	
	sub load_info {
				
		### transcript 
		
		## load amino acide change
		if (! -e $file_rate) {print "$file_rate cannot find\t";exit;}	
		open IN, $file_rate or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]/)){next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<3){die "incomplete information in $line ";}
			$trichanges{$sets[0].">".$sets[1]}=$sets[2];
		}
		close IN;
		
		open IN, $file_aminoacid or die $file_rate." cannot find!\n";
		while(my $line=<IN>){
			if(!($line =~ /^[ACGT]{3}/)) {next;}
			chomp($line);
			my @sets=split(/\s+/,$line);
			if(@sets<4){die "uncomplete line : $line "}
			$amino_tables{$sets[0]}=$sets[3];
		}
		close IN;
	}	
	

	
	

print "Input: perl generate_all_variants.pl  transcript outputFile\n";
#	print "Input perl Mutation_rate_20180814.pl transcript_file outname f\n";
	#print "#filter_with_mappability: 0 or 1, 0 means no mappablity filter, otherwise filter the position with mappablity\n";
my $ftranscript=$ARGV[0]; #"$path/Resources/Transcript/GencodeV19_splice.gz";MY $QMAP=1;
my $fout=$ARGV[1];	
my $QMAP=0;
load_info();

main($ftranscript,$fout);
print "output to $fout\n"

