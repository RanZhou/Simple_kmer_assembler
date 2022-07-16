#! /usr/bin/perl
use strict;
use warnings;
my $shortseq_file=shift;
my $kmer=shift;
my $sq_count=1;
my $round_count=1;
my $folder_name="round".$round_count;
mkdir($folder_name);
while($sq_count==1){
	if($round_count==1){
		`perl kmer_asm2.pl $shortseq_file $kmer $folder_name/$round_count.aln $folder_name/$round_count.asm.fa > $folder_name/$round_count.debug.log`;
	}else{
		my $last_rd=$round_count-1;
		my $last_fdn="round".$last_rd;
		my $current_fdn="round".$round_count;
		mkdir($current_fdn);
		`perl kmer_asm2.pl $last_fdn/$last_rd.asm.fa $kmer $current_fdn/$round_count.aln $current_fdn/$round_count.asm.fa > $current_fdn/$round_count.debug.log`;
		my $f1_count=count_seq("$last_fdn/$last_rd.asm.fa");
		my $f2_count=count_seq("$current_fdn/$round_count.asm.fa");
		if($f2_count==$f1_count){
			$sq_count=0;
			last;
		}
	}
	$round_count++;
}

sub count_seq{
	(my $infile)=@_;
	open(O1,$infile)||die $!;
		my $seq1_count=0;
		while(<O1>){
			chomp;
			if($_=~/^>/){
				$seq1_count++;
			}
		}
		close(O1);
	return $seq1_count;
}