#! /usr/bin/perl
use strict;
use warnings;
my $file=shift;
my %seq;
#define search/align kmer size
my $kmer=shift;
my %asm_kmer;
open(IN,$file)||die $!;
my $id=0;
while(<IN>){
    chomp;
    if($_=~/^>/){
	$id++;
    }else{
	$seq{$id}=$_;
	$seq{$_}=$id;
	seqkmerindex($_,$kmer);
    }
}

#
my %merged_path;
my %merged_pathA;
my $mid=0;


my %edge_col;

foreach my $kmerseq (keys %asm_kmer){
    my @askm=@{$asm_kmer{$kmerseq}};
	my $e_id=-1;
	my $alen=scalar(@askm);
	###Create list for edges###
	for(my $x=0;$x<$alen-1;$x++){
		for(my $y=$x+1;$y<$alen;$y++){
			$edge_col{$askm[$x]}{$askm[$y]}=$kmerseq;
			$edge_col{$askm[$y]}{$askm[$x]}=$kmerseq;
		}
	}
    foreach my $ele (@askm){
		if($ele=~/'$/){
			my $nele=substr($ele,0,-1);
			if(exists $merged_path{$nele}){
				$e_id=$merged_path{$nele};
			}
		}else{
			if(exists $merged_path{$ele}){
				$e_id=$merged_path{$ele};
			}
		}
	}
	if($e_id == -1){
		 foreach my $ele (@askm){
			 $merged_path{$ele}=$mid;
		 }
	}else{
		foreach my $ele (@askm){
			 $merged_path{$ele}=$e_id;
		}
	}
	$mid++;
}

foreach my $seqID (keys %merged_path){
	my $smid=$merged_path{$seqID};
	if(exists $merged_pathA{$smid}){
		my @sar=@{$merged_pathA{$smid}};
		push @sar,$seqID;
		$merged_pathA{$smid}=\@sar;
	}else{
		my @sar;
		push @sar,$seqID;
		$merged_pathA{$smid}=\@sar;

	}
}
my %asm_group;
foreach my $s (keys %merged_pathA){
	my @ws=@{$merged_pathA{$s}};
	my @nrws;
	my %rec;
	foreach my $ews (@ws){
		next if(exists $rec{$ews});
		push @nrws,$ews;
		$rec{$ews}=1;
	}
	my $ostr=join(",",@nrws);
	my $wslen=scalar(@nrws);
	print "$s\t$wslen\t$ostr\n";
	$asm_group{$s}=\@nrws;
}

###Building path and assembly ###
my $outputfile=shift;
open(OUT,">$outputfile")||die $!;
my $outputseq=shift;
open(OUT,">$outputfile")||die $!;
open(OSEQ,">$outputseq")||die $!;
foreach my $gid (keys %asm_group){
	my @ws=@{$asm_group{$gid}};
	my %seqidcol;
	foreach my $seqid (@ws){
		$seqidcol{$seqid}=1;
	}
	print  "\n\n\nPATHID:$gid\n";
	print OUT "\nPATHID:$gid\nALN\n";
	my $ostr=join(",",@ws);
	print "$ostr\n";
	my @s=path_search(%seqidcol);
	my %asml=%{$s[0]};
	my %outputasml=%{$s[1]};
	my $conseq="";
	my $prev_len=0;
	foreach my $seqid (sort {$outputasml{$a} <=> $outputasml{$b}} keys %outputasml){
		my $oid=sprintf("%04d",$seqid);
		print OUT "$oid\t$asml{$seqid}\n";
		if(length($conseq)<1 ){
			next if($asml{$seqid}=~/^\./);
			$conseq=$asml{$seqid};
		}else{
			my $exlen=$outputasml{$seqid}-length($conseq);
			next if($exlen ==0);
			my $rightexseq=substr($asml{$seqid},-$exlen);
			next if($rightexseq=~/^\./);
			$conseq=$conseq.$rightexseq;
		}
	}
	my $ogid=sprintf("%03d",$gid);
	print OUT "C$ogid\t$conseq\n";
	my $conlen=length($conseq);
	print OSEQ ">C$ogid\|$conlen\n$conseq\n";
}
close(OUT);
sub path_search{
	(my %seqidcol)=@_;
	my @hash_keys = keys %seqidcol;
	#extension on both sides;
	my %path;
	my $random_seqID = $hash_keys[rand @hash_keys];
	my $lsc=scalar(@hash_keys);
	print "START:$random_seqID\tNumber of seq:$lsc\n";
	my $seedseq = $seqidcol{$random_seqID};
	my %color;
	my $step=0;
	$color{$random_seqID}=$step;
	#$path{0}=$random_seqID;
	my $tosearch=1;
	$step++;
	my %linectrl;
	my %lastcoor;
	my $start_x=0;
	my %asml;
	my %oxswitch;
	my %output_asml;
	$asml{$random_seqID}=$seq{$random_seqID};
	while($tosearch==1){
		my @dotlist;
		my $number_newdot=0;
		foreach my $dot (keys %color){
			if($color{$dot} == $step-1){
				push @dotlist,$dot;
			}
		}
			
		foreach my $dot (@dotlist){
			foreach my $neighbor (keys %{$edge_col{$dot}}){
				next if(exists $color{$neighbor});
				next if(not exists $seqidcol{$neighbor});
				my $kmerseq=$edge_col{$dot}{$neighbor};
				my @dotpos=seqkmerQuery($seq{$dot},$kmerseq);
				my @neighborpos=seqkmerQuery($seq{$neighbor},$kmerseq);
				push @dotpos,@neighborpos;
				$path{$dot}{$neighbor}=\@dotpos;
				$path{$neighbor}{$dot}=\@dotpos;
				$color{$neighbor}=$step;
				$number_newdot++;
				my $t1seq=$seq{$dot};
				my $d1seq=$seq{$neighbor};
				my $t1start=$dotpos[0];
				my $t1ori=$dotpos[1];
				my $d1start=$neighborpos[0];
				my $d1ori=$neighborpos[1];
				if($step>1){
					$t1seq=$asml{$dot};
					my $npad=length($t1seq)-length($seq{$dot});
					my $qseq=substr($t1seq,$npad,length($seq{$dot}));
					if(exists $oxswitch{$t1seq}){
						$t1seq="."x$npad.revcom($t1seq);
					}
					@dotpos=seqkmerQuery($t1seq,$kmerseq);
					$t1start=$dotpos[0];
					$t1ori=$dotpos[1];
					$start_x=0;
				}
				print "#$dot\t$neighbor\t$t1start\t$t1ori\t$d1start\t$d1ori\t$start_x\t$kmerseq\n";
				#print "#$t1start\t$d1start\n";
				if($t1ori+$d1ori == 0){
					$d1seq=revcom($d1seq);
					$d1start=length($d1seq)-($d1start+$kmer);
					$oxswitch{$neighbor}=1;
					
				}
				if($t1start>$d1start){
					#When Reftop on the left,bottom query on the right
					#padding to the bottom query
					## AAAATT
					## ..AATT
					my $str=('.')x($t1start-$d1start);
					my $Ldlpadding=('.')x($start_x+$t1start-$d1start);
					$asml{$neighbor}=$Ldlpadding.$d1seq;
					$lastcoor{$neighbor}=$start_x+$t1start-$d1start;
				}else{
					#When Reftop on the left,bottom query on the right
					#padding to the top seq
					## .....AATT
					## AAAAAAATT
					my $str=('.')x($d1start-$t1start);
					if($d1start-$t1start>$start_x){
					## .....AATT
					## AAAAAAATT
					## .......AATT
					## ..AAAAAAATT
					## CFAAAAAAATT
						my $addx=$d1start-$t1start-$start_x;
						my $i=0;
						foreach my $smid (keys %asml){
							my $astr=('.')x$addx;
							my $smstr=$astr.$asml{$smid};
							$asml{$smid}=$smstr;
						}
						$start_x=$d1start-$t1start;
						$asml{$neighbor}=$d1seq;
						$lastcoor{$neighbor}=0;
					}else{
					## .....AATT
					## AAAAAAATT
					## ...AAAATT
						my $Ldlpadding=('.')x($start_x-$d1start+$t1start);
						$asml{$neighbor}=$Ldlpadding.$d1seq;
						$lastcoor{$neighbor}=$start_x-$d1start+$t1start
					}

				}
				#ALIGN the strings#
			}
		}
		if($number_newdot == 0){
			$tosearch=0;
			print "STEP:$step\n";
			my $asm_seq="";
			my $linec=0;
			foreach my $colk (sort {$color{$a} <=> $color{$b}} keys %color){
				my $ss=sprintf("%04d",$colk);
				print "$ss\t$color{$colk}\t$seq{$colk}\t";
				
				if(exists $linectrl{$colk}){
					print "Unexpected reuse of seq$colk,exit\n";
					exit;
				}
				print "$asml{$colk}\n";
				$output_asml{$colk}=length($asml{$colk});
				$linectrl{$colk}=$linec;
				$linec++;
				
			}
			foreach my $colk (sort {$color{$a} <=> $color{$b}} keys %color){
				print "$colk->";
				foreach my $bb (keys %{$path{$colk}}){
					next if($color{$colk}>$color{$bb});
						print "$bb,";
				}
				print "\n";
			}
			
		}
		$step++;
	}
	return \%asml,\%output_asml;
}

sub backtrackseq{
	
}
sub seqkmerQuery{
    (my $tseq, my $qkseq)=@_;
	my $k=length($qkseq);
	my $qstart=-10;
	my $ori=1;
    for(my $i=0;$i+$k<=length($tseq);$i++){
		my $kseq=substr($tseq,$i,$k);
		if($kseq eq $qkseq){
			$qstart=$i;
			last;
		}else{
			my $rcseq=revcom($kseq);
			if($rcseq eq $qkseq){
				$qstart=$i;
				$ori=-1;
				last;
			}
		}
	}
	if($qstart==-10){
		print "Couldn't find the kmer for the pair!\n";
		print "$tseq\t$qkseq\n";
		exit;
	}
	return $qstart,$ori;
}

sub seqkmerindex{
    (my $tseq, my $k)=@_;
    for(my $i=0;$i+$k<=length($tseq);$i++){
		my $kseq=substr($tseq,$i,$k);
		my $seqid=$seq{$tseq};
		if(exists $asm_kmer{$kseq}){
			my @kmerToSeq=@{$asm_kmer{$kseq}};
			push @kmerToSeq,$seqid;
			$asm_kmer{$kseq}=\@kmerToSeq;
		}else{
			my $rcseq=revcom($kseq);
			if(exists $asm_kmer{$rcseq}){
				my @kmerToSeq=@{$asm_kmer{$rcseq}};
				push @kmerToSeq,$seqid;
				$asm_kmer{$rcseq}=\@kmerToSeq;
			}else{
				my @kmerToSeq;
				push @kmerToSeq,$seqid;
				$asm_kmer{$kseq}=\@kmerToSeq;
			}
		}
    }
}



sub revcom{
    (my $seq)=@_;
    $seq=~tr/atcgATCG/tagcTAGC/;
    my $rcseq=reverse $seq;
    #print $rcseq,"\n";
    return $rcseq;
}
