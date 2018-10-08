open(aa,"sequence_for_tcoffee.fasta");


open(result,">sequence_for_tcoffee2.fasta");

my $id;
my $seq;

while(<aa>){
chomp;
 if($_=~/>/ or eof){
	 my @ids=split(/\s+/,$_);
	 print result "$id\n$seq\n";
     $id=$ids[0];
     $seq="";     
	}else{
	$seq="$seq$_";
	}
}