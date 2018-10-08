open(tcdb,"N_crassa_20130412.aa.fasta");
open(results,">protein_length_N_crassa.txt");



my $tag=0;
my $lengths=0;
my $id;
while(<tcdb>){
chomp;
 if($_=~/^>/){
   ## my @tmp=split(/\|/,$_);
   if($tag>0){
   print results "$id\t$lengths\n"; 
   }
	  $id=(split(/>/,$_))[1];
	  ## print "$id\n";
      $tag=0;
      $lengths=0	  
  }else{
      $tag++;
	  $_=~s/\*//g;
	  $lengths=$lengths+length($_);
  }  

    if(eof){
	   print results "$id\t$lengths\n";
    }
  
}

