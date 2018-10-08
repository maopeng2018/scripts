open(geneList,"genes_to_check");
open(allGene,"combined_expression_orthologous_geneName_all.txt");

open(result,">sequence_for_check_exact_CAZY_substrate.fa");
open(result2,">gene_names_for_check_exact_CAZY_substrate.fa");


my %ortholog;

while(<geneList>){
chomp;
my @line=split(/\|\|/,$_);
my $gene=$line[1];
$ortholog{$gene}="";
# print "$gene\n";
}


my %genes;
while(<allGene>){
my $record=$_;
chomp;
my @line=split(/\t/,$_);
if(exists $ortholog{$line[0]}){
## print result "$record";
my @tmp=split(/\t|\s/,$record);
my $ortho=shift(@tmp);
pop(@tmp);

foreach(@tmp){
     if($_=~/^NA$/){}else{
  my $gene="$_\@$ortho";
  $genes{$_}=$gene;
 }
    }
  }
}



foreach(keys %genes){
 print result2 "$genes{$_}\n";
}



$dir=".";
opendir(DIR,   $dir."\/")   ||   die   "Can't   open   directory   $dir";   
   my @files   =   readdir(DIR); 
 
my @file2;
foreach(@files){
if($_=~/.fasta$/){
push(@file2,$_);
print "$_\n";
}
}



foreach(@file2){
my $file=$_;
open(seq,$file);
my $label=0;

while(<seq>){
   my $record=$_;
   if($_=~/^>/){
   my @tmp=split(/\|/,$_);
   $tmp[0]=~s/>//;
   $id="$tmp[0]\|$tmp[1]\|$tmp[2]\|";
   ##print "$id\n";
   if(exists $genes{$id}){
   $record=">$genes{$id}\n";
   $label=1;
   }else{$label=0}
  }

if($label==1){
  print result "$record";
 }
 
}
close(seq);
}






