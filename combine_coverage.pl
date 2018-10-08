#####Created by Mao Peng   12-01-2016, if any question can contact m.peng@cbs.knaw.nl ######

$dir=".";
opendir(DIR,   $dir."\\")   ||   die   "Can't   open   directory   $dir";   
   my @files   =   readdir(DIR); 
 
my @file2;
foreach(@files){
if($_=~/\.coverage.xls$/){
push(@file2,$_);
}
}

open(result,">all_genes_coverage.txt");


my %hash;
my @sample;

foreach(@file2){
print "Process $_\n";
$choose=$_;
my @tmp=split(/\.coverage\.xls/i,$choose);
my $sample=$tmp[0];
push(@sample,$sample);
open(dat,"$choose");
my $line=0;
while(<dat>){
$line++;
if($line>1){
chomp;
my @line=split(/\s+/,$_);
my $sampleId="$sample\_$line[0]";
$hash{$sampleId}=$line[2];
 }
}

close(dat);
}


open(ids,"all.gene.FPKM.xls");
my @ids;
while(<ids>){
my @line=split(/\s+/,$_);
push(@ids,$line[0]);
}
shift (@ids); ### the first row is not a gene ID;

print result "gene_id";
foreach(@sample){
print result "\t$_";
}
print result "\n";

foreach(@ids){
print result "$_";
my $id=$_;
foreach(@sample){
  my $sampleValue=0;
  my $sampleId="$_\_$id";
  if(exists $hash{$sampleId}){$sampleValue=$hash{$sampleId}};
  print result "\t$sampleValue";
 }
print result "\n";
}

close(result);



