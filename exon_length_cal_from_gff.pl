open(file,"Pycco1_GeneCatalog_genes_20140114.gff");
open(result,">Pyccol_exon_length.txt");
open(result2,">tmp_Pyccol_exon_length.txt");

print result "sequencing\tproteinId\ttranscriptId\texon_length\n";

my %hash_tran;
my %hash_pro;
my %hash;

while(<file>){
chomp;
my @line=split(/\t/,$_);

if($_=~/name\s"/){
my $name=(split(/"/,$line[8]))[1];
## print "$name\n";

if($_=~/transcriptId/){
my $id=(split(/transcriptId/,$line[8]))[1];
$id=~s/\s|;//g;
$hash_tran{$name}=$id;

my $len=$line[4]-$line[3]+1;
if(exists $hash{$name}){$hash{$name}=$hash{$name}+$len;
  print result2 "$id\t$line[3]\t$line[4]\t$len\n";
}else{
    $hash{$name}=$len;	
  print result2 "$id\t$line[3]\t$line[4]\t$len\n";
  }
}

if($_=~/proteinId/){
my $id=(split(/proteinId|exonNumber/,$line[8]))[1];
$id=~s/\s|;//g;
$hash_pro{$name}=$id;
}

}


}




foreach(keys %hash){
my $name=$_;
print result "$name\t$hash_pro{$name}\t$hash_tran{$name}\t$hash{$name}\n";

}

