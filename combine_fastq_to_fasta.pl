#####Created by Mao Peng   12-01-2016, if any question can contact m.peng@cbs.knaw.nl ######
#!/usr/bin/env perl

#!/usr/bin/perl
#!/usr/bin/perl


			  
my $dir=".";
opendir(DIR,   $dir."\\")   ||   die   "Can't   open   directory   $dir";   
   my @files   =   readdir(DIR); 
 
 
 
my @file2;
foreach(@files){
if($_=~/\.fastq$/){
print "$_\n";
push(@file2,$_);
}
}




foreach(@file2){
print "Process $_\n";
my $choose=$_;


# Open input and output files.
open( $in, "<",  $choose)  or die "Can't open $inf: $!";
open( $out, ">>",  "all.fasta")  or die "Can't open $outf.fa: $!";

while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  chomp($temp[2] = <$in>);	# Third line is an id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

  # Substring to inset value.
  $temp[1] = substr($temp[1], $inset);

  # Print to fasta file.
  print $out ">$temp[0]\n";
  print $out "$temp[1]\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# EOF.








}



