#!/usr/bin/env perl

use fgweb::general::general ;
use fgweb::general::numerics ;
use fgweb::general::cmd_input ;

# usage: perl matrix2cybert.pl -i=input.txt


my %ini ;
$ini{'-i'}{desc} = 'gene expression matrix' ;
#$ini{'-i'}{val}  = 'Normalized_expression_RMA_alias_function.txt' ;
$ini{'-i'}{checkisexisting} = 1 ;
$ini{'-i'}{checkisstring} = 1 ;
#$ini{'-i'}{mandatory} = 1 ;

$ini{'-cutoff'}{desc} = 'fdr cutoff for the matrix of number of targets for each comparison of classes.' ;
$ini{'-cutoff'}{val}  = '0.05' ;
$ini{'-cutoff'}{checkisstring} = 1 ;
#$ini{'-cutoff'}{mandatory} = 1 ;

$ini{'-fdrmatrix'}{desc} = 'gene expression matrix' ;
$ini{'-fdrmatrix'}{val}  = 'fdr_matrix.txt' ;
$ini{'-fdrmatrix'}{checkisstring} = 1 ;
#$ini{'-fdrmatrix'}{mandatory} = 1 ;

my %groups ;
my %data ;
my @files ;
my $R_cybert = 'hdarray.R' ;

my $cutoff_fdr = 0.05 ;

# functions
# =====================================

sub parseparams {       # shows the usage information if applicable
    my $title = "\nOriginally made by Sacha van Hijum [sacha.vanhijum\@nizo.nl], revised by Miaomiao Zhou [m.zhou\@cbs.knaw.nl]\n\n" ;
    
    $title   .= "Performs Cyber-T analysis on all combinations of replicated samples\n" ;
    cmd_input::check_cmd_input(\%ini, \@ARGV,$title) ;
} # parseparams
        

sub readmatrix { # reads matrix of gene expressions, col1=gene name, row1=experiment name
    my $line ;
    my @splat ;
    my $col ;
    my $exp ;

    open(FILE,$ini{'-i'}{val}) or general::error('could not load gene expression file '.$ini{'-i'}{val}) ;
    while(<FILE>) {
	$line = $_ ;
	chomp $line ;
	@splat = split(/\t/,$line) ;
	if ($. == 1) {	# parse header
	    for ($col = 1 ; $col < scalar @splat ; $col++) {
		push(@{ $groups{$splat[$col]} }, $col) ;
	    }
	} else {	# parse data
	    @{ $data{$splat[0]} } = @splat ;
	}
    }
    close(FILE) ;

#foreach $exp (sort keys %groups) {
 #   print $exp.'=' ;
#    foreach (@{$groups{$exp}}) {
#	print $_.',' ;
#    }
#    print "\n" ;
#}    

} # readmatrix

sub createexpfiles {
    my $exp1 ;
    my $exp2 ;
    my $col ;
    my $id ;
    my $fname ;

    my %done ;

    my $outdata = '' ;
    my $data_dir = 'data/' ;

    my $rep1 ;
    my $rep2 ;
    
    my $fdr_targets ;
    my %fdr_matrix ;
    
    general::makedir('data') ;
    foreach $exp1 (sort keys %groups) {
	foreach $exp2 (sort keys %groups) {
	    if ($exp1 ne $exp2) {
		$fdr_matrix{$exp1}{$exp2} = 0 ;
		$fdr_matrix{$exp2}{$exp1} = 0 ;
	    } else {
		$fdr_matrix{$exp1}{$exp2} = 'NA' ;
	    }
	    if ( $exp1 ne $exp2 and !defined($done{$exp1.$exp2}) and !defined($done{$exp2.$exp1}) ) {
		$fname = $exp2.'_over_'.$exp1.'.txt' ;
		$outdata = 'id_unique' ;
		
    		$rep1 = 0 ;
		$rep2 = 0 ;
		foreach $col ( @{ $groups{$exp1} }) {
		    $rep1++ ;
		    $outdata .= "\t".$exp1 ;
		}
		foreach $col ( @{ $groups{$exp2} }) {
		    $rep2++ ;
		    $outdata .= "\t".$exp2 ;
		}
		$outdata .= "\n" ;
		print 'processing: '.$exp1.' vs '.$exp2.' -> '.$minrep.' replicates ' ;
                print "\n";
		foreach $id (sort keys %data) {
		    $outdata .= $id ;
		    foreach $col ( @{ $groups{$exp1} }) {
			if ($col >= 1) {
			    $outdata .= "\t".@{ $data{$id} }[$col] ;
			}
		    }
		    foreach $col ( @{ $groups{$exp2} }) {
			if ($col >= 1) {
			    $outdata .= "\t".@{ $data{$id} }[$col] ;
			}
		    }
		    $outdata .= "\n" ;
		}
		general::write2file($data_dir.$fname,$outdata) ;
		
		
		
		
		
	    }
	}
    }

    
} # createexpfiles



sub formatratio {
    my $val = shift() ;

    $val = (- 1 / $val ) if ($val < 0) ;

    return $val ;
} # formatratio



&parseparams() ;
&readmatrix() ;
&createexpfiles() ;
