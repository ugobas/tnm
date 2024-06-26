#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;

# Run: ./script_run_tnm <List_file> <Input_TNM_file>

# path to tnm, modify as needed:
my $tnm="/data/ortizg/software/Linux/ubastolla//tnm"; 
# path to pdb files in CBM cluster:
my $pdbdir="/ngs/databases/pdb/";

my $cluster=0; # run in cluster or in your own machine?

if(scalar(@ARGV)<2){
    print "ERROR, Input list and configuation file must be specified\n";
    print "USAGE: ",$0," <List of PDB codes> <Configuration file>\n";
    die;
}

# List of proteins, modify as needed:
#$List="antibody_dimers.csv";
my $List=$ARGV[0];

# Configuration file:
#$Input="Input_TNM.in";
my $Input=$ARGV[1];

# Working directory, modify as needed:
# $dir="/data/ubastolla/TNM/BINDING";
my $dir=`pwd`; chomp($dir);

open(my $fh, '<:encoding(UTF-8)', $List)
or die "Could not open file '$List' $!";

my $n=0;
while (my $row = <$fh>) {
    # Reading PDB codes and chain identifiers in List
    chomp $row;
   if(substr($row,0,1) eq "#"){next;}
    chomp($row);
    my @res = split(/\s+/, $row);
    my $P1=$res[0];
    # chain codes, modify as needed:
    my $C1=$res[1]; my $P2=$res[2]; my $C2=$res[3]; 
    $n++;

    my $name=$P1;
    if($C1){$name=sprintf("%s%s",$name,$C1);}
    if($P2){$name=sprintf("%s_%s",$name,$P2);}
    if($C2){$name=sprintf("%s%s", $name,$C2);}
    print "Analyzed protein: ",$n," ",$name, "\n";

    my $tmp=sprintf("tmp_%s.in", $name);
    open(my $fo, '>', $tmp);

    print $fo "PDB1= ",$pdbdir,$P1,".pdb\n";
    if($C1){print $fo "CH1= ",$C1,"\n";}
    if($P2){print $fo "PDB2= ",$pdbdir,$P2,".pdb\n";}
    if($C2){print $fo "CH2= ",$C2,"\n";}

    close $fo;

    `cat $Input >> $tmp`;

    # Running the TNM program directly
    # $command=sprintf("%s %s\n", $tnm, $tmp);
    # `$command`; 

    # Sending the program to the cluster
    my $script=sprintf("run_tnm_%s", $name);
    my $dirtmp;
    if($cluster){
	$dirtmp=sprintf("/junk/ubastolla/%s", $name);
	`echo "mkdir " $dirtmp > $script`;
	`echo "mv " $tmp $dirtmp >> $script`;
	`echo "cp Mutation_para.in " $dirtmp >> $script`;
	`echo "cd " $dirtmp >> $script`;
    }
    `echo $tnm $tmp ">>" $name.log >> $script`;
    # Comment the output files that you want to maintain
    `echo "rm -f $P1*MSF*dat" >> $script`;
    `echo "rm -f $P1*Modes*dat" >> $script`;
    `echo "rm -f $P1*Force*dat" >> $script`;
    `echo "rm -f $P1*Bfact*dat" >> $script`;
    `echo "rm -f $P1*coupling.dat" >> $script`;
    `echo "rm -f $P1*sites.in" >> $script`;
    `echo "rm -f $P1*names.dat" >> $script`;
    `echo "rm -f $P1*binding_residues.dat" >> $script`;
    `echo "rm -f $P1*binding_sites.dat" >> $script`;
    `echo "rm -f $P1*pairs_sites.dat" >> $script`;
    `echo "rm -f $name*Confchange.pdb" >> $script`;
    `echo "rm -f $name*regression.dat" >> $script`;
    `echo "rm -f tmp_$name*" >> $script`;
    if($cluster){
	`echo "mv *" $dir >> $script`; # copy everything to main directory
	`echo "cd " $dir >> $script`;
	`echo "rm -rf " $dirtmp >> $script`;
    }
    `echo "rm -f run_tnm_$name*" >> $script`;
    `chmod u+x $script`;
    if($cluster){
	`qsubmit.pl -s "$script" -q x86_64 -n 1`;
    }else{
	`$script`;
    }
}

