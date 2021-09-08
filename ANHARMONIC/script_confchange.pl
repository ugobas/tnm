#!/usr/bin/env perl 
use strict;
use warnings;

# This script runs the TNM program over a set of proteins that undergo
# conformational change, whose initial and final strucure are specified in
# the file $List. The performed analysis includes the assessment of the
# anharmonicity of each normal mode, which slows down the computations and can
# be omitted by setting ANHARMONIC=0 in the file Input_TNM.in.
# The choice of degrees of freedom can be changed by setting DOF=PHI below
# (allowed: PHI, PHIPSI, OMEGA, SIDECHAIN)
# The torsional regularization can be changed setting K_OMEGA and KPHI
# to any non-negative value (KPHI: all torsions except OMEGA).
# The contact cut-off can be changed by setting CONT_THR
# Other parameters and settings can be changed in the Input_TNM.in file
#  
# For each set of parameters, you have to create a folder (f.i. named
# $DOF/KPHI$KPHI) and run the script there.
#
# These are the paths that you have to modify:
#
my $mydir="/data/ubastolla/TNM/ANHARMONIC";  # Local folder with input files
my $List=sprintf("%s/List_pairs.in", $mydir);# File with list of proteins
my $Input=sprintf("%s/Input_TNM.in", $mydir);# File with TNM parameters
my $tnm=sprintf("%s//tnm", $mydir);          # Path to TNM program
my $pdbpath="/data/ortizg/databases/pdb";    # Folder with PDB structures
my $submit="qsubmit.pl -q x86_64 -s ";       # For submitting to cluster
# my $submit="";                             # Leave empty for running locally
#
# Parameters to modify:
my $DOF="OMEGA";  # Degrees of freedom: PHI PHIPSI OMEGA SIDECHAIN
my $KPHI=0.0;     # Torsional regulariation for phi, psi, chi angles: 0.0 0.2 2
my $K_OMEGA=2.0;  # Torsional regularization for omega angle
my $CONT_THR=4.5; # Cut-off for contact definition
# end of modifications

open(my $fh, '<:encoding(UTF-8)', $List)
or die "Could not open file '$List' $!";

while (my $row = <$fh>) {
    # Reading PDB codes and chain identifiers in List
    chomp $row;
    my @res = split(/\s+/, $row); 
    if(substr($res[0],0,1) eq "#"){next;}
    my $P1; my $P2; my $C1; my $C2;

    for(my $k=1; $k<=2; $k++){
	if($k==1){
	    $P1=$res[0]; $C1=$res[1];
	    $P2=$res[2]; $C2=$res[3];
	}else{
	    $P1=$res[2]; $C1=$res[3];
	    $P2=$res[0]; $C2=$res[1];
	}

	my $pair=sprintf("%s%s-%s%s", $P1, $C1, $P2, $C2);
	print "Praparing pair ",$pair,"\n";

	# Modifying the Input file
    
	my $tmp=sprintf("tmp_%s.in", $pair);
	open(my $fo, '>', $tmp); my $out;
	$out=sprintf("PDB1= %s/%s.pdb\n", $pdbpath, $P1); print $fo $out;
	$out=sprintf("CH1= %s\n", $C1); print $fo $out;
	$out=sprintf("PDB2= %s/%s.pdb\n", $pdbpath, $P2); print $fo $out;
	$out=sprintf("CH2= %s\n", $C2); print $fo $out;

	if($DOF eq "PHI"){
	    $out=sprintf("PSI= 0\nOMEGA=0\nSIDECHAIN=0\n"); 
	}elsif($DOF eq "PHIPSI"){
	    $out=sprintf("PSI= 1\nOMEGA=0\nSIDECHAIN=0\n"); 
	}elsif($DOF eq "OMEGA"){
	    $out=sprintf("PSI= 1\nOMEGA=1\nSIDECHAIN=0\n"); 
	}elsif($DOF eq "SIDECHAIN"){
	    $out=sprintf("PSI= 1\nOMEGA=1\nSIDECHAIN=1\n");
	} 
	print $fo $out;
 
	#$out=sprintf("K_PHI= %.2f\nK_PSI= %.2f\nK_OMEGA= %.2f\nK_CHI= %.2f\nK_TOR= %.2f\n", $KPHI, $KPHI, $KPHI, $KPHI, $KPHI);
	$out=sprintf("K_PHI= %.2f\nK_PSI= %.2f\nK_CHI= %.2f\nK_TOR= %.2f\n", $KPHI, $KPHI, $KPHI, $KPHI);
	print $fo $out; 
	$out=sprintf("K_OMEGA= %.2f\n", $K_OMEGA); print $fo $out;
	$out=sprintf("CONT_THR= %.2f\n", $CONT_THR); print $fo $out;

	close $fo;

	`cat $Input >> $tmp`;
	#`mv $tmp $DIR`;

	# Running the TNM program 
	my $script=sprintf("script_tnm_%s", $pair);
	`echo $tnm $tmp > $script`;
	`chmod u+x $script`;
	#`mv $script $DIR`;

	#chdir $DIR or die "Cannot change to dir '$DIR' $!";
	`$submit $script`;
	print "Submitting ",$script,"\n";
	#chdir "../../"
    }
}
close $fh;


