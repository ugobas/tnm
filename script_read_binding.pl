#!/usr/bin/env perl 
use strict;

my $List="antibody_dimers.csv";
my $Input="Input_TNM.in";
my $tnm="/data/ortizg/software/Linux/ubastolla//tnm";
my $dir="/data/ubastolla/TNM/BINDING";
my $pdbdir="/data/ortizg/databases/pdb/";
my $Output="Binding_energies.txt";
my $Fits="Binding_energies_fits.txt";

open(my $fh, '<:encoding(UTF-8)', $List)
or die "Could not open file '$List' $!";
open(my $fo, '>', $Output);
print $fo "### name\tDG_exp1\tdof\tDG_cont\tDG_NM\tDSASA\tG_NM\tCoordination\tDirectionality\tDeformation\tMin_Coord\tMin_Dir\tMax_Def\tChain_Coord\tChain_Dir\tChain_Def\n";
print $fo "### 0 0 0 1 1 1 0 1 1 1 1 1 1 1 1 1\n";
#\tDG_exp2

my @dof; my @DG_NM; my @DG_cont; my @DG_exp1; my @DG_exp2; my @name;
my @DG_NM_dof; my @G_NM; my @G_NM_apo; my @D_SASA;
my @Coordination; my @Directionality; my @Deformation;
my @Min_Coord; my @Min_Dir; my @Max_Def;
my @Chain_Coord; my @Chain_Dir; my @Chain_Def;
my @x; my @y; my $n=0;

while (my $row = <$fh>) {
    # Reading PDB codes and chain identifiers in List
    chomp $row;
    if(substr($row,0,1) eq "#"){next;}
    my @res = split(/;/, $row);
    my $P1=$res[0]; my $C1=$res[1]; my $C2=$res[2]; my $C3=$res[3];
    $DG_exp1[$n]=$res[6]; $DG_exp2[$n]=$res[14]; 

    print "Analyzed protein: ",$P1," ",$C1,$C2,$C3," ",$n+1,"\n";

    $name[$n]=sprintf("%s%s%s%s", $P1,$C1,$C2,$C3);

    my $tmp=`ls $P1*.summary.dat`;
    chomp $tmp;
    unless(-e $tmp){
	print "Could not open file '$tmp'\n";
	print $fo "# ", $name[$n],"\n";
	next;
    }
    open(my $f2, '<:encoding(UTF-8)', $tmp);

    while (my $row = <$f2>) {
	# Reading results in summary
	if(substr($row,0,1) eq "#"){next;}
	chomp $row;
	@res = split(/\s+/, $row);
	if($res[0] eq "Degrees"){
	    $dof[$n]= $res[3];
	}elsif($res[0] eq "Binding"){
	    if($res[3] eq "NM"){$DG_NM[$n]=$res[4];}
	    elsif($res[3] eq "cont"){$DG_cont[$n]=$res[4];}
	}elsif($res[0] eq "Free"){
	    if($res[4] eq "apo"){$G_NM_apo[$n]= $res[5];}
	    else{$G_NM[$n]= $res[4];}
	}elsif($res[0] eq "SASA"){
	    $D_SASA[$n]= $res[2];
	}
    }
    close $f2;
    $DG_NM_dof[$n]=$DG_NM[$n]; if($dof[$n]){$DG_NM_dof[$n]/=$dof[$n];}

    $tmp=`ls $P1*_interface.dat`;
    chomp $tmp;
    unless(-e $tmp){
	print "Could not open file '$tmp'\n";
	print $fo "# ", $name[$n],"\n";
	next;
    }
    open(my $f2, '<:encoding(UTF-8)', $tmp);

    while (my $row = <$f2>) {
	# Reading results in interface
	if(substr($row,0,1) eq "#"){next;}
	chomp $row;
	@res = split(/\s+/, $row);
	if($res[0] eq "Directionality_coupling:"){
	    $Directionality[$n]=$res[2]; $Min_Dir[$n]=$res[4];
	}elsif($res[0] eq "Coordination_coupling:"){
	    $Coordination[$n]=$res[2]; $Min_Coord[$n]=$res[4];
	}elsif($res[0] eq "Deformation_coupling:"){
	    $Deformation[$n]=$res[2]; $Max_Def[$n]=$res[4];
	}
    }
    close $f2;

    $tmp=`ls $P1*_couplings.dat`;
    chomp $tmp;
    unless(-e $tmp){
	print "Could not open file '$tmp'\n";
	print $fo "# ", $name[$n],"\n";
	next;
    }
    open(my $f2, '<:encoding(UTF-8)', $tmp);

    my $m=0;
    $Chain_Coord[$n]=0; $Chain_Dir[$n]=0; $Chain_Def[$n]=0; 
    while (my $row = <$f2>) {
	# Reading results in couplings
	if(substr($row,0,1) eq "#"){next;}
	chomp $row;
	@res = split(/\s+/, $row);
	if((($res[0] eq $C3)&&($res[1] ne $C3))|
	   (($res[1] eq $C3)&&($res[0] ne $C3))){
	    $m++;
	    $Chain_Dir[$n]+=$res[2];
	    $Chain_Coord[$n]+=$res[3];
	    $Chain_Def[$n]+=$res[4];
	}
    }
    close $f2;

    my $out=sprintf("%s\t%.4g\t%d\t%.4g\t%.4g\t%.4g\t%.4g",
		    $name[$n],$DG_exp1[$n],$dof[$n],$DG_cont[$n],
		    $DG_NM[$n],
		    $D_SASA[$n],$G_NM[$n]);
    $out=sprintf("%s\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g",
		 $out,$Coordination[$n],$Directionality[$n],$Deformation[$n],
		 $Min_Coord[$n],$Min_Dir[$n],$Max_Def[$n]); 
    $out=sprintf("%s\t%.4g\t%.4g\t%.4g\n", $out,
		 $Chain_Coord[$n],$Chain_Dir[$n],$Chain_Def[$n]);

    #,$DG_exp2[$n]
    print $fo $out;
    $n++;

}
close $fo;
print "Writing ",$Output,"\n";

open(my $fo, '>', $Fits); 
print $fo
    "# a Corr(DG_exp,a*DG_NM+DG_cont) Corr(DG_exp,a*D_SASA+DG_cont) Corr(DG_exp,a*D_SASA+DG_NM)\n";

for(my $i=0; $i<$n; $i++){$y[$i]=$DG_exp1[$i];}
for(my $a=-100; $a<100; $a+=0.1){
    for(my $i=0; $i<$n; $i++){
	#$x[$i]=$a*$DG_NM[$i]+$DG_cont[$i];
	$x[$i]=$a*$G_NM[$i]+$DG_cont[$i];
    }
    my $r1=Corr_coeff();

    for(my $i=0; $i<$n; $i++){
	$x[$i]=$a*$D_SASA[$i]+$DG_cont[$i];
    }
    my $r2=Corr_coeff();

    for(my $i=0; $i<$n; $i++){
	$x[$i]=$a*$D_SASA[$i]+$G_NM[$i];
    }
    my $r3=Corr_coeff();

    my $out=sprintf("%.3g\t%.4g\t%.4g\t%.4g\n", $a, $r1, $r2, $r3);
    print $fo $out;
}
close $fo;

print "Writing ",$Fits,"\n";

sub Corr_coeff{
    my $y1=0; my $y2=0; my $x1=0; my $x2=0; my $xy=0;
    for(my $i=0; $i<$n; $i++){
	$x1+=$x[$i]; $x2+=$x[$i]*$x[$i];
	$y1+=$y[$i]; $y2+=$y[$i]*$y[$i];
	$xy+=$x[$i]*$y[$i];
    }
    my $r=($n*$xy - $x1*$y1)/
	sqrt(($n*$x2 - $x1*$x1)*($n*$y2 - $y1*$y1));
    return $r;
}
