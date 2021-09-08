#!/usr/bin/env perl 
#use strict;
#use warnings;

# This script reads the output of the TNM program run over a set of 
# conformational changes, with initial and final strucure specified in
# the file $List, and prints two tables with selected data for individual
# proteins and averaged over all proteins.
#
# IMPORTANT: The script assumes that TNM results are stored in a folder named
# $DOF/$KPHI where $DOF are the degrees of freedom and $KPHI is the tors.
# regularization. These variables are also used to name the output table.
#
# Output files:
# Pairs_summary_ALL_$DOF_$KPHI_$type.txt $type= harmonic, anharmonic
# (results for individual proteins)
# Pairs_average_all_ALL_$DOF_$KPHI_$type.txt
# (results averaged over all protein, for harmonic, anharmonic and difference
# anharmonic-harmonic)
#
# These are the parameters that you have to modify:
my $List="List_pairs.in";
# File with list of proteins
my $DOF="PHI";  # Degrees of freedom: PHI PHIPSI OMEGA SIDECHAIN
my $KPHI="KPHI0.2";     # Torsional regularization, e.g. 0.0 0.2 2.0
# end modifications

my $dir=sprintf("%s/%s/", $DOF, $KPHI);
my $Ref="ALL";
my $ANHAR=1;

#$name_type[0]="Fit(Lambda=0)";
#$name_type[1]="Fit(RRR_MP)";
#$name_type[2]="Fit(RRR_Cv)";
#$name_type[3]="Delta_theta";
#$name_type[4]="Iterative-RRR";
#$ntype=5;

#$name_type[0]="Fit(RRR_Cv)";

my $name_type;
$name_type[0]="harmonic";
$name_type[1]="anharmonic";
my $ntype=2;

my $name_o;
for(my $i=0; $i<40; $i++){$name_o[$i]="";}

my $name; my $i=0;
# Do not depend on normal modes nor torsional fit
$name[$i]="Nres"; $i++; # Num res
#$name[$i]="Native interactions"; $i++;
$name[$i]="Degrees of freedom";
$name_o[$i]="DOF"; $i++;
$name[$i]="RMSD_crystal"; $i++;
$name[$i]="RMSD_Cart"; $i++;
$name[$i]="RMSD_reconstruct_1"; $i++; # Warning, depends on torsional fit
my $ndata_ini=$i;

# Depend on normal modes but not on torsional fit
$name[$i]="Recp.Coll(therm)"; $i++;
$name[$i]="Ave_collectivity_NM"; $i++;
$name[$i]="RMSD_Norm.Modes"; $i++;
$name[$i]="Internal_fluct";  $name_o[$i]="Inferred_internal_fluct"; $i++;
$name[$i]="Force constant factor";  $name_o[$i]="Force_constant";$i++;
$name[$i]="r(B_pred,B_exp)"; $name_o[$i]="Bfactor_correlation"; $i++;
$name[$i]="r(therm,change)_Cart"; $name_o[$i]="Cartesian_correlation"; $i++;
#$name[$i]="RMSD_ini"; $i++;
#$name[$i]="RMSD_min"; $i++;
#$name[$i]="RMSD_std_ini(1+2)"; $i++;
#$name[$i]="RMSD_std_min(1+2)"; $i++;
$ndata=$i;

# Depend on normal modes and torsional fit
my $iZ=-1; my $iDE=-1;
#$name[$i]="RMSD_reconstruct_1"; $i++;
$name[$i]="Recp.Coll(cc)"; $i++;
$name[$i]="r(therm,change)_Tors"; $name_o[$i]="Torsional_correlation"; $i++;
$name[$i]="r[c^2,1/w^2]"; $name_o[$i]="Normal_mode_correlation"; $i++;
$name[$i]="r[(c*w)^2,1/w^2]"; $name_o[$i]="Null_mode_dev_corr"; $i++;
#$name[$i]="Recp.Coll((cc*w)^2)"; $i++;
#$name[$i]="Significance(Z)"; $iZ=$i; $i++;
$name[$i]="Delta_E0_ave/Delta_E"; $name_o[$i]="Null_mode_dev_barrier"; $i++;
$name[$i]="P(Delta_E0<Delta_E)";  $name_o[$i]="Null_mode_dev_signif"; $iDE=$i;
$i++;
my $ndata_summ=$i;
#$name_t[9]="RMSD_WTors";
#$name_t[10]="Tors_frac";

my $name1; my $value;
for($i=0; $i<$ndata_summ; $i++){
    my @res = split(/\s+/, $name[$i]);
    $name1[$i]=$res[0];
    if($name_o[$i] eq ""){$name_o[$i]=$name1[$i];}
    $value[$i]=scalar(@res);
}

$i=0; my $name_anh;
#$name_anh[$i]="1/omega"; # col.2
$name_anh[$i]="log(dKL)"; # col.4
$name_anh[$i]="log(E_anharmonic)"; # col.6
$name_anh[$i]="Coll_Cart"; # col.3
#$name_anh[$i]="max.rmsd";  # col.7
$name_anh[$i]="Corr(dKL,coll)";
$name_anh[$i]="Corr(dKL,1/omega_har^2)";
$name_anh[$i]="Corr(dKL,1/omega_anh^2)";
$name_anh[$i]="Corr(1/omega_anh^2,1/omega_har^2)";
my $ndata_anh=$i;

my $ndata_all=$ndata_summ+$ndata_anh;
my $type; my $data_1; my $data_2; my $diff_1; my $diff_2;
my $Significant_t; my $Significant_DE;
for($type=0; $type< $ntype; $type++){
    for($i=0; $i<$ndata_all; $i++){
	$data_1[$type][$i]=0; $data_2[$type][$i]=0;
    }
    $Significant_t[$type]=0;
    $Significant_DE[$type]=0;
}
for($i=0; $i<$ndata_all; $i++){
    $diff_1[$i]=0; $diff_2[$i]=0;
}
my $data_anh_1; my $data_anh_2;
for($i=0; $i<$ndata_anh; $i++){
    $data_anh_1[$i]=0; $data_anh_2[$i]=0;
}


# File with results for all pairs, direct and inverse each in a line
my $header="#1=pair";
for($i=0; $i<$ndata_summ; $i++){
    $header=sprintf("%s\t%d=%s", $header, $i+2, $name_o[$i]);
}
my $k=$i;
for($i=0; $i<$ndata_anh; $i++){
    $header=sprintf("%s\t%d=%s", $header, $k+$i, $name_anh[$i]);
}

my $Output1;
for($type=0; $type<$ntype; $type++){
    $Output1[$type]=sprintf("Pairs_summary_%s_%s_%s_%s.txt",
			    $Ref, $DOF, $KPHI, $name_type[$type]);
    open(my $fo, '>', $Output1[$type]);
    print $fo $header,"\n";
    close $fo;
}


my $npair=0;
my $nanha=0;
open(my $fh, '<:encoding(UTF-8)', $List)
or die "Could not open file '$List' $!";
while (my $row = <$fh>) {
    # Reading PDB codes and chain identifiers in List
    chomp $row;
    my @res = split(/\s+/, $row);
    if(substr($res[0],0,1) eq "#"){next;}

    my $P1; my $P2; my $C1; my $C2;
    for(my $p=1; $p<=2; $p++){
	if($p==1){
	    $P1=$res[0]; $C1=$res[1];
	    $P2=$res[2]; $C2=$res[3];
	}else{
	    $P1=$res[2]; $C1=$res[3];
	    $P2=$res[0]; $C2=$res[1];
	}
	
	my $pair=sprintf("%s%s_%s%s", $P1, $C1, $P2, $C2);
	print "Analyzing pair: ",$pair,"\n";

	my $count; my $data; my $data_anh;
	for($i=0; $i<$ndata_summ; $i++){
	    $count[$i]=0;
	    for($type=0; $type<$ntype; $type++){$data[$type][$i]=0;} 
	}
	for($i=0; $i<$ndata_anh; $i++){$data_anh[$i]=0;}

	# file summary
	#$input=sprintf("Summary_%s_MIN4.5_TNM_%s%s.dat", $pair, $Ref, $DOF);
	my $input=sprintf("%s%s.summary.dat", $dir, $pair);
	if(open($f2, '<:encoding(UTF-8)', $input)==0){
	    print "WARNING, Could not open file '$input' $!\n"; next;
	}
	my $type=-1;
	while (my $row2 = <$f2>) {
	    chomp $row2;
	    my @res2 = split(/\s+/, $row2);
	    if($res2[0] eq "#####"){next;}
	    #if($res2[0] eq "#####"){$type++; next;}	    
	    for($i=0; $i<$ndata_summ; $i++){
		my $s= substr($res2[0],0,length($name1[$i]));
		if($s eq $name1[$i]){
		    $type=$count[$i]; $count[$i]++;
		    $data[$type][$i]=$res2[$value[$i]];
		    if(($data[$type][$i] eq "nan")||
		       ($data[$type][$i] eq "-nan")){
			$data[$type][$i]=0;
		    }
		    last;
		}
	    }
	}
	close($f2);

	# Check
	my $all=1;
	for($i=$ndata_ini; $i<$ndata_summ; $i++){
	    if($count[$i]!=$ntype){
		print "ERROR, only ",$count[$i]," instance of ",
		$name[$i]," found in ",$input,"\n"; $all=0; last;
	    }
	}
	if($all==0){next;}
	$npair++;

	if(Read_anharmonic($P1.$C1)){
	    for($i=0; $i<$ndata_anh; $i++){
		$data_anh_1[$i]+=$data_anh[$i];
		$data_anh_2[$i]+=$data_anh[$i]*$data_anh[$i];
	    }
	    $nanha++;
	}

	for($type=0; $type<$ntype; $type++){
	    open(my $fo, '>>', $Output1[$type]);
	    my $out=sprintf("%s",$pair); my $y;
	    for($i=0; $i<$ndata_summ; $i++){
		if($i>=$ndata_ini){$y=$data[$type][$i];}
		else{$y=$data[0][$i];}
		$out=sprintf("%s\t%.3g", $out, $y);
		$data_1[$type][$i]+=$y;
		$data_2[$type][$i]+=$y*$y;
	    }
	    if($data[$type][$iZ]>2){$Significant_t[$type]++;}
	    if($data[$type][$iDE]<0.05){$Significant_DE[$type]++;}
	    for($i=0; $i<$ndata_anh; $i++){
		$out=sprintf("%s\t%.3g", $out, $data_anh[$i]);
	    }
	    print $fo $out,"\n";
	    close $fo;
	} # end type
	for($i=$ndata_ini; $i<$ndata_summ; $i++){
	    my $diff=$data[1][$i]-$data[0][$i]; # anhar-har
	    $diff_1[$i]+=$diff;
	    $diff_2[$i]+=$diff*$diff;
	}
    }
}
for($type=0; $type<$ntype; $type++){
    print "Writing ", $npair," conf.ch. in ",$Output1[$type],"\n";
}


$header="#1=Para"; $k=2;
for($i=$ndata_ini; $i<$ndata_summ; $i++){
    $header=sprintf("%s\t%d=%s s.e.", $header, $k, $name_o[$i]);
    $k+=2;
}
if($iZ>=0) {$header=sprintf("%s\t%d=F(t>2) s.e.", $header, $k); $k+=2;}
if($iDE>=0){
    $header=sprintf("%s\t%d=F(DE<DE0,p<0.05) s.e.", $header, $k); $k+=2;
}
for($i=0; $i<$ndata_anh; $i++){
    $header=sprintf("%s\t%d=%s s.e.", $header, $k, $name_anh[$i]);
    $k+=2;
}
my $Output=sprintf("Pairs_average_all_%s_%s_%s.txt", $Ref, $DOF, $KPHI);
open(my $fo, '>', $Output);

my $out="# ";
for($i=0; $i<$ndata_ini; $i++){
    Ave_se($data_1[0][$i], $data_2[0][$i], $npair);
    $out=sprintf("%s\t%s=\t%.4g %.2g", $out,
		 $name[$i], $data_1[0][$i], $data_2[0][$i]);
}
print $fo $out,"\n";
print $fo $header,"\n";

for($type=0; $type<$ntype; $type++){

    print $fo "# ",$name_type[$type],":\n";
    # my $out=sprintf("%s",$name_type[$type]);
    my $out=sprintf("%s_%s",$DOF, $KPHI);
    for($i=$ndata_ini; $i<$ndata_summ; $i++){
	Ave_se($data_1[$type][$i], $data_2[$type][$i], $npair);
	$out=sprintf("%s\t%.4g %.2g", $out,
		     $data_1[$type][$i], $data_2[$type][$i]);
    }
    my $p;
    if($iZ>=0){
	$p=$Significant_t[$type]/=$npair;
	$out=sprintf("%s\t%.3g %.2g",$out, $p, sqrt($p*(1-$p)/$npair));
    }
    if($iDE>=0){
	$p=$Significant_DE[$type]/=$npair;
	$out=sprintf("%s\t%.3g %.2g",$out, $p, sqrt($p*(1-$p)/$npair));
    }
    for($i=0; $i<$ndata_anh; $i++){
	if($type==0){Ave_se($data_anh_1[$i], $data_anh_2[$i], $nanha);}
	$out=sprintf("%s\t%.4g %.2g", $out,
		     $data_anh_1[$i], $data_anh_2[$i]);
    }
    print $fo $out,"\n";
}

print $fo "# anharmonic-harmonic:\n";
my $out=sprintf("%s_%s",$DOF, $KPHI);
#my $out="difference";
for($i=$ndata_ini; $i<$ndata_summ; $i++){
    Ave_se($diff_1[$i], $diff_2[$i], $npair);
    $out=sprintf("%s\t%.4g %.2g", $out, $diff_1[$i], $diff_2[$i]);
}
print $fo $out,"\n";

print "Writing ", $Output,"\n";
close($fo);

my $Output=sprintf("Pairs_average_all_line_%s_%s_%s.txt",$Ref,$DOF,$KPHI);
open(my $fo, '>', $Output);
$out="# ";
for($i=0; $i<$ndata_ini; $i++){
    $out=sprintf("%s\t%s=\t%.4g %.2g", $out,
		 $name[$i], $data_1[0][$i], $data_2[0][$i]);
}
print $fo $out,"\n";

$out=sprintf("#Var\t2=%s s.e.\t4=%s s.e.\t6=diff s.e.\n",
	     $name_type[0], $name_type[1]);
print $fo $out;
for($i=$ndata_ini; $i<$ndata_summ; $i++){
    my $out=sprintf("%s",$name1[$i]);
    for($type=0; $type<$ntype; $type++){
	$out=sprintf("%s\t%.4g %.2g", $out,
		     $data_1[$type][$i], $data_2[$type][$i]);
    }
    $out=sprintf("%s\t%.4g %.2g", $out, $diff_1[$i], $diff_2[$i]);
    print $fo $out,"\n";
}
my $p, $dp;
if($iZ>=0){
    $out="Significant_t";
    for($type=0; $type<$ntype; $type++){
	$p[$type]=$Significant_t[$type];
	$dp[$type]=$p[$type]*(1-$p[$type])/$npair;
	$out=sprintf("%s\t%.3g %.2g",$out, $p[$type], sqrt($dp[$type]));
    }
    $out=sprintf("%s\t%.3g %.2g",$out, $p[1]-$p[0], sqrt($dp[1]+$dp[2]));
    print $fo $out,"\n";
}
if($iDE>=0){
    $out="Significant_DE";
    for($type=0; $type<$ntype; $type++){
	$p[$type]=$Significant_DE[$type];
	$dp[$type]=$p[$type]*(1-$p[$type])/$npair;
	$out=sprintf("%s\t%.3g %.2g",$out, $p[$type], sqrt($dp[$type]));
    }
    $out=sprintf("%s\t%.3g %.2g",$out, $p[1]-$p[0], sqrt($dp[1]+$dp[2]));
    print $fo $out,"\n";
}

for($i=0; $i<$ndata_anh; $i++){
    my $out=sprintf("%s",$name_anh[$i]);
    $out=sprintf("%s\t%.4g %.2g", $out,
		 $data_anh_1[$i], $data_anh_2[$i]);
    print $fo $out,"\n";
}

print "Writing ", $Output,"\n";
close($fo);

sub Ave_se{ #($sum_1, $sum_2, $n)
    my $n=$_[2];
    if($n>1){
	$_[0]/=$n;
	$_[1]=($_[1]/$n-$_[0]*$_[0])/($n-1);
	if($_[1] > 0.0){
	    $_[1]=sqrt($_[1]);
	}else{
	    print "WARNING, negative variance: ", $_[1],"\n";
	    $_[1]=0;
	}
    }
}

sub Read_anharmonic{ # $_[0]=$pair

    for(my $i=0; $i<7; $i++){$data_anh[$i]=0;}

    my $input=
	#sprintf("%s/%s%s_MIN4.5_TNM_%s%s_anharmonic.dat",
	#	    $dir[$k], $P1, $C1, $Ref[$k], $DOF[$k]);
	sprintf("%s%s.anharmonic.dat", $dir, $_[0]);
    if(! -s $input){
	print "WARNING, could not open file ",$input,"\n"; return 0;
    }
    open(my $f2, '<:encoding(UTF-8)', $input);

    my $n=0;
    my $wh_w=0;  my $wa_w=0; my $coll_w=0; my $E_w=0; my $max_rmsd_w=0;
    my $dKL_w=0; my $dKL_coll=0; my $dKL_wh=0; my $dKL_wa=0; my $wa_wh=0;
    # my $dKL_1=0; my $coll_1=0; my $wh_1=0; my $wa_1=0;
    my $dKL_2=0; my $coll_2=0; my $wh_2=0; my $wa_2=0;
    my $w_sum=0;
    while (my $row2 = <$f2>) {
	if(substr($row2,0,1) eq "#"){next;}
	chomp $row2;
	my @res2 = split(/\s+/, $row2);
	if(($res2[4] eq "inf")||($res2[4]<=0)){next;} # dKL
	my $dKL=log($res2[4]);
	my $w_har=$res2[1]*$res2[1]; #1/omega^2
	my $w_anhar=$res2[2]*$res2[2];
	my $coll=$res2[3];
	my $E=$res2[6]; if($E>0.05){$E=log($E);}else{$E=log(0.05);}
	my $max_rmsd=$res2[7];
	my $w=$w_har; # Weight with harmonic frequencies
	$w_sum+=$w;
	$wh_w+=$w*$w_har;
	$wa_w+=$w*$w_anhar;
	$E_w+=$w*$E;
	$dKL_w+=$w*$dKL;
	$coll_w+=$w*$coll;
	$max_rmsd_w+=$w*$max_rmsd;
	#$dKL_coll+=$dKL*$coll;
	#$dKL_1+=$dKL; $dKL_2+=$dKL*$dKL;
	#$coll_1+=$coll; $coll_2+=$coll*$coll;
	#$dKL_wh+=$dKL*$w_har;
	#$wh_1+=$w_har; $wh_2+=$w_har*$w_har;
	#$dKL_wa+=$dKL*$w_anhar;
	#$wa_1+=$w_anhar; $wa_2+=$w_anhar*$w_anhar;
	#$wa_wh+=$w_har*$w_anhar;
	$dKL_coll+=$w*$dKL*$coll;
	$dKL_2+=$w*$dKL*$dKL;
	$coll_2+=$w*$coll*$coll;
	$dKL_wh+=$w*$dKL*$w_har;
	$wh_2+=$w*$w_har*$w_har;
	$dKL_wa+=$w*$dKL*$w_anhar;
	$wa_2+=$w*$w_anhar*$w_anhar;
	$wa_wh+=$w*$w_har*$w_anhar;
	$n++;
    }
    close $f2;
    if($n<1){
	print "ERROR in ",$input," modes=",$n,"\n";
	return 0;
    }

    #$data_anh[0]=$w_sum/$norm;
    $norm=$w_sum;
    $E_w/=$norm; 
    $dKL_w/=$norm;
    $wa_w/=$norm;
    $wh_w/=$norm;
    $coll_w/=$norm;
    $data_anh[0]=$dKL_w;
    $data_anh[1]=$E_w;
    $data_anh[2]=$coll_w;
    #$data_anh[4]=$max_rmsd_w/$norm;
    #$dKL_1/=$n; $dKL_2-= $n*$dKL_1*$dKL_1;
    #$coll_1/=$n; $coll_2-= $n*$coll_1*$coll_1;
    #$wa_1/=$n; $wa_2-= $n*$wa_1*$wa_1;
    #$wh_1/=$n; $wh_2-= $n*$wh_1*$wh_1;
    #
    $dKL_2= $dKL_2/$norm-$dKL_w*$dKL_w;
    $coll_2= $coll_2/$norm- $coll_w*$coll_w;
    $wa_2= $wa_2/$norm-$wa_w*$wa_w;
    $wh_2= $wh_2/$norm-$wh_w*$wh_w;
    if(($dKL_2<=0)||($coll_2<=0)||($wa_2<=0)||($wh_2<=0)){
	print "WARNING in ",$input," Var(dKL)=0\n"; # die;
	return 0;
    }else{
	# Compute correlation coefficients
	$data_anh[3]=($dKL_coll/$norm-$dKL_w*$coll_w)/sqrt($dKL_2*$coll_2);
	$data_anh[4]=($dKL_wh/$norm-$dKL_w*$wh_w)/sqrt($dKL_2*$wh_2);
	$data_anh[5]=($dKL_wa/$norm-$dKL_w*$wa_w)/sqrt($dKL_2*$wa_2);
	$data_anh[6]=($wa_wh/$norm-$wa_w*$wh_w)/sqrt($wa_2*$wh_2);
    }
    return 1;
}
