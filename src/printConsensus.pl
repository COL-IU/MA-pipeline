#!/usr/bin/perl

my $snpFile = $ARGV[0];
my $consensusFile = $ARGV[1];
my $MIN_NUM_READS=$ARGV[2]-1;

open (SNP, "$snpFile") || die "cant load file";
open (DEST, ">$consensusFile") || die "cant load file";

$| = 1;
while (<SNP>)
{
	chomp;
	@input = split(" ", $_);	
	$scaffold = $input[0];
	$position = $input[1];
	$ref = $input[2];
	print DEST $scaffold;
	print DEST " ";
	print DEST $position;
	print DEST " ";
	print DEST $ref;
	$MAnum = 1;
	$Acons = 0;
	$Ccons = 0;
	$Gcons = 0;
	$Tcons = 0;
	$allTotal = 0;
	$totalLength = scalar(@input);
	for ($x=3;$x<scalar(@input);$x=$x+8)
	{
		$MAtotal = 0;
		$A{$MAnum}=$input[$x];
		$a{$MAnum}=$input[$x+1];
                $C{$MAnum}=$input[$x+2];
                $c{$MAnum}=$input[$x+3];
                $G{$MAnum}=$input[$x+4];
                $g{$MAnum}=$input[$x+5];
                $T{$MAnum}=$input[$x+6];
                $t{$MAnum}=$input[$x+7];
		
		$Atotal = $A{$MAnum}+$a{$MAnum};
                $Ctotal = $C{$MAnum}+$c{$MAnum};
                $Gtotal = $G{$MAnum}+$g{$MAnum};
                $Ttotal = $T{$MAnum}+$t{$MAnum};

		$MAtotal += $Atotal+$Ctotal+$Gtotal+$Ttotal;
		$allTotal += $Atotal+$Ctotal+$Gtotal+$Ttotal;	

		$Acons += $Atotal;
                $Ccons += $Ctotal;
                $Gcons += $Gtotal;
                $Tcons += $Ttotal;
	
		if ($MAtotal eq 0)
		{
			print DEST " -";
			$MAnum++;
			next;
		}

		if (($Atotal/$MAtotal)>.80)
		{
			if (($A{$MAnum}>$MIN_NUM_READS) and ($a{$MAnum}>$MIN_NUM_READS))
			{
				print DEST " A";
				$MAnum++;
				next;
			}
		}
                if (($Ctotal/$MAtotal)>.80)
                {
                        if (($C{$MAnum}>$MIN_NUM_READS) and ($c{$MAnum}>$MIN_NUM_READS))
                        {
                                print DEST " C";
				$MAnum++;
                                next;
                        }
                }
                if (($Gtotal/$MAtotal)>.80)
                {
                        if (($G{$MAnum}>$MIN_NUM_READS) and ($g{$MAnum}>$MIN_NUM_READS))
                        {
                                print DEST " G";
				$MAnum++;
                                next;
                        }
                }
                if (($Ttotal/$MAtotal)>.80)
                {
                        if (($T{$MAnum}>$MIN_NUM_READS) and ($t{$MAnum}>$MIN_NUM_READS))
                        {
                                print DEST " T";
				$MAnum++;
                                next;
                        }
                }
		print DEST " -";	
		$MAnum++;
		next;
	}
	if ($allTotal eq 0)
        {
		print DEST "\n";
        	next;
        }
	if (($Acons/$allTotal)>.50)
        {
        	print DEST " A";
	}
        if (($Ccons/$allTotal)>.50)
        {
                print DEST " C";
        }
        if (($Gcons/$allTotal)>.50)
        {
                print DEST " G";
        }
        if (($Tcons/$allTotal)>.50)
        {
                print DEST " T";
        }
	print DEST "\n";
}
