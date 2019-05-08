#!/usr/bin/perl

my $snpFile = $ARGV[0];

open (SNP, "$snpFile") || die "cant load file";


$| = 1;
while (<SNP>)
{
	chomp;
	@input = split(" ", $_);	
	$scaffold = $input[0];
	$position = $input[1];
	$ref = $input[2];
	@maArray = ();
	for ($x=3;$x<(scalar(@input)-1);$x++)
	{
		push(@maArray, $input[$x]);
	}
	$cons = $input[$x];

	$line = 0;
	$diff = 0;
	foreach $element (@maArray)
	{
		if ($element =~ /-/)
		{
			next;
		}
		if ($element ne $cons)
		{
			$diff++;	
		}
		$line++;
	}	
	if (($line > 2) and ($diff < $line) and ($diff > 0) and ($diff < 2))
	{
	    print $_, "\n";
	}
}
