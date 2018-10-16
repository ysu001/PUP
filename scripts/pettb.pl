#!/usr/bin/perl
#Includes
use strict;
use warnings;
use IO::File;
use Getopt::Long;
use POSIX;
use List::Util qw[min max];

my %options = (
		infile => undef,
		constant => undef,
               	output => undef
               );
&GetOptions(
            'infile=s'   => \$options{input},
	    'constant=s' => \$options{constant},
            'output=s' => \$options{output},
            'output2=s' => \$options{output2},
            );

#Usage
my $usage = "

USAGE: forYI.pl -infile <infile> -constant <constant> -output <output> -output2 <output2>

";


die $usage unless ((-s  $options{input}) and (defined $options{constant}) and (defined $options{output}) and (defined $options{output2}));


my $file = new IO::File $options{input};
my $constant = $options{constant};
my $outfile = new IO::File ">$options{output}";
my $outfile2 = new IO::File ">$options{output2}";
my @array1; my @array2; my @array3; my @array4; my @array5; my %indexhash;
my $i = 0;

while (my $line =<$file>) {
	chomp($line);
	my @array = split(/\s+/, $line);
	push(@array1, $array[2]);
	push(@array2, $array[4]);
 
}

my $sum;
my $total = @array1;
@array3 = ($array1[0]);
while ($i < $total-1) {
	my $j = $i+1;
	$sum = $array3[$i++] + $array1[$j];
	push(@array3, $sum);
}

my $start = 0;
foreach my $element (@array3){
	my $ceil = $element/$constant;
	$ceil= ceil($ceil);
	if (($ceil - $start) > 1){
		$ceil = $start + 1;
	}
	push(@array4, $ceil);
	print $outfile2 "$ceil\n";
	$start = $ceil;
}

$i = 0;
foreach my $element (@array4) {
	push (@{$indexhash{$element}}, $array2[$i++]);
}

foreach my $key (sort {$a <=> $b} keys %indexhash) {
	my $min = min(@{$indexhash{$key}});
	my $max = max(@{$indexhash{$key}});
	print $outfile "$min\t$max\n";
}

