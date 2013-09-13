#!/usr/bin/perl
use warnings;

@ARGV == 2 || die("Usage: $0 <directory> <result>");

$DIR = $ARGV[0];
$RESULT = $ARGV[1];

@probs = (0.25, 0.5, 0.75);
@correct_alleles = (4, 3, 2, 1);
@parents = (5, 10, 15, 20);
@methods = (
#	"./csodp -m 7",
#	"./csodp --heuristic 1 -n 50", "./csodp --heuristic 1 -n 100", "./csodp --heuristic 1 -n 150",
#	"./csodp --heuristic 2 -n 50", "./csodp --heuristic 2 -n 100", "./csodp --heuristic 2 -n 150",
#	"./csodp --heuristic 3 -n 50", "./csodp --heuristic 3 -n 100", "./csodp --heuristic 3 -n 150",
#	"./csodp --heuristic 4 -n 50", "./csodp --heuristic 4 -n 100", "./csodp --heuristic 4 -n 150",
#	"./csodp --heuristic 5 -n 50", "./csodp --heuristic 5 -n 100", "./csodp --heuristic 5 -n 150",
#	"./csodp --heuristic 1 -n 10 -l 5", "./csodp --heuristic 1 -n 10 -l 10",
#	"./csodp --heuristic 1 -n 50 -l 5", "./csodp --heuristic 1 -n 100 -l 5", "./csodp --heuristic 1 -n 150 -l 5",
#	"./csodp --heuristic 2 -n 10 -l 5", "./csodp --heuristic 2 -n 10 -l 10",
#	"./csodp --heuristic 2 -n 50 -l 5", "./csodp --heuristic 2 -n 100 -l 5", "./csodp --heuristic 2 -n 150 -l 5",
#	"./csodp --heuristic 4 -n 10 -l 5", "./csodp --heuristic 4 -n 10 -l 10",
#	"./csodp --heuristic 4 -n 50 -l 5", "./csodp --heuristic 4 -n 100 -l 5", "./csodp --heuristic 4 -n 150 -l 5",
#	"./csodp --heuristic 5 -n 10 -l 5", "./csodp --heuristic 5 -n 10 -l 10",
#	"./csodp --heuristic 5 -n 50 -l 5", "./csodp --heuristic 5 -n 100 -l 5", "./csodp --heuristic 5 -n 150 -l 5"
#	"./csora -m 500000", "./csoga"
#	"./csodp -g -c 2 -m 7 --heuristic 1 -n 50", "./csodp -g -c 2 -m 7 --heuristic 1 -n 100", "./csodp -g -c 2 -m 7 --heuristic 1 -n 150",
#	"./csodp -g -c 2 -m 7 --heuristic 2 -n 50", "./csodp -g -c 2 -m 7 --heuristic 2 -n 100", "./csodp -g -c 2 -m 7 --heuristic 2 -n 150",
#	"./csodp -g -c 2 -m 7 --heuristic 3 -n 50", "./csodp -g -c 2 -m 7 --heuristic 3 -n 100", "./csodp -g -c 2 -m 7 --heuristic 3 -n 150",
#	"./csodp -g -c 2 -m 7 --heuristic 4 -n 50", "./csodp -g -c 2 -m 7 --heuristic 4 -n 100", "./csodp -g -c 2 -m 7 --heuristic 4 -n 150",
#	"./csodp -g -c 2 -m 7 --heuristic 5 -n 50", "./csodp -g -c 2 -m 7 --heuristic 5 -n 100", "./csodp -g -c 2 -m 7 --heuristic 5 -n 150",
	"./csodp -g -c 2 -m 7 --heuristic 1 -n 10 -l 5", "./csodp -g -c 2 -m 7 --heuristic 1 -n 10 -l 10",
	"./csodp -g -c 2 -m 7 --heuristic 1 -n 50 -l 5", "./csodp -g -c 2 -m 7 --heuristic 1 -n 100 -l 5", "./csodp -g -c 2 -m 7 --heuristic 1 -n 150 -l 5",
	"./csodp -g -c 2 -m 7 --heuristic 2 -n 10 -l 5", "./csodp -g -c 2 -m 7 --heuristic 2 -n 10 -l 10",
	"./csodp -g -c 2 -m 7 --heuristic 2 -n 50 -l 5", "./csodp -g -c 2 -m 7 --heuristic 2 -n 100 -l 5", "./csodp -g -c 2 -m 7 --heuristic 2 -n 150 -l 5",
	"./csodp -g -c 2 -m 7 --heuristic 3 -n 10 -l 5", "./csodp -g -c 2 -m 7 --heuristic 3 -n 10 -l 10",
	"./csodp -g -c 2 -m 7 --heuristic 3 -n 50 -l 5", "./csodp -g -c 2 -m 7 --heuristic 3 -n 100 -l 5", "./csodp -g -c 2 -m 7 --heuristic 3 -n 150 -l 5",
	"./csodp -g -c 2 -m 7 --heuristic 4 -n 10 -l 5", "./csodp -g -c 2 -m 7 --heuristic 4 -n 10 -l 10",
	"./csodp -g -c 2 -m 7 --heuristic 4 -n 50 -l 5", "./csodp -g -c 2 -m 7 --heuristic 4 -n 100 -l 5", "./csodp -g -c 2 -m 7 --heuristic 4 -n 150 -l 5",
	"./csodp -g -c 2 -m 7 --heuristic 5 -n 10 -l 5", "./csodp -g -c 2 -m 7 --heuristic 5 -n 10 -l 10",
	"./csodp -g -c 2 -m 7 --heuristic 5 -n 50 -l 5", "./csodp -g -c 2 -m 7 --heuristic 5 -n 100 -l 5", "./csodp -g -c 2 -m 7 --heuristic 5 -n 150 -l 5",
#	"./csora -m 500000", "./csoga"
	"./csodp -g -c 2 -m 7 --heuristic 10 -n 10", "./csodp -g -c 2 -m 7 --heuristic 10 -n 20", "./csodp -g -c 2 -m 7 --heuristic 10 -n 30",
	"./csodp -g -c 2 -m 7 --heuristic 10 -n 40", "./csodp -g -c 2 -m 7 --heuristic 10 -n 50",
	"./csodp -g -c 2 -m 7 --heuristic 10 -n 60", "./csodp -g -c 2 -m 7 --heuristic 10 -n 70", "./csodp -g -c 2 -m 7 --heuristic 10 -n 80",
	"./csodp -g -c 2 -m 7 --heuristic 10 -n 90", "./csodp -g -c 2 -m 7 --heuristic 10 -n 100",
	"./csodp -g -c 2 -m 7 --heuristic 11 -n 10", "./csodp -g -c 2 -m 7 --heuristic 11 -n 20", "./csodp -g -c 2 -m 7 --heuristic 11 -n 30",
	"./csodp -g -c 2 -m 7 --heuristic 11 -n 40", "./csodp -g -c 2 -m 7 --heuristic 11 -n 50",
	"./csodp -g -c 2 -m 7 --heuristic 11 -n 60", "./csodp -g -c 2 -m 7 --heuristic 11 -n 70", "./csodp -g -c 2 -m 7 --heuristic 11 -n 80",
	"./csodp -g -c 2 -m 7 --heuristic 11 -n 90", "./csodp -g -c 2 -m 7 --heuristic 11 -n 100",
	"./csodp -g -c 2 -m 7 --heuristic 12 -n 10", "./csodp -g -c 2 -m 7 --heuristic 12 -n 20", "./csodp -g -c 2 -m 7 --heuristic 12 -n 30",
	"./csodp -g -c 2 -m 7 --heuristic 12 -n 40", "./csodp -g -c 2 -m 7 --heuristic 12 -n 50",
	"./csodp -g -c 2 -m 7 --heuristic 12 -n 60", "./csodp -g -c 2 -m 7 --heuristic 12 -n 70", "./csodp -g -c 2 -m 7 --heuristic 12 -n 80",
	"./csodp -g -c 2 -m 7 --heuristic 12 -n 90", "./csodp -g -c 2 -m 7 --heuristic 12 -n 100",
	);
@methodnames = (
#	"dp",
#	"dp1-g-50", "dp1-g-100", "dp1-g-150",
#	"dp2-g-50", "dp2-g-100", "dp2-g-150",
#	"dp3-g-50", "dp3-g-100", "dp3-g-150",
#	"dp4-g-50", "dp4-g-100", "dp4-g-150",
#	"dp5-g-50", "dp5-g-100", "dp5-g-150",
	"dp1-g-10-5", "dp1-g-10-10", "dp1-g-50-5", "dp1-g-100-5", "dp1-g-150-5",
	"dp2-g-10-5", "dp2-g-10-10", "dp2-g-50-5", "dp2-g-100-5", "dp2-g-150-5",
	"dp3-g-10-5", "dp3-g-10-10", "dp3-g-50-5", "dp3-g-100-5", "dp3-g-150-5",
	"dp4-g-10-5", "dp4-g-10-10", "dp4-g-50-5", "dp4-g-100-5", "dp4-g-150-5",
	"dp5-g-10-5", "dp5-g-10-10", "dp5-g-50-5", "dp5-g-100-5", "dp5-g-150-5",
#	"ra", "ga");
	"dp10-g-10", "dp10-g-20", "dp10-g-30", "dp10-g-40", "dp10-g-50",
	"dp10-g-60", "dp10-g-70", "dp10-g-80", "dp10-g-90", "dp10-g-100",
	"dp11-g-10", "dp11-g-20", "dp11-g-30", "dp11-g-40", "dp11-g-50",
	"dp11-g-60", "dp11-g-70", "dp11-g-80", "dp11-g-90", "dp11-g-100",
	"dp12-g-10", "dp12-g-20", "dp12-g-30", "dp12-g-40", "dp12-g-50",
	"dp12-g-60", "dp12-g-70", "dp12-g-80", "dp12-g-90", "dp12-g-100",
	);
$number_of_inputs = 20;

foreach my $prob (@probs)
{
	foreach my $allelecount (@correct_alleles)
	{
		foreach my $parentcount (@parents)
		{
			my $dir = "heterozygous-l6-p$parentcount-max$allelecount-hp$prob";

			for (my $i = 1; $i <= $number_of_inputs; $i++)
			{
				for (my $j = 0; $j < @methods; $j++)
				{
					print "Method $methodnames[$j]. Running $DIR/$dir/$i/input$i.xml...";
					system("$methods[$j] $DIR/$dir/$i/input$i.xml > $DIR/$dir/$i/$methodnames[$j].dot 2>> $RESULT");
					system("dot -Tpng $DIR/$dir/$i/$methodnames[$j].dot -o $DIR/$dir/$i/$methodnames[$j].png");
					print " Done!\n";
				}
			}
		}
	}
}
