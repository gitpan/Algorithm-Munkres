package Algorithm::Munkres;

use 5.008005;
use strict;
use warnings;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

our @EXPORT = qw( assign );

our $VERSION = '0.01';

#Variables global to the package
my @mat = ();
my @mask = ();
my @colcov = ();
my @rowcov = ();
my $Z0_row = 0;
my $Z0_col = 0;
my @path = ();
#my @org_mat = ();

#The exported subroutine.
#Expected Input: Reference to the input matrix (MxN)
#Output: Mx1 matrix, giving the column number of the value assigned to each row. (For more explaination refer perldoc)
sub assign
{
    #reference to the input matrix
    my $rmat = shift;
    my $rsolution_mat = shift;
    
    #variables local to the subroutine
    my $step = 0;
    my ($i, $j) = (0,0);
#    my @out_mat = ();

    #the input matrix
    @mat = @$rmat;

    #copy the orginal matrix, before applying the algorithm to the matrix
#   @org_mat = @mat;

    #initialize mask, column cover and row cover matrices
    @colcov = ( 0 ) x $#mat;
    @rowcov = ( 0 ) x $#mat;
    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    $mask[$i][$j] = 0;
	}
    }

    #The algorithm can be grouped in 6 steps.
    &stepone();
    &steptwo();
    $step = &stepthree();
    while($step == 4)
    {
	$step = &stepfour();
	while($step == 6)
	{
	    &stepsix();
	    $step = &stepfour();	    
	}
	&stepfive();
	$step = &stepthree();
    }

    #create the output matrix
    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($mask[$i][$j] == 1)
	    {
		$rsolution_mat->[$i] = $j;
	    }
	}
    }

#Code for tracing------------------
    print "\nInput Matrix:\n";
    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    print $mat[$i][$j] . "\t";
	}
	print "\n";
    }
	print "\nMask Matrix:\n";
    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    print $mask[$i][$j] . "\t";
	}
	print "\n";
    }
	print "\nOutput Matrix:\n";
    for($i=0;$i<=$#mat;$i++)
    {
	print $rsolution_mat->[$i] . "\n";
    }
#----------------------------------

}

#Step 1 - Find minimum value for every row and subtract this min from each element of the row.
sub stepone
{
    print "Step 1 \n";
    my ($i, $j, @min);

    #Find the minimum value for every row
    for($i=0;$i<=$#mat;$i++)
    {
	$min[$i] = $mat[$i][0];
	for($j=0;$j<=$#mat;$j++)
	{
	    if($min[$i] > $mat[$i][$j])
	    {
		$min[$i] = $mat[$i][$j];
	    }
	    print $mat[$i][$j] . "\t";
	}    
	
        #Subtract the minimum value of the row from each element of the row.
	for($j=0;$j<=$#mat;$j++)
	{
	    $mat[$i][$j] -= $min[$i];
	}	
    }
    print "Step 1 end \n";
}

#Step 2 - Star the zeroes, Create the mask and cover matrices. Re-initialize the cover matrices for next steps.
#To star a zero: We search for a zero in the matrix and than cover the column and row in which it occurs. Now this zero is starred.
#A next starred zero can occur only in those columns and rows which have not been previously covered by any other starred zero.
sub steptwo
{
    print "Step 2 \n";
 
    my ($i, $j) = (0,0);

    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($mat[$i][$j] == 0 && $colcov[$j] == 0 && $rowcov[$j] == 0)
	    {
		$mask[$i][$j] = 1;
		$colcov[$j] = 1;
		$rowcov[$i] = 1;
	    }
	}
    }
    #Re-initialize the cover matrices
    &clear_covers();
    print "Step 2 end\n";
}


#Step 3 - Check if each column has a starred zero. If yes then the problem is solved else proceed to step 4
sub stepthree
{
    print "Step 3 \n";

    my $cnt = 0;
    my ($i, $j) = (0,0);

    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($mask[$i][$j] == 1)
	    {
		$colcov[$j] = 1;
		$cnt++;
	    }
	}
    }
    if($cnt > $#mat)
    {
       print "Step 3 end. Next expected step 7 \n";
       return 7;
    }
    else
    {
       print "Step 3 end. Next expected step 4 \n";
       return 4;
    }

}

#Step 4 - Try to find a zero which is not starred and whose columns and rows are not yet covered. 
#If such a zero found, prime it, try to find a starred zero in its row, 
#                                                 if not found proceed to step 5 
#                                                 else continue
#Else proceed to step 6.
sub stepfour
{
    print "Step 4 \n";
    my $done  = 0;
    my $step = 0;
    my ($row,$col) = (-1,-1);
    my $star_col = -1;

    while($done == 0)
    {
	($row, $col) = &find_a_zero();
	if($row < 0)
	{
	    $done = 1;
	    $step = 6;
	}
	else
	{
	    $mask[$row][$col] = 2;
	    $star_col = &find_star_in_row($row);
	    if($star_col >= 0)
	    {
		$col = $star_col;
		$rowcov[$row] = 1;
		$colcov[$col] = 0;
	    }
	    else
	    {
		$done = 1;
		$Z0_row = $row;
		$Z0_col = $col;
		$step = 5;;
	    }
	}
    }
    print "Step 4 end. Next expected step $step\n";
    return $step;
}

#Tries to find yet uncovered zero
sub find_a_zero
{
    my ($row, $col);
    my ($i, $j);
    my $done = 0;
    $row = -1;
    $col = -1;
    $i = 0;

    do
    {
	$j = 0;
	do
	{
	    if($mat[$i][$j] == 0 && $rowcov[$i] == 0 and $colcov[$j] == 0)
	    {
		$row = $i;
		$col = $j;
		$done = 1;
	    }
	    $j++;
	}while($j<=$#mat);
	$i++;
    }while($i<=$#mat && $done == 0);

    return ($row, $col);
}

#Tries to find starred zero in the given row and returns the column number
sub find_star_in_row
{
    my $row = shift;
    my $col = -1;
    my $j = 0;

    for($j=0;$j<=$#mat;$j++)
    {
	if($mask[$row][$j] == 1)
	{
	    $col = $j;
	}
    }
    return $col;
}

#Step 5 - Try to find a starred zero in the column of the uncovered zero found in the step 4.
#If starred zero found, try to find a prime zero in its row.
#Continue finding starred zero in the column and primed zero in the row until, 
#we get to a primed zero which does not have a starred zero in its column.
#At this point reduce the non-zero values of mask matrix by 1. i.e. change prime zeros to starred zeroes.
#Clear the cover matrices and clear any primes i.e. values=2 from mask matrix.
sub stepfive
{
    print "Step 5 \n";

    my $cnt = 0;
    my $done = 0;
    my ($row, $col) = (-1,-1);

    $path[$cnt][0] = $Z0_row;
    $path[$cnt][1] = $Z0_col;
    
    while($done == 0)
    {
	$row = &find_star_in_col($path[$cnt][1]);
	if($row > -1)
	{
	    $cnt++;
	    $path[$cnt][0] = $row;
	    $path[$cnt][1] = $path[$cnt - 1][1];
	}
	else
	{
	    $done = 1;
	}
	if($done == 0)
	{
	    $col = &find_prime_in_row($path[$cnt][0]);
	    $cnt++;
	    $path[$cnt][0] = $path[$cnt - 1][0];
	    $path[$cnt][1] = $col;
	}
    }
    &convert_path($cnt);
    &clear_covers();
    &erase_primes();

    print "Step 5 end \n";
}

#Tries to find starred zero in the given column and returns the row number
sub find_star_in_col
{
    my $col = shift;
    my $row = -1;
    my $i = 0;

    for($i=0;$i<=$#mat;$i++)
    {
	if($mask[$i][$col] == 1)
	{
	    $row = $i;
	}
    }
    
    return $row;
}

#Tries to find primed zero in the given row and returns the column number
sub find_prime_in_row
{
    my $row = shift;
    my $col = -1;
    my $j = 0;

    for($j=0;$j<=$#mat;$j++)
    {
	if($mask[$row][$j] == 2)
	{
	    $col = $j;
	}
    }
    
    return $col;
}

#Reduces non-zero value in the mask matrix by 1.
#i.e. converts all primes to stars and stars to none.
sub convert_path
{
    my $cnt = shift;
    my ($i, $j) = (0,0);

    for($i=0;$i<=$cnt;$i++)
    {
	if($mask[$path[$i][0]][$path[$i][1]] == 1)
	{
	   $mask[$path[$i][0]][$path[$i][1]] = 0;
	}
	else
	{
	   $mask[$path[$i][0]][$path[$i][1]] = 1;
	}
    }
}

#Clears cover matrices
sub clear_covers
{
    my $i;
    for($i=0;$i<=$#mat;$i++)
    {
	$rowcov[$i] = 0;
	$colcov[$i] = 0;
    }
}

#Changes all primes i.e. values=2 to 0.
sub erase_primes
{
    my ($i, $j);

    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($mask[$i][$j] == 2)
	    {
		$mask[$i][$j] = 0;
	    }
	}
    }

}

#Step 6 - Find the minimun value from the rows and columns which are currently not covered.
#Subtract this minimum value from all the elements of the columns which are not covered.
#Add this minimum value to all the elements of the rows which are covered.
#Proceed to step 4.
sub stepsix
{
    print "Step 6 \n";
    my ($i, $j);
    my $minval = 0;

    $minval = &find_smallest();
    
    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($rowcov[$i] == 1)
	    {
		$mat[$i][$j] += $minval;
	    }
	    if($colcov[$j] == 0)
	    {
		$mat[$i][$j] -= $minval;
	    }
	}
    }

    print "Step 6 end \n";
}

#Finds the minimum value from all the matrix values which are not covered.
sub find_smallest
{
    my $minval = 99999999999;
    my ($i, $j);    

    for($i=0;$i<=$#mat;$i++)
    {
	for($j=0;$j<=$#mat;$j++)
	{
	    if($rowcov[$i] == 0 && $colcov[$j] == 0)
	    {
		if($minval > $mat[$i][$j])
		{
		    $minval = $mat[$i][$j];
		}
	    }
	}
    }
    return $minval;
}

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__

=head1 NAME

Algorithm::Munkres - Perl extension for Munkres' solution to classical Assignment problem

=head1 SYNOPSIS

use Algorithm::Munkres;

    @mat = (
        [ 12, 3, 7, 4, 10],
	[ 5, 10, 6, 2, 4],
	[ 8, 5, 1, 4, 9],
	[ 15, 2, 7, 8, 10],
	[ 7, 2, 8, 1, 12],
	);

assign(\@mat);

=head1 DESCRIPTION

This module implements the Munkres' solution to classical Assignment problem.
If we have a MxN matrix of Workers and the Jobs to be assigned, then this module assigns one job to each worker in such a way that
overall/sum of the cost of jobs is the lowest possible value.
  eg: MxN = [12   3   7   4   10
              5  10   6   2    4
              8   5   1   4    9
             15   2   7   8   10
              7   2   8   1   12]

is the input matrix, then this module will return a matrix of Mx1 dimensions, like:
    solution = [3
                4
                2
                1
                0]

where value '3' says that 4th column in the input matrix is the assigned column/solution for the 1st row.
      value '4' says that 5th column in the input matrix is the assigned column/solution for the 2st row and so on.
      Thus for the given input matrix, the optimal assignment i.e. the assignment which minimizes the sum of all job costs is:
  eg: MxN = [0   0   0   4   0
             0   0   0   0   4
             0   0   1   0   0
             0   2   0   0   0 
             7   0   0   0   0]
     

Note: The 'assign' subroutine expects the input array's reference and not the complete array. eg:assign(\@mat);

=head2 EXPORT

"assign" function by default.

=head1 SEE ALSO

http://216.249.163.93/bob.pilgrim/445/munkres.html

=head1 AUTHOR

Anagha Kulkarni, University of Minnesota, Duluth
kulka020@d.umn.edu

Ted Pedersen, University of Minnesota, Duluth
tpederse@d.umn.edu

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2000-2004, Ted Pedersen and Anagha Kulkarni

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut
