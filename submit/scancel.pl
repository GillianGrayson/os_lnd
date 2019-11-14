use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 18692943; ($val <= 18693042); $val+=1)
{
	system "scancel $val";
}