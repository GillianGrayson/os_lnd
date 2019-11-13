use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 16639843; ($val <= 16639893); $val+=1)
{
	system "scancel $val";
}