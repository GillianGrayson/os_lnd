use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 36187000; ($val <= 36189000); $val+=1)
{
	system "scancel $val";
}