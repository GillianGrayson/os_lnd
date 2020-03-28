use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 30040423; ($val <= 30040433); $val+=1)
{
	system "scancel $val";
}