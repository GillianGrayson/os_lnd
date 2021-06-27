use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 27610248; ($val <= 27610448); $val+=1)
{
	system "scancel $val";
}