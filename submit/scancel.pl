use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 66290714; ($val <= 66291114); $val+=1)
{
	system "scancel $val";
}