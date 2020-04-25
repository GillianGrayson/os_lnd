use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 32390111; ($val <= 32390131); $val+=1)
{
	system "scancel $val";
}