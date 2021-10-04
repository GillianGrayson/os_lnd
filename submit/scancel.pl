use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 34008681; ($val <= 34008704); $val+=1)
{
	system "scancel $val";
}