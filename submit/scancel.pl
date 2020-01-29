use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 23512412; ($val <= 23512545); $val+=1)
{
	system "scancel $val";
}