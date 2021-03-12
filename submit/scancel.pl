use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 16378826; ($val <= 16379127); $val+=1)
{
	system "scancel $val";
}