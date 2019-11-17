use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 18746045; ($val <= 18746145); $val+=1)
{
	system "scancel $val";
}