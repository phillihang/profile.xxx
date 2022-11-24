
use strict;
#use warnings;
use File::Copy;

my $newDir="files";
unless(-d $newDir)
{
	mkdir $newDir or die $!;
}

my @allFiles=glob("*");
foreach my $subDir(@allFiles)
{
	if((-d $subDir) && ($subDir ne $newDir))
	{
		opendir(SUB,"./$subDir") or die $!;
		while(my $file=readdir(SUB))
		{
			if($file=~/\.gz$/)
			{
				#`cp ./$subDir/$file ./$newDir`;
				copy("$subDir/$file","$newDir") or die "Copy failed: $!";
			}
		}
		close(SUB);
	}
}

