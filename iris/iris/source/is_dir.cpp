#include "utils.hpp"
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include <cstdio>
int is_dir( const char * rep_name ) 
{
	struct stat status;
	stat( rep_name, &status );

	if ( status.st_mode & S_IFDIR )
	 return 1;
	else
	 return 0;
}

int file_exist( const char * file_name )
{
	FILE * file = fopen (file_name, "r" );
	if ( !file )
		return 0;
	else
	{
		fclose(file);
		return 1;
	}
}

int create_file( const char * file_name )
{
	if ( file_exist( file_name ) )
		return 1;
	else
	{
		FILE * file = fopen (file_name, "w" );
		if ( !file )
			return -1;
		fclose(file);
		return 0;
	}
	
}

