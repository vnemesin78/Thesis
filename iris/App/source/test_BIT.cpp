#include "lib_iris.hpp"
#include <iostream>
using namespace std;
int main (int argc, char ** argv )
{
	unsigned int success = 0, total = 0;
	stringstream oss, oss2, oss3;
	if ( argc < 4 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	
	c_get_iris_template * obj;
	obj = new c_get_iris_template;
	if ( obj->setup ( argv + 3, argc - 3, &cout ) )
		return 1;
		
	mkdir( argv[2], 014777 );
	for ( unsigned int i = 0; i < 1000; ++ i )
	{

		oss2 << argv[2] << "/" << i;
		mkdir( oss2.str().c_str(), 014777 );
		oss2.str("");
		
		oss3.str("");
		if ( i < 10 )
			oss3 << "00" << i;
		else if ( i < 100 )
			oss3 << "0" << i;
		else
			oss3 << i;
		
		
		//Oeil gauche
		oss2 << argv[2] << "/" << i << "/L";
		mkdir( oss2.str().c_str(), 014777 );

		for ( unsigned int j = 0; j  < 10; ++ j )
		{
			oss << argv[1] << "/" << oss3.str() << "/L/S5" << oss3.str() << "L0" << j << ".jpg";
			cout << oss.str().c_str() << endl;
			IplImage * image = cvLoadImage( oss.str().c_str(), 0 );
			if ( image ) 
			{
				total ++;
				if ( ! obj->segment( image, "NULL") )
				{
					success ++;
				}
				obj->data().save( oss2.str().c_str(), j );	
				
			}
			oss.str("");
			cvReleaseImage(&image);
		}
		oss2.str("");
		
		
		
		oss2 << argv[2] << "/" << i << "/R";
		mkdir( oss2.str().c_str(), 014777 );
		for ( unsigned int j = 0; j  < 10; ++ j )
		{
			oss << argv[1] << "/" << oss3.str() << "/R/S5" << oss3.str() << "R0" << j << ".jpg";
			cout << oss.str().c_str() << endl;
			IplImage * image = cvLoadImage( oss.str().c_str(), 0 );
			if ( image ) 
			{
				total ++;
				if ( ! obj->segment( image, "NULL") )
				{
					success ++;
				}
				obj->data().save( oss2.str().c_str(), j );	
				cvReleaseImage(&image);
			}
			oss.str("");
		}
		oss2.str("");
		oss3.str("");
		
		
		
		
	}
	
	cout << success << "/" << total << endl;
	cout << (10000 * success / total ) / 100.0 << "%" << endl;
	
	
	
	
	
	return 0;
}
