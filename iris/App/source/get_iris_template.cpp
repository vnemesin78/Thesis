#include "lib_iris.hpp"

int main( int argc, char ** argv )
{
	c_get_iris_template obj;

	if ( argc == 1 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}

	if ( 	! strcmp( argv[1], "--help" )	||
			! strcmp( argv[1], "-H" ) )
	{
		cout << "HELP" << endl;
		cout << " Arg[1] : Image filename " << endl;
		cout << " Arg[2] : Template filename " << endl;
		cout << " Arg[3 - n] : Parameters files (opt)" << endl;
		cout << "Author: Valérian Némesin." << endl;
		return 0;
	}

	if ( argc < 4 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}

	int q = 0;

	IplImage * image = cvLoadImage( argv[1], 0 );


	//Image location
	if ( image == NULL )
	{
		cout << "Error : Can't load " << argv[1] << "!" << endl;
		q = 1;
	}

	if ( obj.setup( argv + 3, argc - 3 ) )
		q = 1;
	if (q)
	{
		cout << "Error : Invalid parameter(s)" << endl;
		return 1;
	}

	//Template extraction
	obj.segment( image,  argv[1]);

	c_display_iris_data display( image->width,
								 image->height,
								 image->width / 2,
								 image->height + obj.data().code->height * (image->width / ((double) obj.data().code->width)) +  obj.data().polar_image->height * (image->width / ((double) obj.data().code->width)),
								 10,
								 INT_RGB(255,255,255),
								 3,
								 INT_RGB(0,255,0),
								 INT_RGB(255,0,0),
								 INT_RGB(255,0,0),
								 INT_RGB(0,0,255),
								 INT_RGB(255,0,0),
								 obj.data().code->width,
								 obj.data().code->height,
								 obj.data().polar_image->width,
								 obj.data().polar_image->height );


	display.display( obj.data() );

	IplImage * tmp = cvCloneImage( display.polar_image() );
	cvSetImageROI( tmp, cvRect(0,0, display.polar_image()->width, (display.polar_image()->height * obj.data().nb_samples_iris ) / obj.data().nb_samples ) );
//~ 
	cvShowImage( "polar_image", tmp );
	cvSaveImage( "texture.png", tmp );
	cvShowImage( "Iris Code", display.iris_code() );
	cvSaveImage( "code.png", display.iris_code() );
	cvShowImage( "Image iris", display.iris_image() );
	cvSaveImage( "iris.png", display.iris_image() );
	cvShowImage( "Data", display.data_image() );
	
	cvWaitKey(0);	


	//Save
	obj.data().save(argv[2], 42 );

	return 0;
}
