#include "c_image_thread.hpp"

c_image_thread :: c_image_thread( void )
{
	initialize();
}

int c_image_thread :: setup ( void )
{
	free();
	initialize();
	return 0;
}


int c_image_thread :: setup( 	char ** argv,
								unsigned int argc,
								ostream * _err_stream )
{
	free();
	initialize();
	
	if ( argc == 0 )
	{
		if (_err_stream )
			*_err_stream << "Error : Missing argument(s)" << endl;
		return 1;
	}

	err_stream = _err_stream;
	
	api_parameters params;
	for ( unsigned int i = 0; i < argc; ++ i )
		params.load(argv[i]);

	//Variables à charger
	int q = 0;

	unsigned int 	d_width,
				    d_height;

	if ( api_get_positive_integer( 	params,
									"width",
									&d_width ) )
		d_width = DEFAULT_WIDTH;

	if ( api_get_positive_integer( 	params,
									"height",
									&d_height) )
		d_height = DEFAULT_HEIGHT;

	if ( api_get_integer( 	params,
							"data_type",
							&type,
							err_stream ) )
		q = 1;
		
	_data = new image_data;
	_data->setup( cvSize( d_width, d_height ) );

	if (q)
		return 1;

	return 0;
}

int c_image_thread :: process(	const IplImage * image, 
									const char * img_name )
{
	_data->img_ok = false;
	if ( image == NULL )
	{
		if ( err_stream )
			*err_stream << "Error : (int c_image_thread :: process) : image == NULL" << endl;
		return 1;
	}
	if ( image->nChannels != 1 && image->nChannels != 3 )
	{
		if ( err_stream )
			*err_stream << "Error : (int c_image_thread :: process) : ( image->nChannels != 1 && image->nChannels != 3 )" << endl;
		return 1;
	}
	if ( image->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
			*err_stream << "Error : (int c_image_thread :: process) : ( image->depth != IPL_DEPTH_8U )" << endl;
		return 1;
	}
	//Conversion de l'image en niveau de gris
	_data->set_dim( GET_IMAGE_WIDTH(image), GET_IMAGE_HEIGHT(image) );
	
	if ( image->nChannels == 3 )
		cvConvertImage(	image,
						_data->image,
						CV_BGR2GRAY );
	else
		cvCopyImage(	image,
						_data->image );
						
	switch(type)
	{
		case(1): //Entralecement à corriger
		{
			unsigned int 	width = GET_IMAGE_WIDTH(_data->image),
							height = GET_IMAGE_HEIGHT(_data->image),
							h_m1 = height - 1;
			//Suppression de l'entrelacement
			for ( unsigned int i = 1; i < height; i += 2 )
			{
				if ( i < h_m1 )
				{
					for ( unsigned int j = 0; j < width; ++ j )
					{
						((unsigned char*) _data->image->imageData )[ i * _data->image->widthStep + j ]
							= (
								((unsigned char*) _data->image->imageData )[ ( i + 1 ) * _data->image->widthStep + ( j ) ]
								+ ((unsigned char*) _data->image->imageData )[ ( i - 1 ) * _data->image->widthStep + ( j ) ]
								+ 1 )
							/ 2;
					}
				}
				else
				{
					memcpy( _data->image->imageData + i * _data->image->widthStep,
							_data->image->imageData + (i - 1) * _data->image->widthStep,
							width * sizeof(char) );
				}
			}
			
			_data->img_ok = true;
		}
		break;
		default:
		{
			_data->img_ok = true;
		}
		break;

	}
	
	cvSmooth( _data->image, _data->image, CV_MEDIAN, 3,3 );
	c_histogram hist;
	hist.compute(	(unsigned char*) _data->image->imageData,
					_data->image->width,
					_data->image->height,
					_data->image->widthStep );
	hist.stretch(	(unsigned char*) _data->image->imageData,
					_data->image->width,
					_data->image->height,
					_data->image->widthStep,
					0.01,
					0.2 );
	
	
	_data->name = img_name;
	++ ( _data->frame_id );

	return 0;
}

c_image_thread :: ~c_image_thread()
{
	free();
	initialize();
}

void c_image_thread :: free()
{
	if ( _data )
		delete _data;
}

void c_image_thread :: initialize()
{
	_data = 0;
	err_stream = NULL;
}

