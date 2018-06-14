#include "c_mbgc_image_processing.hpp"
c_MBGC_image_processing :: c_MBGC_image_processing ( 	unsigned int width,
															unsigned int height, 
															ostream * _err_stream )
{
	initialize();
	setup ( width, height, _err_stream );
}

int c_MBGC_image_processing :: setup ( 	unsigned int width,
											unsigned int height, 
											ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	if ( !width || !height )
	{
		if ( err_stream )
		{
			*err_stream << "Error in int c_MBGC_image_processing :: setup) " << endl;
			*err_stream << "Invalid image dimensions!" << endl;
		}
		
		return 1;
	}	
	_width = width;
	_height = height;
	label_obj = new c_label( width, height );
	_image = cvCreateImage( cvSize(width, height), IPL_DEPTH_8U, 1 );
	return 0;
	
}

int c_MBGC_image_processing :: process( const IplImage * src_image )
{
	cvResetImageROI( 	_image );
	
	if ( 	(unsigned int) src_image->width < _width	||
			(unsigned int) src_image->height < _height	)
	{
		if ( err_stream )
		{
			*err_stream << "Error in int process( const IplImage * src_image ) " << endl;
			*err_stream << "Invalid image dimensions!" << endl;
		}
		return 1;
	}
	
	//Recherche des pixels avec une valeur de 127
	memset ( 	_image->imageData,
				0,
				sizeof(char) * _image->widthStep * _height );
				
	for ( unsigned int i = 0; i < _height; ++ i )
	{
		for ( unsigned int j = 0; j < _width; ++ j )
		{
			if ( ( (unsigned char * )  src_image->imageData )[ i * src_image->widthStep + j ] == 127 )
				( (unsigned char * ) _image->imageData )[ i * _image->widthStep + j ] = 255;
		}
	}

	//Etiquettage des régions à 127
	label_obj->label( 	(unsigned char *) _image->imageData, 
						_image->widthStep, 
						0,
						_width,
						_height);
	
	//Destruction des régions trop petites
	label_obj->erase_small_component( _height + _width );
	

	//Recherche de la zone d'intérêt
	unsigned int x_roi = 0,
				 y_roi = 0,
				 w_roi = _width,
				 h_roi = _height;
	
	
	unsigned int r_id = label_obj->get_biggest_region_id();
	if ( r_id != 0 )
	{
		//Position
		double x = 0, 
			   y = 0;
		for ( unsigned int i = 0; i < _height; ++ i )
		{
			for ( unsigned int j = 0; j < _width; ++ j )
			{
				if ( label_obj->label_map()[ i * label_obj->width_step() + j ] == r_id )
				{
					x += j;
					y += i;
				}
			}
		}

		x /= label_obj->surfaces()[r_id];
		y /= label_obj->surfaces()[r_id];
		//
		// +----
		// |
		// |
		if ( x < _width / 2 && y < _height / 2 )
		{
			unsigned int x_s = 1,
						 y_s = 1;
						 
			for ( ; x_s < _width; ++ x_s ) 
			{
				if ( label_obj->label_map()[ ( _height / 2 ) * label_obj->width_step() + x_s ] != r_id )
					break;
			}
			for ( ; y_s < _height; ++ y_s ) 
			{
				if ( label_obj->label_map()[ y_s * label_obj->width_step() + _width / 2 ] != r_id )
					break;
			}
			
			x_roi = x_s;
			y_roi = y_s + 6;
			w_roi = _width - x_s;
			h_roi = _height - y_roi;
			
			
		}
		//
		// ----+
		//     |
		//     |
		else if ( x >= _width / 2 && y < _height / 2 )
		{
			unsigned int x_s = _width - 2,
						 y_s = 1;
						 
			for ( ; x_s != (unsigned int ) -1; -- x_s ) 
			{
				if ( label_obj->label_map()[ ( _height / 2 ) * label_obj->width_step() + x_s ] != r_id )
					break;
			}
			for ( ; y_s < _height; ++ y_s ) 
			{
				if ( label_obj->label_map()[ y_s * label_obj->width_step() + _width / 2 ] != r_id )
					break;
			}
			x_roi = 1;
			y_roi = y_s + 6;
			w_roi = x_s - 2;
			h_roi = _height - y_roi;
		}
		// |
		// |
		// +----
		//
		else if ( x < _width / 2 && y >= _height / 2 )
		{
			unsigned int x_s = 1,
						 y_s = _height - 2;
						 
			for ( ; x_s < _width; ++ x_s ) 
			{
				if ( label_obj->label_map()[ ( _height / 2 ) * label_obj->width_step() + x_s ] != r_id )
					break;
			}
			for ( ; y_s != (unsigned int ) -1; -- y_s ) 
			{
				if ( label_obj->label_map()[ y_s * label_obj->width_step() + _width / 2 ] != r_id )
					break;
			}
			
			x_roi = x_s;
			y_roi = 6;
			w_roi = _width - x_s;
			h_roi = y_s - 7;
			

		}
		//      |
		//      |
		// -----+
		//
		else
		{
			unsigned int x_s = _width - 2,
						 y_s = _height - 2;
			for ( ; x_s != (unsigned int ) -1; -- x_s ) 
			{
				if ( label_obj->label_map()[ ( _height / 2 ) * label_obj->width_step() + x_s ] != r_id )
					break;
			}
			for ( ; y_s != (unsigned int ) -1; -- y_s ) 
			{
				if ( label_obj->label_map()[ y_s * label_obj->width_step() + _width / 2 ] != r_id )
					break;
			}
			x_roi = 1;
			y_roi = 6;
			w_roi = x_s - 2;
			h_roi = y_s - 7;
		}
	}
	
	cvCopyImage( 	src_image,
					_image );
	cvSetImageROI( 	_image,
					cvRect ( x_roi,
							 y_roi,
							 w_roi,
							 h_roi ) );
	
	//Suppression de l'entrelacement
	for ( unsigned int i = 1; i < h_roi; i += 2 )
	{
		if ( i < h_roi - 1 )
		{
			for ( unsigned int j = 0; j < w_roi; ++ j )
			{
				((unsigned char*) _image->imageData )[ ( i + y_roi ) * _image->widthStep + ( j + x_roi ) ] = ( ((unsigned char*) _image->imageData )[ ( i + y_roi + 1 ) * _image->widthStep + ( j + x_roi ) ] + ((unsigned char*) _image->imageData )[ ( i + y_roi - 1 ) * _image->widthStep + ( j + x_roi ) ] + 1) / 2;
			}
		}
		else
		{
			for ( unsigned int j = 0; j < w_roi; ++ j )
			{
				((unsigned char*) _image->imageData )[ ( i + y_roi ) * _image->widthStep + ( j + x_roi ) ] = ((unsigned char*) _image->imageData )[ ( i + y_roi - 1) * _image->widthStep + ( j + x_roi ) ];
			}
		}
	}
	
	cvSmooth( _image, _image, CV_MEDIAN, 3, 3 );
	
	return 0;
}

c_MBGC_image_processing :: ~c_MBGC_image_processing()
{
	free();
	initialize();
}

void c_MBGC_image_processing :: free()
{
	if ( _image )
		cvReleaseImage( &_image );
	if ( label_obj )
		delete label_obj;
	
}

void c_MBGC_image_processing :: initialize()
{
	label_obj = 0;
	_image = 0;
	_width = 0;
	_height = 0;
	err_stream = NULL;
}
