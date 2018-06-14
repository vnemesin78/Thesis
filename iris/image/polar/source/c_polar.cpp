#include "c_polar.hpp"

c_polar :: c_polar ( 	unsigned int nb_radii_max,
						unsigned int nb_dir_max )
{
	initialize();
	setup( nb_radii_max, nb_dir_max );
}

void c_polar :: setup (	unsigned int nb_radii_max,
						unsigned int nb_dir_max )
{
	free();
	initialize();
	
	polar_image = cvCreateImage ( 	cvSize( nb_dir_max, nb_radii_max ),
									IPL_DEPTH_64F,
									1 );
									
	polar_mask = cvCreateImage ( 	cvSize( nb_dir_max, nb_radii_max ),
									IPL_DEPTH_8U,
									1 );
									
	
	data = (double*) polar_image->imageData;
	mask = (unsigned char*) polar_mask->imageData;
			
	_nb_dir_max = nb_dir_max;
	_nb_radii_max = nb_radii_max;
	
}

template <class type> int c_polar :: compute (	const double & x_center,
												const double & y_center,
												unsigned int nb_radii,
												unsigned int nb_dir,
												unsigned int r_min,
												unsigned int r_max,
												const type * img_data,
												const unsigned char * img_mask,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												unsigned int width_step_mask  )
{
	if ( nb_dir > _nb_dir_max || nb_radii > _nb_radii_max )
		return 1;
	
	_nb_dir = nb_dir;
	_nb_radii = nb_radii;
	
	image_2_polar (	(double *) polar_image->imageData ,
					(unsigned char *) polar_mask->imageData ,
					x_center,
					y_center, 
					nb_radii,
					nb_dir,
					(unsigned int) polar_image->widthStep / sizeof( double ), //width_step_p_img
					(unsigned int) polar_mask->widthStep / sizeof( char ), // widthStep(mask)
					r_min,
					r_max,
					img_data, 
					img_mask,
					width,
					height,
					width_step,
					width_step_mask );
	cvSetImageROI( 	polar_image, 
					cvRect ( 0, 0, nb_dir, nb_radii ) );
	return 0;
}

int c_polar :: resize (	double * p_data,
						unsigned char * p_mask,
						unsigned int nb_radii,
						unsigned int nb_dir,
						unsigned int width_step,
						unsigned int mask_width_step )
{


	image_resize ( 	p_data,
					p_mask,
					nb_dir,
					nb_radii,
					width_step, // WS image
					mask_width_step, // WS mask
					data, 
					mask,
					_nb_dir,
					_nb_radii,
					polar_image->widthStep / sizeof( double ), //WS polar image
					polar_mask->widthStep / sizeof( char ) ); //WS polar mask
	return 0;
}

int c_polar :: resize( 	IplImage * p_image,
						IplImage * p_mask )
{
	if ( p_image == 0 || p_mask == 0 )
	{
		*err_stream << "Error : Args in  int c_polar :: resize( IplImage * p_image, IplImage * p_mask );" << endl;
		return -1;
	}
	if ( p_image->depth != IPL_DEPTH_64F ||  p_mask->depth != IPL_DEPTH_8U )
	{
		*err_stream << "Error : Image(s) format in int c_polar :: resize( IplImage * p_image, IplImage * p_mask );" << endl;
		return -2;
	}
	
	
	unsigned int w_1, 
				 w_2, 
				 h_1, 
				 h_2;
	double * image_data;
	unsigned char * mask_data;
	unsigned int ws_1 = p_image->widthStep / sizeof( double ), //width_step
				 ws_2 = p_mask->widthStep / sizeof( char ) ;
	if ( p_image->roi )
	{
		w_1 = p_image->roi->width;
		h_1 = p_image->roi->height;
		image_data = ( ( double * ) p_image->imageData ) + p_image->roi->yOffset * ws_1 + p_image->roi->xOffset;
	}
	else
	{
		w_1 = p_image->width;
		h_1 = p_image->height;
		image_data = ( double * ) p_image->imageData;
	}
	
	if ( p_mask->roi )
	{
		w_2 = p_mask->roi->width;
		h_2 = p_mask->roi->height;
		mask_data = ( ( unsigned char * ) p_mask->imageData ) + p_mask->roi->yOffset * ws_2 + p_mask->roi->xOffset;
	}
	else
	{
		w_2 = p_mask->width;
		h_2 = p_mask->height;
		mask_data = ( unsigned char * ) p_mask->imageData;
	}
	
	if ( w_1 != w_2 || h_1 != h_2 )
	{
		*err_stream << "Error : Different image and mask dimensions in int c_polar :: resize( IplImage * p_image, IplImage * p_mask );" << endl;
		return -3;
	}
	
	return resize( 	image_data,
					mask_data,
					h_1,
					w_1,
					ws_1,
					ws_2 );
}

c_polar :: ~c_polar()
{
	free();
	initialize();
}


void c_polar :: initialize()
{
	
	_nb_dir_max = 0;
	_nb_dir = 0;
	_nb_radii_max = 0;
	_nb_radii = 0;
	polar_image = 0;
	polar_mask = 0;
	data = 0;
	mask = 0;
	err_stream = &cout;
}

void c_polar :: free()
{
	if ( polar_image )
	{
		cvReleaseImage ( &polar_image );
	}
	
	if ( polar_mask )
		cvReleaseImage ( &polar_mask );
}

//~ 
int c_polar :: compute ( 	const double & x_center,
							const double & y_center,
							unsigned int nb_radii,
							unsigned int nb_dir,
							unsigned int r_min,
							unsigned int r_max,
							const IplImage * image,
							const IplImage * mask )
{
	
	if ( image == 0 || mask == 0 )
	{
		*err_stream << "Error : Args in  int c_polar :: compute ( 	const double & x_center, const double & y_center, unsigned int nb_radii, unsigned int nb_dir, unsigned int r_min, unsigned int r_max, const IplImage * image, const IplImage * mask );" << endl;
		return -1;
	}
	if ( mask->depth != IPL_DEPTH_8U )
	{
		*err_stream << "Error : Mask ( Bad format) in  int c_polar :: compute ( 	const double & x_center, const double & y_center, unsigned int nb_radii, unsigned int nb_dir, unsigned int r_min, unsigned int r_max, const IplImage * image, const IplImage * mask );" << endl;
		return -2;
	}

	unsigned int width;
	unsigned int height;
	unsigned int width_step;
	//Mask
	const unsigned char * mask_data;
	unsigned int width_step_mask = mask->widthStep;
	if ( mask->roi )
		mask_data = (unsigned char*) mask->imageData + mask->roi->yOffset * width_step_mask + mask->roi->xOffset;
	else
		mask_data = (unsigned char*) mask->imageData;

	switch( image->depth )
	{
		case ( IPL_DEPTH_8U ):
		{
			const unsigned char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const unsigned char *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return c_polar :: compute (	x_center,
										y_center,
										nb_radii,
										nb_dir,
										r_min,
										r_max,
										image_data,
										mask_data,
										width,
										height,
										width_step,
										width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_8S ):
		{
			const char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const char *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return c_polar :: compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_16U ):
		{
			const unsigned short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const unsigned short int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_16S ):
		{
			const short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const short int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_32S ):
		{
			const long int * image_data;
			width_step = image->widthStep / sizeof( long int );
			if ( image->roi != NULL )
			{
				image_data = ( (const long int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const long int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_32F ):
		{
			const float * image_data;
			width_step = image->widthStep / sizeof( float );
			if ( image->roi != NULL )
			{
				image_data = ( (const float *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const float *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		case ( IPL_DEPTH_64F ):
		{
			const double * image_data;
			width_step = image->widthStep / sizeof( double );
			if ( image->roi != NULL )
			{
				image_data = ( (const double *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const double *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute (	x_center,
								y_center,
								nb_radii,
								nb_dir,
								r_min,
								r_max,
								image_data,
								mask_data,
								width,
								height,
								width_step,
								width_step_mask  );
		}
		break;
		default:
			*err_stream << "Error :Image ( Bad format) in  int c_polar :: compute ( 	const double & x_center, const double & y_center, unsigned int nb_radii, unsigned int nb_dir, unsigned int r_min, unsigned int r_max, const IplImage * image, const IplImage * mask );" << endl;
			return -4;
		break;
	}
}



C_POLAR_DEF(unsigned char)
C_POLAR_DEF(char)

C_POLAR_DEF(unsigned short int)
C_POLAR_DEF(short int)

C_POLAR_DEF(unsigned long int)
C_POLAR_DEF(long int)

C_POLAR_DEF(float)
C_POLAR_DEF(double)






