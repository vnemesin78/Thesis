#include "c_focus_score.hpp"
c_focus_score :: c_focus_score(	unsigned int width,
								unsigned int height,
								const gsl_matrix * kernel,
								double c,
								ostream * _err_stream )
{
	initialize();
	setup (	width,
			height,
			kernel, 
			c,
			_err_stream );
}

int c_focus_score :: setup (	unsigned int width,
								unsigned int height,
								const gsl_matrix * _kernel,
								double c,
								ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	if (	_kernel == 0	|| 
			c == 0			||
			width == 0		||
			height == 0		)
	{
		if ( err_stream )
			*err_stream << "Error : invalid argument(s) in int c_focus_score :: setup" << endl;
		return 1;
	}
	_c = c;
	_width = width;
	_height = height;
	//Kernel
	kernel = cvCreateMat(_kernel->size1,
						 _kernel->size2,
						 CV_32FC1);
	
	for ( unsigned int i = 0; i < _kernel->size1; ++ i )
		for ( unsigned int j = 0; j < _kernel->size2; ++ j )
		kernel->data.fl[ i * kernel->step / sizeof(float) + j ] = _kernel->data[ i * _kernel->tda + j ];
		
	//Objet d'interpolation + image
	tmp_img = cvCreateImage (	cvSize(	width,
										height),
								IPL_DEPTH_32F,
								1 );
	tmp_mask = cvCreateImage (	cvSize(	width,
										height),
								IPL_DEPTH_8U,
								1 );
													
	f_img = cvCreateImage (	cvSize(	width,
									height),
							IPL_DEPTH_32F,
							1 );
	r_img = cvCreateImage (	cvSize(	width,
									height),
							IPL_DEPTH_32F,
							1 );							
	r_mask = cvCreateImage (	cvSize(	width,
									height),
							IPL_DEPTH_8U,
							1 );							
							
							
							
	float factor[4];
	factor[0] = 1.0;
	factor[1] = 1.0;
	factor[2] = 1.0;
	factor[3] = 1.0;
	
	interpol_2d.setup (	width,
						height,
						factor );
							
	return 0;
}

int c_focus_score :: setup(	api_parameters & params,
							ostream * _err_stream,
							const char * prefix,
							const char * width_name,
							const char * height_name,
							const char * kernel_name,
							const char * c_name)
{
	int q = 0;
	double c;
	gsl_matrix _kernel;
	stringstream oss;
	err_stream = _err_stream;
	oss << prefix << "::" << c_name;
	if ( api_get_double( 	params, 
							oss.str().c_str(),
							&c,
							err_stream ) )
		q = 1;
	oss.str("");
	oss << prefix << "::" << kernel_name;
	if ( api_get_matrix ( 	params,
							oss.str().c_str(),
							&_kernel,
							err_stream ) )
		q = 1;
	oss.str("");
	unsigned int width;
	oss << prefix << "::" << width_name;
	if ( api_get_positive_integer ( params,
									oss.str().c_str(),
									&width,
									err_stream ) )
		q = 1;
	oss.str("");
	
	unsigned int height;
	oss << prefix << "::" << height_name;
	if ( api_get_positive_integer ( params,
									oss.str().c_str(),
									&height,
									err_stream ) )
		q = 1;
	oss.str("");
	
								
		
		
		
	if ( q ) 
		return 1;
	return ( setup (	width, 
						height, 
						&_kernel, 
						c,
						err_stream ) );
}


template<class type> double c_focus_score :: get_score ( 	const type * image,
															unsigned int width,
															unsigned int height,
															unsigned int width_step)
{
	//Check
	if ( 	! width 	||
			! height	||
			! image		)
	{
		if ( err_stream )
			*err_stream << "Error : Invalid argument(s) in template<class type> double c_focus_score :: get_score!" << endl;
		return nan("");		
	}
	
	//Check des dims.
	memset( tmp_img->imageData,
			0,
			tmp_img->widthStep * tmp_img->height );
	cvSetImageROI( tmp_img,
				   cvRect( 0, 0, width, height ) );
	cvSetImageROI( 	f_img,
					cvRect(	0,
							0,
							width,
							height ) );	  
				   
	double q = 1 / 255.0;	   
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( (float*) ( tmp_img->imageData + tmp_img->widthStep * i ) )[j]
				= image[ i * width_step + j] * q;
		}		
	}

	cvFilter2D( tmp_img,  
				f_img, 
				kernel );
				
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			float & v = ( (float*) ( f_img->imageData + i * f_img->widthStep ))[j];
			if ( v >= 3.0 )
				v = 3;
		}
	}
	
	double m_0 = get_image_mean( tmp_img ),
			m_f = get_image_mean( f_img ),
			v_f = get_image_variance( f_img, m_f ); 
	if (isnan(m_0) || isnan(m_f) || isnan(v_f) || m_0 == 0 )
		return 0;
	return ( v_f + m_f * m_f ) / ( m_0 * m_0 );
	
	
	
}




template<class type> double c_focus_score :: get_score ( 	const type * image,
															const unsigned char * mask,
															unsigned int width,
															unsigned int height,
															unsigned int width_step,
															unsigned int mask_width_step,
															int x,
															int y,
															unsigned int w,
															unsigned int h,
															unsigned int new_w,
															unsigned int new_h )
{
	//Check
	if ( 	! width 	||
			! height	||
			! image		||
			! mask      ||
			! w			||
			! h			||
			! new_h		||
			! new_w		)
	{
		if ( err_stream )
			*err_stream << "Error : Invalid argument(s) in template<class type> double c_focus_score :: get_score!" << endl;
		return nan("");		
	}
	
	//Check des dims.
	
	
	
	
	//Recopie de l'image
	memset( tmp_img->imageData,
			0,
			tmp_img->widthStep * tmp_img->height );
	cvSetImageROI( tmp_img,
				   cvRect( 0, 0, w, h ) );
				   
	memset( tmp_mask->imageData,
			0,
			tmp_mask->widthStep * tmp_mask->height );				   
	cvSetImageROI( tmp_mask,
				   cvRect( 0, 0, w, h ) );				   
	
	int 	y_start = MAX(0, (int) -y),
			y_end	= MIN(h, (int) y - height + 1),
			x_start = MAX(0, (int) -x),
			x_end 	= MIN(w, (int) x - width + 1);
	
	//~ if ( x_end != (int) w || y_end != (int) h)
		//~ return 0;
		
		
	for ( int i = y_start; i < y_end; ++ i )
	{
		for ( int j = x_start; j < x_end; ++ j )
		{
			( (float*) ( tmp_img->imageData + tmp_img->widthStep * i ) )[j]
				= image[ (i + y) * width_step + (j + x)] / 255.0;
			( (unsigned char *) ( tmp_mask->imageData + tmp_mask->widthStep * i ) )[j]
				= mask[ (i + y) * mask_width_step + (j + x)];	
		}		
	}
	//Redimensionnement de l'image
	memset ( r_img->imageData,
			 0,
			 r_img->widthStep * r_img->height );
	memset ( r_mask->imageData,
			 0,
			 r_mask->widthStep * r_mask->height );
			 			 
			 
	
	cvSetImageROI( r_img,
				   cvRect( 0, 0, new_w, new_h ) );
	cvSetImageROI( r_mask,
				   cvRect( 0, 0, new_w, new_h ) );
				   
	image_resize ( 	(float*) r_img->imageData,
					(unsigned char*) r_mask->imageData,
					new_w,
					new_h,
					r_img->widthStep / sizeof(float), 
					r_mask->widthStep / sizeof(unsigned char), 
					(float*) tmp_img->imageData,
					(unsigned char*) tmp_mask->imageData,
					w,
					h,
					tmp_img->widthStep / sizeof(float),
					tmp_mask->widthStep / sizeof(unsigned char) );
	
	
	
	
	
	//Calcul du score
	if ( interpol_2d.interpolate( 	r_img, 
									r_mask, 
									0 ) )
	{
		if ( err_stream )
			*err_stream << "Error in double c_focus_score :: get_score" << endl;
		return nan("");
	}	
	
	//Filtrage
	interpol_2d.get_image( r_img );	
	cvSetImageROI( 	f_img,
					cvRect(	0,
							0,
							new_w,
							new_h ) );	
	
	cvFilter2D( r_img,  
				f_img, 
				kernel );
				
	for ( unsigned int i = 0; i < new_h; ++ i )
	{
		for ( unsigned int j = 0; j < new_w; ++ j )
		{
			float & v = ( (float*) ( f_img->imageData + i * f_img->widthStep ))[j];
			if ( v >= 1.0 )
				v = 1.0;
		}
	}
				
				
									
	double m_0 = get_image_mean( r_img, r_mask, 255 ),
			m_f = get_image_mean( f_img, r_mask, 255 ),
			v_0 = get_image_variance( r_img, r_mask, 255, m_0 ),
			v_f = get_image_variance( f_img, r_mask, 255, m_f );
	double n = get_image_mean( r_mask ) * new_w * new_h /255.0;
	if (isnan(m_0) || isnan(m_f) || isnan(v_0) || isnan(v_f)  || isnan(n) || m_0 == 0 )
		return 0;
	return ( n * n * ( m_f * m_f + v_f ) / ( m_0 * m_0 + v_0) );
}

double c_focus_score :: get_score (	const IplImage * image,
										const IplImage * mask,
										const CvRect & roi,
										unsigned int new_w,
										unsigned int new_h )
{
	unsigned int 	width,
					height,
					x_offset,
					y_offset,
					width_step;
	
	unsigned int 	m_width,
					m_height,
					m_x_offset,
					m_y_offset,
					m_width_step;
	
	
	
	if ( image == NULL || mask == NULL)
	{
		if ( err_stream )
			*err_stream << "Error : image is NULL!" << endl;
		return nan("");
	}
	
	GET_IMAGE_DIM( 	image, 
					width,
					height,
					width_step,
					x_offset,
					y_offset );
	GET_IMAGE_DIM( 	mask, 
					m_width,
					m_height,
					m_width_step,
					m_x_offset,
					m_y_offset );		
	
	if ( 	m_width 	!= width		||
			m_height	!= height		||
			mask->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
			*err_stream << "Error : MASK!" << endl;
		return nan("");
	}
	
	switch ( image->depth )
	{
		case(IPL_DEPTH_8S):
			return ( 
				get_score ( ( (char *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(char),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);
		break;
		case(IPL_DEPTH_8U):
			return ( 
				get_score ( ( (unsigned char *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(char),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);
		break;
		case(IPL_DEPTH_16S):
			return ( 
				get_score ( ( (short int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(short int),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);		
		break;		
		case(IPL_DEPTH_16U):
			return ( 
				get_score ( ( (unsigned short int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(short int),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);		
		break;		
		case(IPL_DEPTH_32S):
			return ( 
				get_score ( ( (long int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(long int),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);		
		break;		
		case(IPL_DEPTH_32F):
			return ( 
				get_score ( ( (float *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(float),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);		
		break;		
		case(IPL_DEPTH_64F):
			return ( 
				get_score ( ( (double *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							( (unsigned char *) ( mask->imageData + m_width_step * m_x_offset ) ) + m_y_offset,
							width,
							height,
							width_step / sizeof(double),
							m_width_step / sizeof(char),
							roi.x,
							roi.y,
							roi.width,
							roi.height,
							new_w,
							new_h )
					);			
		break;		
		default:
			if ( err_stream )
				*err_stream << "Error : image format!" << endl;
			return nan("");
		break;
		
		
	}

}


double c_focus_score :: get_score (	const IplImage * image )
{
	unsigned int 	width,
					height,
					x_offset,
					y_offset,
					width_step;

	if ( image == NULL )
	{
		if ( err_stream )
			*err_stream << "Error : image is NULL!" << endl;
		return nan("");
	}
	
	GET_IMAGE_DIM( 	image, 
					width,
					height,
					width_step,
					x_offset,
					y_offset );		
	switch ( image->depth )
	{
		case(IPL_DEPTH_8S):
			return ( 
				get_score ( ( (char *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(char) )
					);
		break;
		case(IPL_DEPTH_8U):
			return ( 
				get_score ( ( (unsigned char *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(char) )
					);
		break;
		case(IPL_DEPTH_16S):
			return ( 
				get_score ( ( (short int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(short int) )
					);		
		break;		
		case(IPL_DEPTH_16U):
			return ( 
				get_score ( ( (unsigned short int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(short int) )
					);		
		break;		
		case(IPL_DEPTH_32S):
			return ( 
				get_score ( ( (long int *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(long int)  )
					);		
		break;		
		case(IPL_DEPTH_32F):
			return ( 
				get_score ( ( (float *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(float) )
					);		
		break;		
		case(IPL_DEPTH_64F):
			return ( 
				get_score ( ( (double *) ( image->imageData + width_step * x_offset ) ) + y_offset,
							width,
							height,
							width_step / sizeof(double) )
					);			
		break;		
		default:
			if ( err_stream )
				*err_stream << "Error : image format!" << endl;
			return nan("");
		break;
		
		
	}
}

c_focus_score :: ~c_focus_score()
{
	free();
	initialize();
}

void c_focus_score :: initialize()
{
	tmp_img = 0;
	f_img = 0;
	r_img = 0;
	tmp_mask = 0;
	r_mask = 0;
	_c = 0;
	kernel = 0;
	err_stream = NULL;
	_width = 0;
	_height = 0;
}

void c_focus_score :: free()
{
	if ( kernel )
		cvReleaseMat( &kernel );
	if ( tmp_img )
		cvReleaseImage( &tmp_img );
	if ( tmp_mask )
		cvReleaseImage( &tmp_mask );
	if ( r_img )
		cvReleaseImage( &r_img );
	if ( r_mask )
		cvReleaseImage( &r_mask );
	if ( f_img )
		cvReleaseImage( &f_img );
}

C_FOCUS_SCORE_GET_SCORE(unsigned char)
C_FOCUS_SCORE_GET_SCORE(char)
C_FOCUS_SCORE_GET_SCORE(unsigned short int)
C_FOCUS_SCORE_GET_SCORE(short int)
C_FOCUS_SCORE_GET_SCORE(long int)
C_FOCUS_SCORE_GET_SCORE(float)
C_FOCUS_SCORE_GET_SCORE(double)
