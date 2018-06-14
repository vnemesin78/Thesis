#include "c_iris_code.hpp"
/**@fn 
 * @param width : largeur de l'image
 * @param height : hauteur de l'image
 * @param nb_directions : nombre
 * @param nb_samples : nombre d'échantillons sur chq rayons
 * @param wavelenght 
 * @param sigma
 * @brief
 * Constructeur
**/
c_iris_code :: c_iris_code (	unsigned int nb_directions,
								unsigned int nb_samples,
								double wavelenght,
								double sigma,
								ostream * _err_stream )
{
	initialize();
	if ( setup ( 	nb_directions,
					nb_samples,
					wavelenght,
					sigma,
					_err_stream ) )
		throw ( invalid_argument("Arguments of c_iris_code  !") );
}

/**@fn 
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param nb_directions : nombre
			 * @param nb_samples : nombre d'échantillons sur chq rayons
			 * @param wavelenght 
			 * @param sigma
			 * @brief
			 * Setup
			 **/
int c_iris_code :: setup (	unsigned int nb_directions,
							unsigned int nb_samples,
							double wavelenght,
							double sigma,
							ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	//CHeck des arg.
	if ( !nb_directions || !nb_samples || wavelenght <= 0 || sigma <= 0 )
	{
		if ( err_stream )
			*err_stream << "Error in int c_iris_code :: setup( ... ) : Invalid argument(s)!" << endl;
		return 1;
	}
	
	_nb_directions = nb_directions;
	_nb_samples = nb_samples;
	_wavelenght = wavelenght;
	_sigma = sigma;
	
	//Images
	_code_img = cvCreateImage ( cvSize( _nb_directions, 
										_nb_samples * 2),
								IPL_DEPTH_8U,
								1 );
	
	_mask_img = cvCreateImage ( cvSize( _nb_directions,
										_nb_samples ),
								IPL_DEPTH_8U,
								1 );
	
	_g_img_im = cvCreateImage ( cvSize( _nb_directions, 
										_nb_samples ),
								IPL_DEPTH_64F,
								1 );

	_g_img_re = cvCreateImage ( cvSize( _nb_directions, 
										_nb_samples ),
								IPL_DEPTH_64F,
								1 );	
								
	//Tmp
	_array = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * _nb_directions );
	_f_array = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * _nb_directions );
	tmp = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * _nb_directions );
	
	//Filtre log-gabor
	gabor_filter = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * _nb_directions );
	create_gabor_filter ( 	gabor_filter, 
							_nb_directions, 
							_wavelenght, 
							_sigma );
							
	//Plans
	convolve_build_fftw_plans (	plan_in,
								plan_out,
								_f_array,
								_array,
								_nb_directions,
								tmp );
								
	structuring_element_1 = 
		cvCreateStructuringElementEx( 	(_nb_samples / 4) * 2 + 1,
										(_nb_samples / 4) * 2 + 1,
										(_nb_samples / 4),
										(_nb_samples / 4),
										CV_SHAPE_ELLIPSE,
										NULL);		
	tmp1 = 
		cvCreateImage ( cvSize( _nb_directions, _nb_samples), 
						IPL_DEPTH_64F, 
						1 );		
	tmp2 = 
		cvCreateImage ( cvSize( _nb_directions, _nb_samples), 
						IPL_DEPTH_64F, 
						1 );				
	float tmp[4];
	tmp[0] = 1;
	tmp[1] = 1;
	tmp[2] = 1;
	tmp[3] = 1;
	interpol_2d_obj.setup( 	_nb_directions,
							_nb_samples,
							tmp );
	
	return 0;
}

int c_iris_code :: setup(	api_parameters & params,
							ostream * _err_stream,
							const char * prefix,
							const char * nb_directions_name,
							const char * nb_samples_name,
							const char * wavelenght_name,
							const char * sigma_name )
{
	int q = 0;
	unsigned int nb_directions; //Largeur
	unsigned int nb_samples; //Hauteur
	double wavelenght;
	double sigma;
	err_stream = _err_stream;
	stringstream oss;
	oss << prefix << "::" << nb_directions_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&nb_directions, 
									err_stream) )
		q = 1;
	oss.str("");
	
	oss << prefix << "::" << nb_samples_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&nb_samples,
									err_stream ) )
		q = 1;
	oss.str("");
		
	oss << prefix << "::" << wavelenght_name;
	if (api_get_double(	params, 
						oss.str().c_str(),  
						&wavelenght,
						err_stream ) )
		q = 1;
	oss.str("");
		
	oss << prefix << "::" << sigma_name;
	if (api_get_double(	params, 
						oss.str().c_str(),  
						&sigma,
						err_stream) )
		q = 1;
	oss.str("");
		
	if ( q )
		return 1;
		
	return ( setup (	nb_directions, 
						nb_samples, 
						wavelenght, 
						sigma,
						err_stream) );
}


int c_iris_code :: compute_iris_code (	const IplImage * image,
											const IplImage * mask )
{
	unsigned int 	img_width,
					img_height,
					img_x_offset,
					img_y_offset,
					img_width_step,
					mask_width,
					mask_height,
					mask_x_offset,
					mask_y_offset,
					mask_width_step;	
	
	
	if ( !image || !mask )
	{
		if ( err_stream)
			*err_stream << "Error 404: Image or mask not found!" << endl;
		return -1;
	}

	
	GET_IMAGE_DIM( image, 
				   img_width, 
				   img_height,
				   img_width_step,
				   img_x_offset,
				   img_y_offset );
				   
	GET_IMAGE_DIM( mask, 
				   mask_width, 
				   mask_height,
				   mask_width_step,
				   mask_x_offset,
				   mask_y_offset );				   
				   
				   
				   
				   
	
	
	if ( 	img_width != _nb_directions 	|| 
			mask_width != _nb_directions	||
			img_height != _nb_samples		||
			mask_height != _nb_samples		||
			mask->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
		{
			*err_stream << "Error : Image(s) dimension(s)!" << endl;
			*err_stream << "width" << endl;
			*err_stream << img_width << "\t" << mask_width << "\t" << _nb_directions << endl;
			*err_stream << "height" << endl;
			*err_stream << img_height << "\t" << mask_height << "\t" << _nb_samples << endl;		
		}
		return 1;
	}
	
	switch ( image->depth )
	{
		case (IPL_DEPTH_8U):
			compute_iris_code ( (unsigned char* ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(char),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(char),
								mask_width_step );
		break;
		case (IPL_DEPTH_8S):
			compute_iris_code ( (char* ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(char),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(char),
								mask_width_step );
		break;
		case (IPL_DEPTH_16U):
			compute_iris_code ( (unsigned short int * ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(unsigned short int),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(unsigned short int),
								mask_width_step );
		break;
		case (IPL_DEPTH_16S):
			compute_iris_code ( (short int * ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(short int),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(short int),
								mask_width_step );
		break;
		case (IPL_DEPTH_32S):
			compute_iris_code ( (long int * ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(long int),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(long int),
								mask_width_step );
		break;
		case (IPL_DEPTH_32F):
			compute_iris_code ( (float * ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(float),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(float),
								mask_width_step );
		break;
		case (IPL_DEPTH_64F):
			compute_iris_code ( (double * ) image->imageData + img_x_offset + img_width_step * img_y_offset / sizeof(double),
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width_step / sizeof(double),
								mask_width_step );
		break;
		default:
			if ( err_stream )
				*err_stream << "Error: Unsupported image format!" << endl;
			return 1;
		break;
	}
	
	

	return 0;
}


template <class type> void c_iris_code :: compute_iris_code ( 	const type * img_data,
																	const unsigned char * mask,
																	unsigned int img_width_step,
																	unsigned int mask_width_step )
{
	double occ_f = 0;
	
	memset( _g_img_im->imageData, 0, _g_img_im->widthStep * _g_img_im->height );
	memset( _g_img_re->imageData, 0, _g_img_re->widthStep * _g_img_re->height );	
	
	//Recopie du masque
	for ( unsigned int i = 0; i < _nb_samples; ++ i )
	{
		memcpy( ( (unsigned char*)  _mask_img->imageData ) + i * _mask_img->widthStep,
				mask + i * mask_width_step,
				_nb_directions );
	}
	cvMorphologyEx( _mask_img,
					_mask_img,
					NULL,
					structuring_element_1,
					CV_MOP_OPEN,
					1);	
	
	if ( interpol_2d_obj.interpolate(	 img_data, 
										 _mask_img, 
										 img_width_step, 
										 -1 ) )
		return;

	double m = get_image_mean ( interpol_2d_obj.image(), _mask_img, 255 );
	double v = get_image_variance ( interpol_2d_obj.image(), _mask_img, 255, m );
	
	for ( unsigned int i = 0; i < _nb_samples; ++ i )
	{
		//Recopie
		for ( unsigned int j = 0; j < _nb_directions; ++ j )
		{
			_array[j][1] = 0;
			_array[j][0] =( ( (double*) ( interpol_2d_obj.image()->imageData + i * interpol_2d_obj.image()->widthStep ) )[ j ] - m ) / sqrt(v);	
		}	
	}
	
	cvSmooth( 	interpol_2d_obj.image(), 
				tmp1, 
				CV_BLUR, 
				5, 
				1);
	cvSub(	tmp1, 
			interpol_2d_obj.image(), 
			tmp1);
			
	double nrj = 0;
	unsigned int n_pixel = 0;
	//Pour chaque rayon
	for ( unsigned int i = 0; i < _nb_samples; ++ i )
	{
		//Recopie
		for ( unsigned int j = 0; j < _nb_directions; ++ j )
		{
			if ( mask[ i * mask_width_step + j ])
			{
				double v = GET_IMAGE_PIXEL(tmp1, i, j, double);	
				if ( v < 0 )
					v = -v;
				if ( v > 2 )
					v = 2;
				nrj += v;
				n_pixel ++;
			}

		}	
				
	}
	
	_nrj_ratio = nrj / n_pixel * get_image_mean( _mask_img ) / 255.0 / 2;

	
	//Pour chaque rayon
	for ( unsigned int i = 0; i < _nb_samples; ++ i )
	{
		//Recopie
		for ( unsigned int j = 0; j < _nb_directions; ++ j )
		{
			_array[j][1] = 0;
			_array[j][0] =( ( (double*) ( interpol_2d_obj.image()->imageData + i * interpol_2d_obj.image()->widthStep ) )[ j ] - m ) / sqrt(v);	
		}	

		//Filtrage par log-gabor
		convolve( gabor_filter,
				  _nb_directions,
				  tmp,
				  plan_in,
				  plan_out );
		
		//Calcul de l'iris code
		for ( unsigned int j = 0; j < _nb_directions; ++ j )
		{
			double mod = sqrt(_f_array[j][0] * _f_array[j][0] + _f_array[j][1] * _f_array[j][1]);
			( (double*) (_g_img_re->imageData + _g_img_re->widthStep * i ) )[j] = _f_array[j][0] / mod;
			( (double*) (_g_img_im->imageData + _g_img_im->widthStep * i ) )[j] = _f_array[j][1] / mod;	
			
			
					
			int h1 = _f_array[j][0] > 0, //Binarisation de la partie réelle
				h2 = _f_array[j][1] > 0, //Binarisation de la phase
				h3 = ( ( _f_array[j][0] * _f_array[j][0] + _f_array[j][1] * _f_array[j][1] ) > 1.0e-8); //Validité de la phase.
			//Masque
			( (unsigned char*)  _mask_img->imageData )[ i * _mask_img->widthStep + j ] = 
				255 * ( h3 && ( ( (unsigned char*)  _mask_img->imageData )[ i * _mask_img->widthStep + j ] ) ); //A vérifier
			
			if ( ( _mask_img->imageData )[ i * _mask_img->widthStep + j ] == 0 )
			{
				( (unsigned char*)  _code_img->imageData )[ i * _code_img->widthStep + j ] = 
					128; //A vérifier			
				
				( (unsigned char*)  _code_img->imageData )[ ( i + _nb_samples ) * _code_img->widthStep + j ] = 
					128; //A vérifier
			}
			else
			{
				( (unsigned char*)  _code_img->imageData )[ i * _code_img->widthStep + j ] = 
					255 * h1; //A vérifier			
				
				( (unsigned char*)  _code_img->imageData )[ ( i + _nb_samples ) * _code_img->widthStep + j ] = 
					255 * h2; //A vérifier
			}	
				
				
				
				
			if ( h3 )
				occ_f += 1;
		}
	}

}

c_iris_code :: ~c_iris_code()
{
	free();
	initialize();
}

void c_iris_code :: free()
{
	if ( _code_img )
		cvReleaseImage ( &_code_img );
	if ( _mask_img )
		cvReleaseImage( & _mask_img );
	if ( tmp )
	{
		fftw_destroy_plan(plan_out);
		fftw_destroy_plan(plan_in);
	}
	
	if ( _g_img_im )
		cvReleaseImage ( &_g_img_im );
	if ( _g_img_re )
		cvReleaseImage ( &_g_img_re );
	
	if ( _array )
		fftw_free ( _array );
	if ( _f_array )
		fftw_free( _f_array );
	if ( tmp )
	{
		fftw_free( tmp );
	}
	if ( gabor_filter )
		fftw_free ( gabor_filter );
	if ( structuring_element_1 )
			cvReleaseStructuringElement( & structuring_element_1 );
	if (tmp1)
		cvReleaseImage( &tmp1);
	if (tmp2)
		cvReleaseImage( &tmp1);
}
			
void c_iris_code :: initialize()
{
	_code_img = 0;
	_mask_img = 0;
	_nb_directions = 0;
	_nb_samples = 0;
	_wavelenght = 0;
	_nrj_ratio = 0;
	_sigma = 0;
	_array = 0;
	_g_img_im = 0;
	_g_img_re = 0;
	_f_array = 0;
	tmp = 0;
	tmp2 = 0;
	gabor_filter = 0;
	tmp1 = 0;
	err_stream = NULL;
	structuring_element_1 = 0;
}

C_IRIS_CODE_COMPUTE_IRIS_CODE(char)
C_IRIS_CODE_COMPUTE_IRIS_CODE(unsigned char)
C_IRIS_CODE_COMPUTE_IRIS_CODE(short int)
C_IRIS_CODE_COMPUTE_IRIS_CODE(unsigned short int)


C_IRIS_CODE_COMPUTE_IRIS_CODE(long int)
C_IRIS_CODE_COMPUTE_IRIS_CODE(float)
C_IRIS_CODE_COMPUTE_IRIS_CODE(double)

