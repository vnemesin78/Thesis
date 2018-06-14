#include "c_eyelids_segmentation.hpp"

c_eyelids_segmentation :: c_eyelids_segmentation ( void )
{
	initialize();
}

c_eyelids_segmentation :: c_eyelids_segmentation ( 	unsigned int width,
													unsigned int height,
													unsigned int nb_samples,
													unsigned int nb_directions,
													unsigned int nb_iter_gem,
													unsigned int nb_iter_icm,
													unsigned int nb_iter_em,
													double delta,
													unsigned int closing,
													unsigned int opening,
													ostream * _err_stream
													)
{
	initialize();
	setup( 	 width, height, nb_samples, nb_directions,
			nb_iter_gem, nb_iter_icm, nb_iter_em, delta, closing, opening,
			_err_stream );
}
						
int c_eyelids_segmentation :: setup ( void )
{
	free();
	initialize();
	return 0;
}

int c_eyelids_segmentation :: setup( 	api_parameters & params, 
											unsigned int width,
											unsigned int height,
											ostream * _err_stream,
											const char * n_space, 
											const char * nb_directions_name,
											const char * nb_samples_name,
											const char * delta_name,
											const char * nb_iter_gem_name,
											const char * nb_iter_icm_name,
											const char * nb_iter_em_name,
											const char * closing_name,
											const char * opening_name
										)
{
	free();

	initialize();
	err_stream = _err_stream;
	stringstream oss;
	int q = 0;
	
	unsigned int nb_iter_gem;
	oss << n_space << "::" << nb_iter_gem_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_iter_gem,
									err_stream	) )
		q = 1;
	oss.str("");	
	
	
	unsigned int closing;
	oss << n_space << "::" << closing_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&closing,
									err_stream	) )
		q = 1;
	oss.str("");		
	
	unsigned int opening;
	oss << n_space << "::" << opening_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&opening,
									err_stream	) )
		q = 1;
	oss.str("");		
	
	
	unsigned int nb_iter_icm;
	oss << n_space << "::" << nb_iter_icm_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_iter_icm,
									err_stream	) )
		q = 1;
	oss.str("");	
	
	
	
	unsigned int nb_iter_em;
	oss << n_space << "::" << nb_iter_em_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_iter_em,
									err_stream	) )
		q = 1;
	oss.str("");		
	
	
	double delta;
	oss << n_space << "::" << delta_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&delta,
							err_stream	) )
		q = 1;
	oss.str("");	
	
	unsigned int nb_directions;
	oss << n_space << "::" << nb_directions_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_directions,
									err_stream	) )
		q = 1;
	oss.str("");	
	
	unsigned int nb_samples;
	oss << n_space << "::" << nb_samples_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_samples,
									err_stream	) )
		q = 1;
	oss.str("");	
	
	
	
	
	if ( q )
		return -1;
	
	return (
		setup(	width, height, nb_samples, nb_directions, nb_iter_gem, 
				nb_iter_icm, nb_iter_em, delta, closing, opening, err_stream ) 
			);
	
}

int c_eyelids_segmentation :: setup ( 		unsigned int width,
											unsigned int height,
											unsigned int nb_samples,
											unsigned int nb_directions,
											unsigned int nb_iter_gem,
											unsigned int nb_iter_icm,
											unsigned int nb_iter_em,
											double delta,
											unsigned int closing,
											unsigned int opening,
											ostream * _err_stream)
{
	free();
	initialize();
	err_stream = _err_stream;
	_nb_samples = nb_samples;
	_nb_radii = nb_directions;
	_width = width;
	_height = height;
	_nb_iter_gem = nb_iter_gem;
	_nb_iter_icm = nb_iter_icm;
	_nb_iter_em = nb_iter_em;
	_closing = closing;
	_opening = opening;
	_delta = delta;
	
	
	
	gibbs_params gibbs;
	gibbs.a = 1;
	gibbs.delta = 0; // _delta;
	
	Mouloud_obj_2.setup (	1,
							1,
							log_gaussian_pdf,
							gaussian_estimation,		
							NULL,
							0,
							log_gaussian_pdf,
							gaussian_estimation,
							NULL,
							0,			
							log_gibbs_pdf,
							&gibbs,
							sizeof( gibbs_params ),
							_nb_radii,
							_nb_samples );
	
	gibbs.delta = _delta;
	Mouloud_obj.setup (	_nb_iter_gem,
						_nb_iter_icm,
						log_gaussian_pdf,
						gaussian_estimation,		
						NULL,
						0,
						log_gaussian_pdf,
						gaussian_estimation,
						NULL,
						0,			
						log_gibbs_pdf,
						&gibbs,
						sizeof( gibbs_params ),
						_nb_radii,
						_nb_samples );
						
	float toto[4];
	toto[0] = 0;
	toto[1] = 0;
	toto[2] = 1;
	toto[3] = 1;				
	alloc();
	inter->setup( _nb_radii, _nb_samples, toto );
	return 0;
}

c_eyelids_segmentation ::  ~c_eyelids_segmentation()
{
	free();
	initialize();
}	

int c_eyelids_segmentation :: segment_eyelids ( 	const IplImage * polar_image,
													const IplImage * polar_mask,
													const IplImage * mask,
													unsigned int nb_samples_iris,
													const double & x_pupil,
													const double & y_pupil,
													const double * pupil_radii,
													const double * iris_radii,
													double mean_pupil,
													double sigma_pupil )
{	

	
	if (	!polar_image		||
			!polar_mask			||
			!mask				||
			!nb_samples_iris	||
			!pupil_radii		||
			!iris_radii			)
	{
		if ( err_stream )
			*err_stream << "Error: Invalid argument(s) in int c_eyelids_segmentation :: segment_eyelids" << endl;
		return -1;
	}
	if ( polar_image->depth != IPL_DEPTH_8U )
	{
		IplImage * tmp_img = 
			cvCreateImage (	cvSize( GET_IMAGE_WIDTH(polar_image), 
									GET_IMAGE_HEIGHT(polar_image) ), 
							IPL_DEPTH_8U, 
							1 );
		
		//Conversion
		if ( convert_image (	tmp_img, 
								polar_image ) )
		{
			cvReleaseImage ( &tmp_img );
			return 1;
		}
								
		image_resize (	(unsigned char *) _p_img->imageData,
						(unsigned char *) _p_mask->imageData,
						_nb_radii,
						_nb_samples,
						_p_img->widthStep, 
						_p_mask->widthStep, 
						(unsigned char *) tmp_img->imageData, 
						(unsigned char *) get_image_data(polar_mask), 
						tmp_img->width,
						tmp_img->height,
						tmp_img->widthStep,
						polar_mask->widthStep );
						
		cvReleaseImage ( &tmp_img );
		
	}
	else
	{
		image_resize (	(unsigned char *) _p_img->imageData,
						(unsigned char *) _p_mask->imageData,
						_nb_radii,
						_nb_samples,
						_p_img->widthStep, 
						_p_mask->widthStep, 
						(unsigned char *) get_image_data(polar_image), 
						(unsigned char *) get_image_data(polar_mask), 
						GET_IMAGE_WIDTH( polar_image ),
						GET_IMAGE_HEIGHT( polar_mask ),
						polar_image->widthStep,
						polar_mask->widthStep );		
	}
	
	
	if ( structuring_element_2 )
		cvMorphologyEx( _p_img,
						_p_img,
						NULL,
						structuring_element_2,
						CV_MOP_OPEN,
						1);	
	if ( structuring_element_1 )
		cvMorphologyEx( _p_img,
						_p_img,
						NULL,
						structuring_element_1,
						CV_MOP_CLOSE,
						1);			
	inter->interpolate( _p_img, _p_mask, 10);					
	inter->get_image( _p_img );
	memset( _p_mask->imageData, 255, _p_mask->widthStep * _p_mask->height);
	
	gaussian_params obj_params,
					bg_params;
	for ( unsigned int i = 0; i < _nb_radii; ++ i )
	{
		_radii[i] = ( nb_samples_iris * _nb_samples )/ GET_IMAGE_HEIGHT(polar_image);
	}
	
	
	if ( Mouloud_obj.compute_edges( 	_radii,
										(void *) &bg_params,
										(void *) &obj_params,
										_p_img,
										_p_mask ) )
	{
		//*err_stream << "Warning : Eyelids segmentation failed!" << endl;
		for ( unsigned int i = 0; i < _nb_radii; ++ i )
		{
			_radii[i] = ( nb_samples_iris * _nb_samples )/ GET_IMAGE_HEIGHT(polar_image);
		}
	}
	

	for ( unsigned int i = 0; i < _nb_radii; ++ i )
	{
		if ( isnan( _radii[i]) )
			return 1;
		_radii[i] *= ( GET_IMAGE_HEIGHT(polar_image) / ( ( double) _nb_samples ) );	
		if ( _radii[i] > nb_samples_iris )
			_radii[i] = nb_samples_iris;
	}	

	
	//Construction du masque
	compute_masks( 	pupil_radii, 
					iris_radii, 
					x_pupil, 
					y_pupil, 
					mask, 
					polar_mask, 
					nb_samples_iris );
	
	return 0;
}

void c_eyelids_segmentation :: compute_masks(	const double * pupil_radii,
												const double * iris_radii,
												const double & x_pupil,
												const double & y_pupil,
												const IplImage * mask,
												const IplImage * polar_mask,
												unsigned int nb_samples_iris )
{
	unsigned int nb_directions = GET_IMAGE_WIDTH ( polar_mask );
	if ( _polar_mask )
	{
		if (	(unsigned int)_polar_mask->width < GET_IMAGE_WIDTH ( polar_mask ) ||
				(unsigned int)_polar_mask->height < nb_samples_iris )
		{
			cvReleaseImage ( &_polar_mask );
			polar_mask = 0;
		}
		
	}
	if (! _polar_mask )
	{
		_polar_mask = cvCreateImage ( 	cvSize( GET_IMAGE_WIDTH ( polar_mask ), nb_samples_iris  ),
										IPL_DEPTH_8U,
										1 );
	}
	
	//Masque polaire
	cvSetImageROI( _polar_mask, 
				   cvRect( 0, 0, nb_directions, nb_samples_iris ) );
	
	for ( unsigned int i = 0; i < nb_samples_iris; ++ i )
	{
		memcpy( _polar_mask->imageData + i *_polar_mask->widthStep,
				 polar_mask->imageData + i * polar_mask->widthStep,
				 nb_directions * sizeof(char) );
	}





	double dx = ( (double) _nb_radii ) / nb_directions;
	
	for ( unsigned int i = 0; i < nb_directions; ++ i )
	{

		double 	x = i * dx;
		int y = (int) x;
		double r =  ( 1 - x + y ) * _radii[y % _nb_radii] + _radii[ (y + 1) % _nb_radii] * ( x - y );
		for ( unsigned int j = r; j < nb_samples_iris; ++ j )
		{
			( (unsigned char*) _polar_mask->imageData + j * _polar_mask->widthStep )[i] = 0;
		}

	}	
	for ( unsigned int i = 0; i < _nb_radii; ++ i )
	{
		double angle;		
		angle = ( 2 * i - 1.0) * M_PI / _nb_radii;
		_points[i].x = x_pupil + cos( angle ) * ( pupil_radii[i] + _radii[i] * ( iris_radii[i] - pupil_radii[i]  ) / nb_samples_iris );
		_points[i].y = y_pupil + sin( angle ) * ( pupil_radii[i] + _radii[i] * ( iris_radii[i] - pupil_radii[i]  ) / nb_samples_iris );
	}

	//Masque carthésien
	memset ( _mask->imageData,
			 0,
			 _mask->height * _mask->widthStep );
	
	cvSetImageROI( _mask,
				   cvRect( 0, 
						   0,
						   GET_IMAGE_WIDTH(mask),
						   GET_IMAGE_HEIGHT(mask) ) );
	
	int nb_radii = _nb_radii;   

	cvFillPoly(	_mask,
				&_points,
				&nb_radii,
				1,
				CV_RGB(255, 255, 255) );

	
	for ( unsigned int i = 0; i < (unsigned int) mask->height; ++ i )
	{
		for ( unsigned int j = 0; j < (unsigned int) mask->width; ++ j )
		{
			
			if ( ( (unsigned char*) (_mask->imageData + _mask->widthStep * i ) )[j] == 255 )
			{
				( (unsigned char*) (_mask->imageData + _mask->widthStep * i ) )[j]
					= ( (unsigned char*) (mask->imageData + mask->widthStep * i ) )[j];
			}
		}

	}
					
}

void c_eyelids_segmentation :: free()
{
	if ( _mask )
		cvReleaseImage ( &_mask );
	if ( _polar_mask )
		cvReleaseImage ( &_polar_mask );
	if ( _p_mask )
		cvReleaseImage ( &_p_mask );
	if ( _p_img )
		cvReleaseImage ( &_p_img  );
		
	if ( _radii )
		delete[] _radii;
	if ( _points )
		delete[] _points;
	if ( structuring_element_1 )
		cvReleaseStructuringElement(&structuring_element_1);
	if ( structuring_element_2 )
		cvReleaseStructuringElement(&structuring_element_2);
	if ( inter )
		delete inter;
}

void c_eyelids_segmentation :: initialize()
{
	_width = 0;
	_height = 0;
	_nb_samples = 0;
	_nb_radii = 0;
	structuring_element_1 = 0;
	structuring_element_2 = 0;
	_closing = 0;
	_opening = 0;
	_nb_iter_gem = 0;
	_nb_iter_icm = 0;
	_nb_iter_em = 0;
					  
	_delta = 0;
	_mask = 0;
	_polar_mask = 0;
	_p_img = 0;
	_p_mask = 0;
	_radii = 0;
	_points = 0;
	err_stream = NULL;
	inter = 0;
}

void c_eyelids_segmentation :: alloc()
{
	if ( _width && _height )
		_mask = cvCreateImage ( cvSize( _width, _height),
								IPL_DEPTH_8U,
								1 );
	if ( _nb_radii && _nb_samples )
	{
		_p_img = cvCreateImage ( 	cvSize( _nb_radii, _nb_samples ),
										IPL_DEPTH_8U,
										1 );
		_p_mask = cvCreateImage ( 	cvSize( _nb_radii, _nb_samples ),
										IPL_DEPTH_8U,
										1 );
		
		_radii = new double[_nb_radii];
		_points = new CvPoint[_nb_radii];
	}
	if ( _opening )
	{
		structuring_element_1 = cvCreateStructuringElementEx( 	_opening,
																_opening,
																_opening / 2,
																_opening / 2,
																CV_SHAPE_ELLIPSE,
																NULL);
	}
	
	if ( _closing )
	{
		structuring_element_2 = cvCreateStructuringElementEx( _closing,
															  _closing,
															  _closing / 2,
															  _closing / 2,
															  CV_SHAPE_CROSS,
															  NULL);
	}
	inter = new c_interpol_2d;
	
}




