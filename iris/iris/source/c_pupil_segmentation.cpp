#include "c_pupil_segmentation.hpp"
void draw_histogram( const char * windows_name, const c_histogram * hist, unsigned int * minima = NULL, unsigned int nb_minima = 0 )
{
	IplImage * tmp = cvCreateImage(  cvSize( 3 * hist->nb_bins(), 480 ), IPL_DEPTH_8U, 3 );
	memset(tmp->imageData, 0, tmp->height * tmp->widthStep * sizeof(char) );
	double max = 0;
	for ( unsigned int i = 0; i < hist->nb_bins(); ++ i )
	{
		if ( hist->s_data()[i] > max )
			max = hist->s_data()[i];
		
	}
	
	for ( unsigned int i = 0; i < hist->nb_bins(); ++ i )
	{
		unsigned char r = 0, g = 0, b = 0;
		r = 255;
		
		cvRectangle( tmp, cvPoint( 3*i, 479), cvPoint( 3 * (i +1), 480 * (1 - hist->s_data()[i] / max) ), CV_RGB( r, g , b ), CV_FILLED );
		
	}
	
	if ( nb_minima > 0 )
	{
		for ( unsigned int i = 0; i <  nb_minima; ++ i )
		{
			unsigned char r = 0, g = 0, b = 0;
			g = 255;	
			cvRectangle( tmp, cvPoint( 3*minima[i], 479), cvPoint( 3 * (minima[i] +1), 0), CV_RGB( r, g , b ), CV_FILLED );
		}
		
	}
	
	cvShowImage( windows_name, tmp);
	
	
	
	
	cvReleaseImage(&tmp);
	
	
}

void c_pupil_segmentation :: get_abs_diff_image( const unsigned char * img_data,
												 unsigned int x,
												 unsigned int y, 
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step )
{
	
	IplImage* image = cvCreateImageHeader( 	cvSize(width, height), 
											IPL_DEPTH_8U,
											1 );
	image->widthStep = width_step;
	image->imageData = (char*) ( img_data + x + y * width_step );
	
	
	//ROIs

	cvSetImageROI( _abs_diff_image_8b, cvRect( 0,0, width, height ) );
	
	//Filtrage de l'image
	cvSmooth( 	image, 
				_abs_diff_image_8b,
				CV_MEDIAN,
				_window_size,
				_window_size  );
	
	//Calcul de la diff
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( (unsigned char*) (_abs_diff_image_8b->imageData + _abs_diff_image_8b->widthStep * i ))[j] = 
				abs ( 
					( (int) ( (unsigned char*) (_abs_diff_image_8b->imageData + _abs_diff_image_8b->widthStep * i ))[j] ) 
					- ( (int) ( img_data + x + y * width_step )[i * width_step + j] ) 
					);
		}
	}

	
	//Debug
	//~ cvShowImage( "Image", image );
	//~ cvShowImage( "Diff", _abs_diff_image_8b );
	//~ cvWaitKey(0);
	
	cvReleaseImageHeader( &image );
}

void c_pupil_segmentation :: compute_histograms( const unsigned char * img_data,
												 unsigned int x,
												 unsigned int y, 
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step )
{
	hist_1.compute(	img_data + x + y * width_step,
					width,
					height,
					width_step,
					_sigma_hist_1 );
	
	hist_2.compute(	(unsigned char*) (_abs_diff_image_8b->imageData),
					width,
					height,
					_abs_diff_image_8b->widthStep,
					_sigma_hist_2 );

	//~ draw_histogram( "Hist image", &hist_1 );
	//~ draw_histogram( "Hist diff", &hist_2 );
	//~ cvWaitKey(0);
}

void c_pupil_segmentation :: mode_detection(	unsigned int width,
												unsigned int height )
{
	hist_1.search_minima(	m_image,
							nb_modes_image,
							256 );
	
	

	for ( unsigned int i = 0; i < nb_modes_image; ++ i )
	{
		m_image[i] += _sigma_hist_1/2;
	}
	double ratio = _width * _height / ( (double) width * height ) * 20;
	if ( ratio >= 100 )
		ratio = 99;
	
	m_d_image[0] = hist_2.get_percentile( ratio );
	//~ m_d_image[1] = hist_2.get_percentile( 20 * ratio );
	//~ m_d_image[2] = hist_2.get_percentile( 30 * ratio );

	nb_modes_d_image = 1;
	
		
	//~ hist_2.search_minima(	m_d_image,
							//~ nb_modes_d_image,
							//~ 256 );

	//~ draw_histogram( "Hist image", &hist_1, m_image, nb_modes_image );
	//~ draw_histogram( "Hist diff", &hist_2, m_d_image, nb_modes_d_image );
	//~ cvWaitKey(0);
}


void c_pupil_segmentation :: mask( 	 const unsigned char * img_data,
									 unsigned int x,
									 unsigned int y, 
									 unsigned int width,
									 unsigned int height,
									 unsigned int width_step,
									 unsigned int s1, 
									 unsigned int s2  )
{
	IplImage* image = cvCreateImageHeader( 	cvSize(width, height), 
											IPL_DEPTH_8U,
											1 );
	image->widthStep = width_step;
	image->imageData = (char*) ( img_data + x + y * width_step );
	
	//ROIs
	cvSetImageROI( mask_1, cvRect( 0,0, width, height ) );
	cvSetImageROI( mask_2, cvRect( 0,0, width, height ) );
	
	cvThreshold( 	image, 
					mask_1, 
					s1, 
					255, 
					CV_THRESH_BINARY_INV);	 
	cvThreshold( 	_abs_diff_image_8b, 
					mask_2, 
					s2, 
					255, 
					CV_THRESH_BINARY_INV);	 

	cvMorphologyEx( mask_1,
					mask_1,
					NULL,
					structuring_element_2,
					CV_MOP_CLOSE,
					1);	

	
	
	cvAnd( mask_2, mask_1, mask_1);


	cvMorphologyEx( mask_1,
					mask_1,
					NULL,
					structuring_element_1,
					CV_MOP_CLOSE,
					1);	
	

	//~ cvShowImage( "Mask image", mask_1 );
	//~ cvWaitKey(0);

	

	cvReleaseImageHeader( &image );
}

int c_pupil_segmentation :: label( )
{
	unsigned int region_id;
	
	//Label des régions noires
	bg_regions->label(	(unsigned char*) mask_1->imageData,
						mask_1->widthStep,
						255,
						GET_IMAGE_WIDTH(mask_1),
						GET_IMAGE_HEIGHT(mask_1));
	//Pas d'objet
	if ( ( region_id = bg_regions->get_biggest_region_id() ) == 0 )
		return 1;
	
	bg_regions->copy_region(	(unsigned char*) mask_1->imageData,
								region_id,
								mask_1->widthStep );
	cvMorphologyEx( mask_1,
					mask_1,
					NULL,
					structuring_element_2,
					CV_MOP_CLOSE,
					1);	


		
	//Objets connexes sombres
	connex_regions->label(	(unsigned char*) mask_1->imageData,
							mask_1->widthStep,
							255,
							GET_IMAGE_WIDTH(mask_1),
							GET_IMAGE_HEIGHT(mask_1));
		
	
	
	return 0;
}

int c_pupil_segmentation :: erase_small_regions()
{
	_nb_valid_regions = 0;
	for (unsigned int j = 1; j < connex_regions->nb_regions(); ++ j)
	{
		unsigned int 	x_region,
						y_region,
						w_region,
						h_region;
		connex_regions->get_bounding_box(	x_region,
											y_region,
											w_region,
											h_region,
											j );
		if (	w_region >= 2 * _radius_min &&
				h_region >= 2 * _radius_min &&
				connex_regions->surfaces()[j] > _radius_min * _radius_min * M_PI	)
		{
			_valid_region_labels[_nb_valid_regions] = j;
			_nb_valid_regions ++;
			
		}					 	
	}

	if ( ! _nb_valid_regions )
		return 1;
		
	return 0;
}


void c_pupil_segmentation :: fit_and_select_best_region()
{

	for (unsigned int j = 0; j < _nb_valid_regions; ++ j)
	{
		
	connex_regions->copy_region( 	(unsigned char*) mask_1->imageData, 
									 _valid_region_labels[j], 
									mask_1->widthStep );


		
		
		//Contour
		contour_obj->compute_8c( connex_regions->label_map(),
								 connex_regions->width(),
								 connex_regions->height(),
								 connex_regions->width_step(),
								 _valid_region_labels[j] );
		//Réchantillonnage euclidien
		digitalize_euclide_contour(	_contour,
									_nb_points_contour,
									contour_obj->contour(),
									contour_obj->nb_points_contour() ); 	
	
		//Ransac
		double 	a_tmp, 
				b_tmp, 
				theta_tmp, 
				x_tmp, 
				y_tmp;
				
				
				
		double d = fit_ellipse( x_tmp, y_tmp, a_tmp, b_tmp, theta_tmp );
		if ( d < _score ) 
		{
			_score = d;
			_x = x_tmp;
			_y = y_tmp;
			_a = a_tmp;
			_b = b_tmp;
			_theta = theta_tmp;
		}

		//~ cvEllipse( mask_1, cvPoint( x_tmp, y_tmp), cvSize(a_tmp, b_tmp), theta_tmp / M_PI * 180, 0, 360, CV_RGB(128,128,128), 2 );
		//~ cvShowImage( "Mask", mask_1);
		//~ cvWaitKey(0);
			
	}
}

double c_pupil_segmentation :: fit_ellipse( double & x, 
											double & y, 
											double & a, 
											double & b, 
											double & theta )
{
	double d = 10000000;
	for ( unsigned int i = 0; i < _nb_iter_ransac; ++ i )
	{
		unsigned int nb_points = 0;
		//Génération du masque
		unsigned int p = _nb_points_ransac;
		memset( _contour_mask, 0, sizeof(char) * _nb_points_contour );
		for ( unsigned int j = 0; j < _nb_points_contour && p > 0; ++ j )
		{
			unsigned int q;
			if ( _nb_points_contour == j )
				q = 0;
			else
				q = gsl_rng_uniform_int (rng, _nb_points_contour - j);
			if ( q < p )
			{
				_contour_mask[j] = 255;
				p --;
			}
		}
	
		//Construction du contour convex		
		_convex_boundary->compute( _contour,	
								   _nb_points_contour,
								   _contour_mask );
	   
		for ( unsigned int k = 0; k < _convex_boundary->nb_points(); ++ k )
		{
			unsigned label = 2 * _convex_boundary->labels()[k];
			_convex_contour[ 2 * k ] = _contour[ label ];
			_convex_contour[ 2 * k + 1 ] = _contour[ label + 1 ];
			
		}
		
		digitalize_euclide_contour(	_convex_contour_r,
									_nb_points_contour,
									_convex_contour,
									_convex_boundary->nb_points() ); 				
		//Estimation des paramètres de l'ellipse (1ere étape)
		ellipse.fit( _convex_contour_r,
					 _nb_points_contour,
					 1,
					 NULL );
		double tmp_r = 2 * _ransac_threshold * sqrt(0.5 * ( ellipse.a() * ellipse.a() + ellipse.b() * ellipse.b() ) );

		
		//Suppression des points trop loins du modèle
		for ( unsigned int i = 0; i < _nb_points_contour; ++ i )
		{
			if ( _contour_mask[ i ] )
			{
				nb_points ++;
				double 	dx, 
						dy, 
						mod, 
						angle, 
						r, 
						q;
				
				dx = _contour[2 * i] - ellipse.x_center();
				dy = _contour[2 * i + 1] - ellipse.y_center();
				mod = sqrt(dx * dx + dy * dy);
				angle = atan2( dy, dx );
				r = ellipse.a() * ellipse.b() / ( sqrt( pow(ellipse.a() * sin(angle + ellipse.theta()), 2) + pow(ellipse.b() * cos(angle + ellipse.theta()), 2) ) );
				q =  2 * abs( r - mod );
				if ( q > tmp_r )
				{
					_contour_mask[ i ] = 0;
					nb_points --;
				}
			}

		}
		

		if ( nb_points > 0.5 * _nb_points_ransac )
		{
			double v = 0;
			//Construction du contour convex		
			_convex_boundary->compute( _contour,	
									   _nb_points_contour,
									   _contour_mask );
									   
			for ( unsigned int k = 0; k < _convex_boundary->nb_points(); ++ k )
			{
				unsigned label = 2 * _convex_boundary->labels()[k];
				_convex_contour[ 2 * k ] = _contour[ label ];
				_convex_contour[ 2 * k + 1 ] = _contour[ label + 1 ];
				
			}
						
			digitalize_euclide_contour(	_convex_contour_r,
										_nb_points_contour,
										_convex_contour,
										_convex_boundary->nb_points() ); 				
						
			//Estimation des paramètres de l'ellipse (1ere étape)
			ellipse.fit( _convex_contour_r,
						 _nb_points_contour,
						 1,
						 NULL );
			
			//Calcul de la distance des points au modèle (Approximation)
			v = get_distance( ellipse.x_center(), 
							  ellipse.y_center(), 
							  ellipse.a(), 
							  ellipse.b(), 
							  ellipse.theta() );
			if ( v < d && 
				 ellipse.a() > _radius_min && 
				 ellipse.a() < _radius_max && 
				 ellipse.b() > _radius_min && 
				 ellipse.b() < _radius_max  )
			{

				d = v;
				x = ellipse.x_center();
				y = ellipse.y_center();
				a = ellipse.a();
				b = ellipse.b();
				theta = ellipse.theta();
			}
		}
	}
	return d;

}

double c_pupil_segmentation :: get_distance( 	const double & x, 
												const double & y, 
												const double & a, 
												const double & b, 
												const double & theta ) const
{
	double tmp_r = sqrt(0.5 * ( a * a + b * b ) );
	double d = 0;
	unsigned int nb_points = 0;
	for ( unsigned int i = 0; i < _nb_points_contour; ++ i )
	{
		if ( _contour_mask[ i ] )
		{
			nb_points ++;
			double 	dx = _contour[2 * i] - x,
					dy = _contour[2 * i + 1] - y,
					mod = sqrt(dx * dx + dy * dy),
					angle = atan2( dy, dx ),
					r = a * b / ( sqrt( pow(a * sin(angle + theta), 2) + pow(b * cos(angle + theta), 2) ) );
			double q =  2 * abs( r - mod );
			d += q;
		}

	}

	d /= nb_points * tmp_r;
	return d;
}

c_pupil_segmentation :: c_pupil_segmentation ( void )
{
	initialize();
}








c_pupil_segmentation :: c_pupil_segmentation (	unsigned int width, //Largeur
												unsigned int height, //Hauteur
												unsigned int window_size,
												double sigma_hist_image,
												double sigma_hist_diff_image,
												unsigned int pupil_limit_threshold,
												double radius_min,
												double radius_max,
												unsigned int nb_points_pupil, // Nombre de points du contour rééchantillonné de la pupille
												unsigned int nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
												unsigned int nb_points_ransac,
												unsigned int nb_iter_ransac,
												double ransac_threshold,
												ostream * _err_stream )
{
	initialize();
	if ( setup (	width, //Largeur
					height, //Hauteur
					window_size,
					sigma_hist_image,
					sigma_hist_diff_image,
					pupil_limit_threshold,
					radius_min,
					radius_max,
					nb_points_pupil, // Nombre de points du contour rééchantillonné de la pupille
					nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
					nb_points_ransac,
					nb_iter_ransac,
					ransac_threshold,
					_err_stream ) )
		throw ( invalid_argument("Arguments of c_pupil_segmentation!") );
}

int c_pupil_segmentation :: setup (	unsigned int width, //Largeur
									unsigned int height, //Hauteur
									unsigned int window_size,
									double sigma_hist_image,
									double sigma_hist_diff_image,
									unsigned int pupil_limit_threshold,
									double radius_min,
									double radius_max,
									unsigned int nb_points, // Nombre de points du contour rééchantillonné de la pupille
									unsigned int nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
									unsigned int nb_points_ransac,
									unsigned int nb_iter_ransac,
									double ransac_threshold,
									ostream * _err_stream )
{
	
	//Reset de l'objet
	free();
	initialize();
	
	err_stream = _err_stream;
	
	
	//Check des entrées
	if ( 	! width								||
			! height							||
			! window_size						||
			! pupil_limit_threshold				||
			! sigma_hist_image					||
			! sigma_hist_diff_image				||
			! radius_min						||
			! radius_max						||
			! nb_points							||
			! nb_iterations_ellipse				||
			! nb_points_ransac					||
			nb_points_ransac > nb_points		||
			! nb_iter_ransac 					)
	{
		if ( err_stream )
		{
			*err_stream << 
			"Error: Invalid argument in int c_pupil_segmentation :: setup :" << endl << 
			"( 	! width							||" << endl << 
			"! height							||" << endl << 
			"! window_size						||" << endl << 
			"! sigma_image_hist					||" << endl << 
			"! sigma_hist_diff_image			||" << endl << 
			"! radius_min						||" << endl << 
			"! radius_max						||" << endl << 
			"! nb_points_pupil					||" << endl << 
			"! nb_iteratons_ellipse				||" << endl << 
			"! nb_points_ransac					||" << endl << 
			"nb_points_ransac > nb_points_pupil	||" << endl << 
			"! nb_iter_ransac 					)" << endl;
		}
		return 1;
	}
	
	_width = width;
	_height = height;
	
	
	_window_size = window_size;
	_sigma_hist_1 = sigma_hist_image;
	_sigma_hist_2 = sigma_hist_diff_image;
	_radius_min = radius_min;
	_radius_max = radius_max;
	_pupil_limit_threshold = pupil_limit_threshold;
	_nb_points_contour = nb_points;
	_nb_iterations_ellipse = nb_iterations_ellipse;
	_nb_points_ransac = nb_points_ransac;
	_nb_iter_ransac = nb_iter_ransac;
	_ransac_threshold = ransac_threshold;
	
	alloc();
	return 0;
}

int  c_pupil_segmentation :: setup (	api_parameters & params,
										unsigned int width,
										unsigned int height,
										ostream * _err_stream,
										const char * n_space,
										const char * window_size_name,
										const char * sigma_image_name,
										const char * sigma_diff_name,
										const char * pupil_limit_threshold_name,
										const char * radius_min_name,
										const char * radius_max_name,
										const char * nb_points_name,
										const char * nb_iter_ellipse_name,
										const char * nb_points_ransac_name,
										const char * nb_iter_ransac_name,
										const char * ransac_threshold_name
										)
{

	int q = 0; //error
	stringstream oss;
	err_stream = _err_stream;
	
	unsigned int window_size;
	double sigma_hist_image;
	double sigma_hist_diff_image;
	double radius_min;
	double radius_max;
	unsigned int nb_points; // Nombre de points du contour rééchantillonné de la pupille
	unsigned int nb_iterations_ellipse; //Nombre d'itérations de l'algorithme de fitting d'ellipse 
	unsigned int nb_points_ransac;
	unsigned int nb_iter_ransac;
	
	double ransac_threshold;
	unsigned int pupil_limit_threshold;
	
	oss << n_space << "::" << window_size_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&window_size,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << nb_points_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&nb_points,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << nb_iter_ellipse_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&nb_iterations_ellipse,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << nb_points_ransac_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&nb_points_ransac,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << nb_iter_ransac_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&nb_iter_ransac,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << pupil_limit_threshold_name;
	if ( api_get_positive_integer( 	params,
									oss.str().c_str(),
									&pupil_limit_threshold,
									err_stream ) )
		q = 1;
	oss.str("");
	
	
	oss << n_space << "::" << sigma_image_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&sigma_hist_image,
							err_stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << sigma_diff_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&sigma_hist_diff_image,
							err_stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << radius_min_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&radius_min,
							err_stream	) )
		q = 1;
	oss.str("");

	oss << n_space << "::" << radius_max_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&radius_max,
							err_stream	) )
		q = 1;
	oss.str("");
		
	oss << n_space << "::" << ransac_threshold_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&ransac_threshold,
							err_stream	) )
		q = 1;
	oss.str("");
		
	if ( q )
		return 1;


	return setup (	width, //Largeur
					height, //Hauteur
					window_size,
					sigma_hist_image,
					sigma_hist_diff_image,
					pupil_limit_threshold,
					radius_min,
					radius_max,
					nb_points, // Nombre de points du contour rééchantillonné de la pupille
					nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
					nb_points_ransac,
					nb_iter_ransac,
					ransac_threshold,
					_err_stream );
}

void c_pupil_segmentation :: get_labels( )
{
	
	unsigned int end = 0, size;
	_nb_labels = 0;
	for (; end < nb_modes_image && m_image[end] < _pupil_limit_threshold; ++ end );

	//~ m_image[end] = _pupil_limit_threshold;
	//~ end ++;
	//~ m_d_image[nb_modes_d_image] = 255;
	//~ nb_modes_d_image ++;	

	size = end + nb_modes_d_image - 1;
	bool find = false;
	for ( unsigned int k = 0; k < size; ++ k )
	{
		unsigned int s1, s2;
		if ( k >= end )
		{
			for ( unsigned int l = k - end + 1; l < k && l < nb_modes_d_image; ++ l )
			{
				unsigned int i = k - l,
							 j = l;
							 
				labels_1[_nb_labels] = i;
				labels_2[_nb_labels] = j;
				_nb_labels ++;
			}
		}
		else
		{
			for ( unsigned int l = 0; l <= k && l < nb_modes_d_image; ++ l )
			{
				unsigned int i = k - l,
							 j = l;	 
				labels_1[_nb_labels] = i;
				labels_2[_nb_labels] = j;
				_nb_labels ++;
				
			}
		}



		
	}
	
	
	
}


int c_pupil_segmentation :: segment(	const unsigned char * img_data,
										unsigned int x,
										unsigned int y,
										unsigned int width,
										unsigned int height,
										unsigned int width_step)
{
	if ( width == 0 || 
		 height == 0 )
	{
		if ( err_stream )
			*err_stream << "Error: Tracking!" << endl;
		return -1;
		
	}

	if (	width > c_pupil_segmentation :: _width 	|| 
			height > c_pupil_segmentation :: _height 	)
	{
		if ( err_stream )
		{
			*err_stream << "Warning : Object re-allocation" << endl;
		}
		if ( setup( 	_width, //Largeur
						_height, //Hauteur
						_window_size,
						_sigma_hist_1,
						_sigma_hist_2,
						_pupil_limit_threshold,
						_radius_min,
						_radius_max,
						_nb_points_contour, // Nombre de points du contour rééchantillonné de la pupille
						_nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
						_nb_points_ransac,
						_nb_iter_ransac,
						_ransac_threshold ) )
			{
				if ( err_stream )
				{
					*err_stream << "ERROR: NOOB user!" << endl;
				}
				return 1;
			}
		
		
	}
	get_abs_diff_image( img_data,
						x,
						y, 
						width,
						height,
						width_step );
	compute_histograms( img_data,
						x,
						y, 
						width,
						height,
						width_step );
	
	mode_detection( width, height );
	get_labels();
	bool find = false;
	_score = 10000000;

	for ( unsigned int k = 0; k < _nb_labels && k < 5; ++ k )
	{

		mask(	img_data,
				x,
				y, 
				width,
				height,
				width_step,
				m_image[labels_1[k]], 
				m_d_image[labels_2[k]]  );
		if ( ! label() )
		{
			if ( ! erase_small_regions() )
			{
				fit_and_select_best_region();

			}
			
			
		}
	}
	if ( _score < _ransac_threshold )
	{
		find = true;
	}
	if (!find )
		return 1;
		
	//Masque de la pupille
	//Calcul du seuil de la pupil 90%
	_a += 2;
	_b += 2;
	_r = sqrt( (_a*_a + _b*_b) / 2 );
	
	
	memset( mask_1->imageData, 0, mask_1->widthStep * mask_1->height );
	cvEllipse( mask_1, cvPoint( _x, _y), cvSize(_a, _b), _theta / M_PI * 180, 0, 360, CV_RGB(255,255,255), CV_FILLED );
	hist_1.compute( img_data + y * width_step + x, 
				   (unsigned char*) mask_1->imageData, 
				   width, 
				   height, 
				   width_step, 
				   (unsigned int) mask_1->widthStep);
	_pupil_threshold = hist_1.get_percentile(85);

	if (_pupil_threshold > _pupil_limit_threshold )
		return 1;
	
	_x += x;
	_y += y;
	

	return 0;
}



int c_pupil_segmentation :: segment( const IplImage * image )
{
	unsigned int x,
				 y,
				 width,
				 height,
				 width_step;


	if ( image == NULL )
	{
		////*err_stream << "Error : image is NULL in int c_pupil_segmentation :: segment( const IplImage * image );" << endl;
		return -1;
	}

	if ( image->roi != NULL )
	{
		x = image->roi->xOffset;
		y = image->roi->yOffset;
		width = image->roi->width;
		height = image->roi->height;
	}
	else
	{
		x = 0;
		y = 0;
		width = image->width;
		height = image->height;
	}

	width_step = image->widthStep;


	switch( image->depth)
	{
		case( IPL_DEPTH_8U ):
			return segment(	(unsigned char *) image->imageData,
							x,
							y,
							width,
							height,
							width_step );
		break;
		default:
			if ( err_stream )
				*err_stream << "Error : No supported pixel format" << endl;
			return -1;
		break;
	}
}

c_pupil_segmentation :: ~c_pupil_segmentation()
{
	free();
	initialize();
}

void c_pupil_segmentation :: alloc()
{
	_abs_diff_image_8b = cvCreateImage( cvSize( _width, _height), IPL_DEPTH_8U, 1 );
	mask_1 = cvCreateImage( cvSize( _width, _height), IPL_DEPTH_8U, 1 );
	mask_2 = cvCreateImage( cvSize( _width, _height), IPL_DEPTH_8U, 1 );
	bg_regions = new c_label(_width, _height );
	connex_regions = new c_label(_width, _height );
	
	
	_valid_region_labels = new unsigned int[_width * _height ];
	memset( _valid_region_labels, 0, sizeof(int) * _width * _height  );
	
	
	contour_obj = new c_contour(_width * _height );
	
	
	_contour = new double[ 2 * _nb_points_contour + 2];
	memset( _contour, 0, sizeof(double) * (2 * _nb_points_contour + 2)  );
	
	
	_contour_mask = new unsigned char[ _nb_points_contour];
	memset( _contour_mask, 0, sizeof(unsigned char) * _nb_points_contour  );
	
	rng = gsl_rng_alloc (gsl_rng_taus);
	
	structuring_element_2 = 	
		cvCreateStructuringElementEx( 	2 * (_radius_min / 4) + 1,
										2 * (_radius_min / 4) + 1,
										(_radius_min / 4),
										(_radius_min / 4),
										CV_SHAPE_ELLIPSE,
										NULL);		
	
	structuring_element_1 = 	
		cvCreateStructuringElementEx( 	2 * (_radius_min / 8) + 1,
										2 * (_radius_min / 8) + 1,
										(_radius_min / 8),
										(_radius_min / 8),
										CV_SHAPE_ELLIPSE,
										NULL);		
										
										
	m_image = new unsigned int[256];
	memset( m_image, 0, sizeof(unsigned int) * 256  );
	
	m_d_image = new unsigned int[256];
	memset( m_d_image, 0, sizeof(unsigned int) * 256 );
	
	
	
	labels_1 = new unsigned int[65536];
	memset( labels_1, 0, sizeof(unsigned int) * 65536 );
	
	labels_2 = new unsigned int[65536];
	memset( labels_2, 0, sizeof(unsigned int) * 65536 );
									
	_convex_boundary = new convex_boundary<double>( _nb_points_contour );
	
	_convex_contour = new double[ _nb_points_contour * 2 + 2 ];
	memset( _convex_contour, 0, sizeof(double) * ( _nb_points_contour * 2 + 2 ) );
	
	_convex_contour_r = new double[ _nb_points_contour * 2 + 2];
	memset( _convex_contour_r, 0, sizeof(double) * ( _nb_points_contour * 2 + 2 ) );
}




void c_pupil_segmentation :: free()
{
	
	if (_abs_diff_image_8b)
		cvReleaseImage(&_abs_diff_image_8b);
	if (mask_1)
		cvReleaseImage(&mask_1);
		
	if (mask_2)
		cvReleaseImage(&mask_2);
		
	if ( bg_regions )
		delete bg_regions;
		
	if ( connex_regions )
		delete connex_regions;
	if ( _valid_region_labels )
		delete[] _valid_region_labels;
	if (contour_obj )
		delete contour_obj;
	if ( _contour )
		delete[] _contour;
	
	if ( _contour_mask )
		delete[] _contour_mask;
	if ( m_image )
		delete[] m_image;
		
	if ( m_d_image )
		delete[] m_d_image;
		
	if ( labels_1 )
		delete[] labels_1;
		
	if ( labels_2 )
		delete[] labels_2;
		
	if ( rng )
		gsl_rng_free(rng);
		
	if (structuring_element_1)
		cvReleaseStructuringElement(&structuring_element_1);

	if (structuring_element_2)
		cvReleaseStructuringElement(&structuring_element_2);	
	if ( _convex_boundary )
		delete _convex_boundary;
		
	if( _convex_contour )
		delete[] _convex_contour;
		
	if ( _convex_contour_r )
		delete[] _convex_contour_r;
}

void c_pupil_segmentation :: initialize()
{
	structuring_element_1 = 0;
	structuring_element_2 = 0;
	_width = 0;
	_height = 0;
	_abs_diff_image_8b = 0;
	_window_size = 0;
	_sigma_hist_1 = 0;
	_sigma_hist_2 = 0;
	nb_modes_image = 0;
	nb_modes_d_image = 0;
	_pupil_limit_threshold = 0;
	_nb_labels = 0;
	mask_1 = 0;
	mask_2 = 0;
	bg_regions = 0; 
	connex_regions = 0;
	_valid_region_labels = 0;
	_nb_valid_regions = 0;
	_radius_min = 0;
	_radius_max = 0;
	contour_obj = 0;
	_contour = 0;
	_nb_points_contour = 0;
	_contour_mask = 0;
	_a = 0;
	_b = 0;
	_x = 0;
	_y = 0;
	_theta = 0;
	_r = 0;
	_score = 0;
	_nb_points_ransac = 0;
	_nb_iter_ransac = 0;
	rng = 0;
	_nb_iterations_ellipse = 0;
	err_stream = 0;
	labels_1 = 0;
	labels_2 = 0;
	m_image = 0;
	m_d_image = 0;
	_convex_boundary = 0;
	_convex_contour = 0;
	_convex_contour_r = 0;
}
