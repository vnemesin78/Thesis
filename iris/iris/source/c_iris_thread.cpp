#include "c_iris_thread.hpp"
#include <cmath>



void draw_curve( 	IplImage * image, 
					double t_min, 
					double t_max, 
					unsigned int nb_t, 
					const function_2d & f, 
					const double * params,
					unsigned char r, 
					unsigned char g, 
					unsigned char b )
{		
	
	double dt;
	if ( nb_t == 0 )
	{
		dt = 0;
	}
	else
		dt = (t_max -  t_min ) / nb_t;
	
	for (unsigned int i = 0; i < nb_t; ++ i )
	{
		double x,y,t;
		t = i * dt + t_min;
		f.get_value( x,y, t, params );
		cvCircle( image, 
				  cvPoint(x,y),
				  2,
				  CV_RGB(r,g,b),
				  CV_FILLED );
	}
		
	
	
	
}


c_iris_thread :: c_iris_thread()
{
	initialize();
}

int c_iris_thread :: setup ( void )
{
	free();
	initialize();
	return 0;
}


int c_iris_thread :: setup( 	char ** argv,
								unsigned int argc,
								ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	if ( argc == 0 )
	{
		if (err_stream )
			*err_stream << "Error : Missing argument(s)" << endl;
		return 1;
	}

	api_parameters params;
	for ( unsigned int i = 0; i < argc; ++ i )
		params.load(argv[i]);

	unsigned int 	d_width,
					d_height;
	int q = 0;
	if ( api_get_positive_integer( 	params,
									"width",
									&d_width ) )
		d_width = DEFAULT_WIDTH;

	if ( api_get_positive_integer( 	params,
									"height",
									&d_height) )
		d_height = DEFAULT_HEIGHT;

	if ( api_get_integer( 	params,
							"score::type",
							&_score_type ) )
		_score_type = 0;
	if ( api_get_double( 	params,
							"iris::validness_factor",
							&_validness_factor,
							err_stream) )
		q = 1;	


	if ( api_get_double( 	params,
							"iris::spot_threshold",
							&_spot_threshold,
							err_stream ) )
		q = 1;


	if ( api_get_double( 	params,
							"iris::n_sigma_up",
							&_n_sigma_up,
							err_stream ) )
		q = 1;


	if ( api_get_double( 	params,
							"iris::n_sigma_low",
							&_n_sigma_low,
							err_stream ) )
		q = 1;

	unsigned int opening, eroding;
	
	if ( api_get_positive_integer(	params, "iris::opening", &opening, err_stream ) )
		q = 1;
		
	if ( api_get_positive_integer(	params, "iris::eroding", &eroding, err_stream ) )
		q = 1;
		
	structuring_element_1 = 
		cvCreateStructuringElementEx( 	opening * 2 + 1,
										opening * 2 + 1,
										opening,
										opening,
										CV_SHAPE_ELLIPSE,
										NULL);

	structuring_element_2 = 
		cvCreateStructuringElementEx( 	eroding * 2 + 1,
										eroding * 2 + 1,
										eroding,
										eroding,
										CV_SHAPE_ELLIPSE,
										NULL);




	iris_segmentation = new c_iris_segmentation;
	eyelid_obj = new c_eyelid_segmentation;
	polar_iris = new c_polar_iris;
	eyelid_segmentation = new c_eyelids_segmentation;
	iris_code = new c_iris_code;
	i_data = new iris_data;



	if ( iris_segmentation->setup( params, d_width, d_height, err_stream ) )
		q = 1;

	if ( polar_iris->setup( params, d_width, d_height, err_stream ) )
		q = 1;
	if ( eyelid_segmentation->setup( params,
									d_width * (1 + 2 * iris_segmentation->r_ratio() ),
									d_height * (1 + 2 * iris_segmentation->r_ratio() ), err_stream ) )
		q = 1;

	if ( iris_code->setup( params, err_stream ) )
		q = 1;

	if ( _score_type == 0 )
	{
		focus_score = new c_focus_score;
		if (  focus_score->setup( params, err_stream ) )
			q = 1;
	}
	
	if ( eyelid_obj->setup( d_width, d_height, params, err_stream ) )
		q = 1;


	iris_data_params args;
		args.img_width = d_width;
		args.img_height = d_height;
		args.polar_width = polar_iris->nb_directions();
		args.polar_height = polar_iris->nb_samples();
		args.nb_samples_iris = polar_iris->nb_samples_iris();
		args.iris_code_width = iris_code->nb_directions();
		args.iris_code_height = iris_code->nb_samples();
		args.iris_width = eyelid_segmentation->width();
		args.iris_height = eyelid_segmentation->height();

	if (i_data->setup( args ))
		q = 1;

	tmp_p_img = cvCreateImage( 	cvSize ( 	iris_code->nb_directions(),
											iris_code->nb_samples() ),
								IPL_DEPTH_64F,
								1 );

	tmp_p_mask = cvCreateImage( cvSize ( 	iris_code->nb_directions(),
											iris_code->nb_samples() ),
								IPL_DEPTH_8U,
								1 );

	label_obj = new c_label ( d_width, d_height );

	if ( q )
		return 1;

	return 0;
}

int c_iris_thread :: segment_iris( const pupil_data & p_data )
{

		i_data->seg_ok = false;
		i_data->p_data = p_data;

		if ( ! p_data.seg_ok )
			return 1;
		if ( !mask_1 || !mask_2 || ! mask_3 )
		{
			mask_1 = cvCreateImage( cvSize( p_data._img_data.width, p_data._img_data.height ),
											IPL_DEPTH_8U,
											1 ); 
			mask_2 = cvCreateImage( cvSize( p_data._img_data.width, p_data._img_data.height ),
											IPL_DEPTH_8U,
											1 );
			mask_3 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
											GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
											IPL_DEPTH_8U,
											1 ); 			
											
		}
		else if ( 	GET_IMAGE_WIDTH(mask_1) < p_data._img_data.width ||
					GET_IMAGE_HEIGHT(mask_2) < p_data._img_data.height )
		{
			cvReleaseImage( &mask_1 );
			cvReleaseImage( &mask_2 );
			cvReleaseImage( &mask_3 );
			mask_1 = cvCreateImage( cvSize( p_data._img_data.width, p_data._img_data.height ),
											IPL_DEPTH_8U,
											1 ); 
			mask_2 = cvCreateImage( cvSize( p_data._img_data.width, p_data._img_data.height ),
											IPL_DEPTH_8U,
											1 ); 
			mask_3 = cvCreateImage( cvSize( p_data._img_data.width, p_data._img_data.height ),
											IPL_DEPTH_8U,
											1 ); 					
											
		}

		cvSetImageROI( mask_1, cvRect( 0, 0, p_data._img_data.width, p_data._img_data.height ) );
		memset( mask_1->imageData, 255, mask_1->widthStep * mask_1->height );

		//Op. intégro diff.


		if ( ! iris_segmentation->segment_iris( 	p_data._img_data.image,
													p_data.smoothed_image,
													mask_1,
													mask_1,
													p_data.x_pupil,
													p_data.y_pupil,
													p_data.a_pupil,
													p_data.b_pupil,
													p_data.theta_pupil ) )
		{
			//Seuillage /p_threhsold et /spot_threshold
			if ( !mask_1 || !mask_2 || ! mask_3)
			{
				mask_1 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 
				mask_2 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 
				mask_3 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 								
												
												
			}
			else if ( 	GET_IMAGE_WIDTH(mask_1) < GET_IMAGE_WIDTH(iris_segmentation->new_image() ) ||
						GET_IMAGE_HEIGHT(mask_2) < GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) )
			{
				cvReleaseImage( &mask_1 );
				cvReleaseImage( &mask_2 );
				cvReleaseImage( &mask_3 );
				mask_1 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 
				mask_2 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 
				mask_3 = cvCreateImage( cvSize( GET_IMAGE_WIDTH(iris_segmentation->new_image() ),
												GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ),
												IPL_DEPTH_8U,
												1 ); 	
			}
			cvSetImageROI( mask_1, cvRect( 0, 0, GET_IMAGE_WIDTH(iris_segmentation->new_image() ), GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ) );
			cvSetImageROI( mask_2, cvRect( 0, 0, GET_IMAGE_WIDTH(iris_segmentation->new_image() ), GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ) );
			cvSetImageROI( mask_3, cvRect( 0, 0, GET_IMAGE_WIDTH(iris_segmentation->new_image() ), GET_IMAGE_HEIGHT(iris_segmentation->new_image() ) ) );	
			
			double 	s_1 = 1.00 * p_data.pupil_threshold,
					s_2 = _spot_threshold * get_image_mean( iris_segmentation->new_image() ) / 128.0,
					s_3 = _spot_threshold * get_image_mean( p_data.smoothed_image ) / 128.0,
					s_4 = s_1 * s_2 / s_3;
					
			memset(	mask_2->imageData, 
					255, 
					mask_2->height * mask_2->widthStep );	
			
			cvEllipse( 	mask_2, 
						cvPoint( iris_segmentation->new_x_pupil(), iris_segmentation->new_y_pupil() ),
						cvSize( iris_segmentation->new_a_pupil() + 2, iris_segmentation->new_b_pupil() + 2 ),
						iris_segmentation->new_theta_pupil() / M_PI * 180.0,
						0,
						360,
						CV_RGB(0,0,0),
						CV_FILLED);
	
			eyelid_obj->compute( 	p_data.smoothed_image, 
									mask_2, 
									iris_segmentation->new_x(), 
									iris_segmentation->new_y(), 
									iris_segmentation->new_r(),
									iris_segmentation->new_x_pupil(), 
									iris_segmentation->new_y_pupil(), 
									iris_segmentation->new_a_pupil(), 
									iris_segmentation->new_b_pupil() );
									
			memset(	mask_2->imageData, 
					0, 
					mask_2->height * mask_2->widthStep );	
					
			cvCircle( 	mask_2, 
						cvPoint( iris_segmentation->new_x(), iris_segmentation->new_y() ),
						iris_segmentation->new_r(),
						CV_RGB(255,255,255),
						CV_FILLED);
			cvEllipse( 	mask_2, 
						cvPoint( iris_segmentation->new_x_pupil(), iris_segmentation->new_y_pupil() ),
						cvSize( iris_segmentation->new_a_pupil() + 2, iris_segmentation->new_b_pupil() + 2 ),
						iris_segmentation->new_theta_pupil() / M_PI * 180.0,
						0,
						360,
						CV_RGB(0,0,0),
						CV_FILLED);		
					
			//~ //Seuillage des spots
			cvThreshold( iris_segmentation->new_image(),
						 mask_1,
						 s_2,
						 255,
						 CV_THRESH_BINARY_INV );
			cvAnd( mask_2, mask_1, mask_1);
						 
						 
			cvThreshold( iris_segmentation->new_image(),
						 mask_2,
						 s_4,
						 255,
						 CV_THRESH_BINARY);	 
						 
			cvAnd( mask_2, eyelid_obj->mask(), mask_2 );
			cvAnd( mask_1, mask_2, mask_1 );

			
			//Transformée polaire normalisée
			if (! polar_iris->compute(	iris_segmentation->new_image(),
										mask_1,
										iris_segmentation->new_x_pupil(),
										iris_segmentation->new_y_pupil(),
										iris_segmentation->new_a_pupil(),
										iris_segmentation->new_b_pupil(),
										iris_segmentation->new_theta_pupil(),
										iris_segmentation->new_x(),
										iris_segmentation->new_y(),
										iris_segmentation->new_r() ) )
			{
				{
					image_resize ( 	(double *) tmp_p_img->imageData,
									(unsigned char *) tmp_p_mask->imageData,
									iris_code->nb_directions(),
									iris_code->nb_samples(),
									tmp_p_img->widthStep / sizeof(double),
									tmp_p_mask->widthStep / sizeof(char),
									(double *) polar_iris->polar_image()->imageData,
									(unsigned char *) polar_iris->polar_mask()->imageData,
									polar_iris->nb_directions(),
									polar_iris->nb_samples_iris(),
									polar_iris->polar_image()->widthStep / sizeof(double),
									polar_iris->polar_mask()->widthStep / sizeof(char) );

					if ( !iris_code->compute_iris_code( 	tmp_p_img,
															tmp_p_mask ) )
					{
						if ( get_image_mean ( iris_code->mask() ) > _validness_factor * 255 )
						{

							i_data->seg_ok = true;
							i_data->x_iris = iris_segmentation->x_ellipse();
							i_data->y_iris = iris_segmentation->y_ellipse();
							i_data->a_iris = iris_segmentation->a_ellipse();
							i_data->b_iris = iris_segmentation->b_ellipse();
							i_data->theta_iris = iris_segmentation->theta_ellipse();
							i_data->new_x_iris = iris_segmentation->new_x();
							i_data->new_y_iris = iris_segmentation->new_y();
							i_data->new_r_iris = iris_segmentation->new_r();
							
							i_data->a_upper = eyelid_obj->a_upper();
							i_data->c_upper = eyelid_obj->c_upper();
							i_data->theta_upper = eyelid_obj->theta_upper();
							
							i_data->a_lower = eyelid_obj->a_lower();
							i_data->c_lower = eyelid_obj->c_lower();
							i_data->theta_lower = eyelid_obj->theta_lower();
							

							cvSetImageROI( 	i_data->polar_image,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( polar_iris->polar_image() ),
													GET_IMAGE_HEIGHT( polar_iris->polar_image() ) ) );

							cvCopyImage( 	polar_iris->polar_image(),
											i_data->polar_image );

							cvSetImageROI( 	i_data->polar_mask,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( polar_iris->polar_mask() ),
													GET_IMAGE_HEIGHT( polar_iris->polar_mask() ) ) );
							cvCopyImage( 	polar_iris->polar_mask(),
											i_data->polar_mask );

							cvSetImageROI( 	i_data->iris_image,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( iris_segmentation->new_image() ),
													GET_IMAGE_HEIGHT( iris_segmentation->new_image() ) ) );

							cvCopyImage( 	iris_segmentation->new_image(),
											i_data->iris_image );

							cvSetImageROI( 	i_data->iris_mask,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( mask_1 ),
													GET_IMAGE_HEIGHT( mask_1 ) ) );
							cvCopyImage( 	mask_1,
											i_data->iris_mask );

							cvSetImageROI( 	i_data->code,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( iris_code->code() ),
													2 * GET_IMAGE_HEIGHT( iris_code->code()  ) ) );
							cvCopyImage( 	iris_code->code(),
											i_data->code );

							cvSetImageROI( 	i_data->code_mask,
											cvRect(	0,
													0,
													GET_IMAGE_WIDTH( iris_code->mask() ),
													GET_IMAGE_HEIGHT( iris_code->mask() ) ) );
							cvCopyImage( 	iris_code->mask(),
											i_data->code_mask );
							
							//Calcul du score
							switch ( _score_type )
							{
								case(1):
									i_data->nrj_ratio = iris_code->nrj_ratio();
								break;
								default:
								{
									CvRect rect = 
										cvRect( 	
											i_data->new_x_iris - i_data->new_r_iris,
											i_data->new_y_iris - i_data->new_r_iris,
											2 * i_data->new_r_iris,
											2 * i_data->new_r_iris );

									if ( 	! ( 	rect.x + rect.width > i_data->iris_image->width ||
											rect.y + rect.height > i_data->iris_image->height 					) )
										i_data->nrj_ratio = 
											focus_score->get_score (
												i_data->iris_image,
												i_data->iris_mask,
												rect,
												rect.width ,
												rect.height );
								}
								break;
							}


						}
						else
						{
							if (err_stream )
								*err_stream << "Warning : Iris code failed : Not enough pixels in iris code (" << get_image_mean ( iris_code->mask() ) / 255.0 * 100 <<"/" << _validness_factor * 100 << ")" << endl;
						}
					}
				}
			}

		}

	return (!i_data->seg_ok );
}

c_iris_thread :: ~c_iris_thread()
{
	free();
	initialize();
}

int c_iris_thread :: get_data ( iris_data & data )
{
		data = *i_data;
	return 0;
}

void c_iris_thread :: free()
{
	if ( tmp_p_img )
		cvReleaseImage ( &tmp_p_img );
	if ( tmp_p_mask )
		cvReleaseImage ( &tmp_p_mask );
	if ( iris_segmentation )
		delete iris_segmentation;
	if ( eyelid_segmentation )
		delete eyelid_segmentation;
	if ( iris_code )
		delete iris_code;
	if (polar_iris )
		delete polar_iris;
	if ( i_data )
		delete i_data;
	if ( focus_score )
		delete focus_score;
	if ( mask_1 )
		cvReleaseImage(&mask_1);
	if ( mask_2 )
		cvReleaseImage(&mask_2);
	if ( mask_3 )
		cvReleaseImage(&mask_3);
	if ( structuring_element_1 )
		cvReleaseStructuringElement( & structuring_element_1 );
	if ( structuring_element_2)
		cvReleaseStructuringElement( & structuring_element_2 );	
	if ( label_obj )
		delete label_obj;
	if (eyelid_obj )
		delete eyelid_obj;
}

void c_iris_thread :: initialize()
{
	tmp_p_img = 0;
	tmp_p_mask = 0;
	structuring_element_1 = 0;
	label_obj = 0;
	structuring_element_2 = 0;
	_score_type = 0;
	_validness_factor = 0;
	iris_segmentation = 0;
	eyelid_segmentation = 0;
	polar_iris = 0;
	iris_code = 0;
	i_data = 0;
	_n_sigma_low = 0;
	_n_sigma_up = 0;
	err_stream = NULL;
	focus_score = 0;
	_spot_threshold = 0;
	mask_1 = 0;
	mask_2 = 0;
	mask_3 = 0;
	eyelid_obj = 0;
}

void c_iris_thread :: reset()
{
	i_data->seg_ok = false;
}







