#include "c_pupil_thread.hpp"

c_pupil_thread :: c_pupil_thread()
{
	initialize();
}

int c_pupil_thread :: setup ( void )
{
	free();
	initialize();
	return 0;
}

int c_pupil_thread :: setup (	char ** argv,	
								unsigned int argc,
								ostream * _err_stream)
{
	setup();
	err_stream = _err_stream;
	
	
	if ( argc == 0 )
	{
		if ( err_stream )
			*err_stream << "Error : Missing argument(s)" << endl;
		return 1;
	}
	
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
							"tracking",
							&tracking ) )
		tracking = 0;
	
	preprocessing = new c_preprocessing;
	pupil_tracking = new c_pupil_tracking;
	pupil_segmentation = new c_pupil_segmentation;
	p_data = new pupil_data;	
	
	
	
	if ( preprocessing->setup( 	params,		
								d_width, 
								d_height, 
								_err_stream ) )
		q = 1;
	if (	pupil_tracking->setup ( 	params,
										d_width,
										d_height,
										err_stream) )
		q = 1;
	if ( pupil_segmentation->setup (	params,
										d_width,
										d_height,
										err_stream ) )
		q = 1;
	
	p_data->setup( cvSize ( d_width, d_height ) );
		
	if ( q )
		return 1;
	pupil_tracking->reset_tracking();
	return 0;
}

int c_pupil_thread :: segment_pupil(	const image_data & img_data )
{
	p_data->_img_data = img_data;	
	p_data->seg_ok = false;
	if ( ! img_data.img_ok )
	{
		pupil_tracking->reset_tracking();
		return 1;		
	}
	const IplImage * img_src = img_data.image;
	
	//Preprocessing
	if ( preprocessing->process( img_src ) )
	{
		if (err_stream )
		{
			*err_stream << "Error: Unknown error!" << endl;
			return 1;
		}	
	}

	
	cvSetImageROI( 	p_data->smoothed_image, 
					cvRect( 0,
							0, 
							GET_IMAGE_WIDTH(preprocessing->smoothed_image()), 
							GET_IMAGE_HEIGHT(preprocessing->smoothed_image() ) ) );
	cvCopyImage( preprocessing->smoothed_image(), p_data->smoothed_image );
	//Segmentation
	int	x_s = MAX( pupil_tracking->x_roi(),  GET_IMAGE_X_OFFSET( preprocessing->smoothed_image() ) ),
		y_s = MAX( pupil_tracking->y_roi(),  GET_IMAGE_Y_OFFSET( preprocessing->smoothed_image() ) ),
		x_e = MIN( pupil_tracking->x_roi() + pupil_tracking->width_roi(), GET_IMAGE_X_OFFSET( preprocessing->smoothed_image() ) + GET_IMAGE_WIDTH(preprocessing->smoothed_image()) ),
		y_e = MIN( pupil_tracking->y_roi() + pupil_tracking->height_roi(), GET_IMAGE_Y_OFFSET( preprocessing->smoothed_image() ) + GET_IMAGE_HEIGHT(preprocessing->smoothed_image()) );

	if ( ! pupil_segmentation->segment( 	(unsigned char*) (preprocessing->smoothed_image()->imageData),
											x_s,
											y_s,
											x_e - x_s,
											y_e - y_s,
											preprocessing->smoothed_image()->widthStep / sizeof(unsigned char) ) )
	{
		//Données de segmentation
		p_data->seg_ok = true;
		//Tracking
		if ( tracking )
		{
			double r = sqrt((pow(pupil_segmentation->a(),2) + pow(pupil_segmentation->b(),2) )/2);

			pupil_tracking->predict_next_position( 	pupil_segmentation->x(), 
													pupil_segmentation->y(), 
													r,
													img_data.frame_id );				
		

		}
		if ( p_data->seg_ok )
		{
			
			//Données de segmentation et de tracking
			p_data->x_pupil = pupil_segmentation->x();
			p_data->y_pupil = pupil_segmentation->y();
			p_data->a_pupil = pupil_segmentation->a();
			p_data->b_pupil = pupil_segmentation->b();
			p_data->theta_pupil = pupil_segmentation->theta();
			p_data->roi = cvRect ( 	pupil_tracking->x_roi(), 
									pupil_tracking->y_roi(), 
									pupil_tracking->width_roi(),
									pupil_tracking->height_roi() );
			p_data->score = pupil_segmentation->score();
			//Moyenne / Variance
			p_data->pupil_threshold = pupil_segmentation->pupil_threshold();

		}
	}
	else
	{

		pupil_tracking->reset_tracking();
	}
		
	return (!p_data->seg_ok );
}

c_pupil_thread :: ~c_pupil_thread()
{
	free();
	initialize();
}


void c_pupil_thread :: free()
{
	if ( pupil_tracking )	
		delete pupil_tracking;
	if ( pupil_segmentation)
		delete pupil_segmentation;
	if ( p_data )
		delete p_data;
	if ( preprocessing ) 
		delete preprocessing;
}

void c_pupil_thread :: initialize()
{
	err_stream = NULL;
	pupil_tracking = NULL;
	pupil_segmentation = NULL;
	p_data = NULL;
	preprocessing = NULL;
}

void c_pupil_thread :: reset ()
{
	pupil_tracking->reset();
	//~ pupil_segmentation->reset();
	p_data->seg_ok = false;
}

