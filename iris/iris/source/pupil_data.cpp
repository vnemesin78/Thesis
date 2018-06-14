#include "pupil_data.hpp"
pupil_data :: pupil_data()
{
	initialize();
}

pupil_data :: pupil_data( const CvSize & image_size )
{
	initialize();
	setup ( image_size );
}

pupil_data :: pupil_data( const pupil_data & data )
{
	initialize();
	setup( data );	
}

void pupil_data :: initialize()
{
	x_pupil = 0;
	y_pupil = 0;
	score = 0;
	a_pupil = 0;
	b_pupil = 0;
	
	pupil_threshold = 0;
	smoothed_image = 0;
	
	seg_ok = false;
	theta_pupil = 0;
	roi.x = 0;
	roi.y = 0;
	roi.width = 0;
	roi.height = 0;
	
}

void pupil_data :: free()
{
	if ( smoothed_image )
	{
		cvReleaseImage( &smoothed_image );
	}
}

void pupil_data :: alloc()
{
	if ( _img_data.width && _img_data.height )
	{
								
		smoothed_image = cvCreateImage( 	cvSize( _img_data.width,
													_img_data.height ),
												IPL_DEPTH_8U,
												1 );
	}
}

int pupil_data :: setup()
{
	free();
	initialize();
	return 0;
}
int pupil_data :: set_dim( unsigned int w, unsigned int h)
{
	if ( w > _img_data.width || h > _img_data.height )
		setup( cvSize( w, h ) );
	cvSetImageROI( 	smoothed_image,
					cvRect(0,0, w, h ) );	
	return 0;
}



int pupil_data :: setup ( const CvSize & image_size )
{
	setup();

	if ( _img_data.setup( image_size ) )
		return 1;
	alloc();
	return 0;
}

int pupil_data :: setup ( const pupil_data & data )
{
	CvSize tmp;
	setup();

	tmp = cvSize( data._img_data.width, data._img_data.height );
	if ( _img_data.setup( tmp ) )
		return 1;
	_img_data = data._img_data;
	alloc();
	cvSetImageROI ( smoothed_image,
					cvRect( 0,
							0,
							GET_IMAGE_WIDTH(data.smoothed_image),
							GET_IMAGE_HEIGHT(data.smoothed_image) ) );
	cvCopyImage ( data.smoothed_image, smoothed_image );
	
	
	
	x_pupil = data.x_pupil;
	y_pupil = data.y_pupil;
	a_pupil = data.a_pupil;
	b_pupil = data.b_pupil;
	score = data.score;
	pupil_threshold = data.pupil_threshold;
	
	seg_ok = data.seg_ok;
	theta_pupil = data.theta_pupil;
	roi = data.roi;		
	return 0;
}

pupil_data & pupil_data :: operator=( const pupil_data & _p_data )
{
	if ( this == &_p_data )
		return *this;
	
	if ( ! ( (*this) == _p_data ) )
		setup( _p_data );

	_img_data = _p_data._img_data;
	cvSetImageROI ( smoothed_image,
					cvRect( 0,
							0,
							GET_IMAGE_WIDTH(_p_data.smoothed_image),
							GET_IMAGE_HEIGHT(_p_data.smoothed_image) ) );
	cvCopyImage ( _p_data.smoothed_image, smoothed_image );
	
	x_pupil = _p_data.x_pupil;
	y_pupil = _p_data.y_pupil;
	a_pupil = _p_data.a_pupil;
	b_pupil = _p_data.b_pupil;
	score = _p_data.score;
	theta_pupil = _p_data.theta_pupil;
	pupil_threshold = _p_data.pupil_threshold;
	roi = _p_data.roi;
	seg_ok = _p_data.seg_ok;
	return *this;
}

pupil_data :: ~pupil_data()
{
	free();
	initialize();
}

bool pupil_data :: operator==( const pupil_data & data ) const
{
	return ( data._img_data == _img_data );
}

int pupil_data :: save ( const char * rep,
						   unsigned int id ) const
{
	int q = 0;
	stringstream oss;
	
	if ( _img_data.save( rep, id ) )
		return 1;
		
	if ( seg_ok == false )
		return 1;
		
	//Données de seg.
	{
		oss << rep << "/" << "pupil_data_" << id << ".txt";
		ofstream file(oss.str().c_str());
		oss.str("");
		
		file << "#Pupil seg. data" << endl;
		file << "x_pupil = " << x_pupil << endl;
		file << "y_pupil = " << y_pupil << endl;		
		file << "a_pupil = " << a_pupil << endl;
		file << "b_pupil = " << b_pupil << endl;		
		file << "theta_pupil = " << theta_pupil << endl;
		file << "threshold = " << pupil_threshold << endl;
		file << "x_roi = " << roi.x << endl;
		file << "y_roi = " << roi.y << endl;
		file << "w_roi = " << roi.width << endl;
		file << "h_roi = " << roi.height << endl;
		file << "score = " << score << endl;
		file.close();	
	}	
	
	{
		
		oss << rep << "/" << "smoothed_image_" << id << ".png";
		cvSaveImage ( oss.str().c_str(),
					  smoothed_image );
		oss.str("");
		
	}
	
	return q;
}

int pupil_data :: load ( const char * rep,
						   unsigned int id )
{
	free();
	initialize();
	
	int q = 0;
	stringstream oss;
	
	if ( _img_data.load( rep, id ) )
		q = 1;
	//Données de seg.
	{
		oss << rep << "/" << "pupil_data_" << id << ".txt";
			api_parameters params;
			if ( params.load( oss.str().c_str() ) )
				q = 1;	
			else
			{
				if ( api_get_integer( params,
									  "x_roi",
									  &roi.x ) )
					q = 1;
				if ( api_get_integer( params,
									  "y_roi",
									  &roi.y ) )
					q = 1;					
				if ( api_get_integer( params,
									  "w_roi",
									  &roi.width ) )
					q = 1;
				if ( api_get_integer( params,
									  "h_roi",
									  &roi.height ) )
					q = 1;					
				if ( api_get_double ( params,
									  "x_pupil",
									  &x_pupil ) )
					q = 1;
				if ( api_get_double ( params,
									  "y_pupil",
									  &y_pupil ) )
					q = 1;					
				if ( api_get_double ( params,
									  "a_pupil",
									  &a_pupil ) )
					q = 1;					
				if ( api_get_double ( params,
									  "score",
									  &score ) )
					q = 1;	
				if ( api_get_double ( params,
									  "b_pupil",
									  &b_pupil ) )
					q = 1;					
				if ( api_get_double ( params,
									  "theta_pupil",
									  &theta_pupil ) )
					q = 1;					
				if ( api_get_double ( params,
									  "threshold",
									  &pupil_threshold ) )
					q = 1;				
			}

		oss.str("");
	}	
	
	oss << rep << "/" << "smoothed_image_" << id << ".png";
	smoothed_image = cvLoadImage( oss.str().c_str(), 0 );
	oss.str("");
	
	if ( ! smoothed_image )
		q = 1;
		
	if ( ! q ) 
		seg_ok = true;
	
	return q;
}


void * pupil_data_alloc( const void * params )
{
	CvSize * data = (CvSize*) params;
	pupil_data * obj = new pupil_data( *data );
	return (void*) obj;
}

void pupil_data_copy( 	void * data_tg,
						const void * data_src )
{
	pupil_data * tg = (pupil_data*) data_tg,
				* src = (pupil_data*) data_src;
				
	(*tg) = (*src);
}

void pupil_data_free( void * p )
{
	pupil_data * tg = (pupil_data*) p;
	delete tg;
}

