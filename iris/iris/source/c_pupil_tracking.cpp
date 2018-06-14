#include "c_pupil_tracking.hpp"

c_pupil_tracking :: c_pupil_tracking  ()
: c_tracking()
{
	c_pupil_tracking :: initialize();
}

c_pupil_tracking :: c_pupil_tracking (	unsigned int image_width,
										unsigned int image_height,
										const gsl_vector * t_0,
										const gsl_matrix * sqrt_q_0,
										const gsl_matrix * f,
										const gsl_matrix * sqrt_q,
										const double & err_sigma,
										ostream * _err_stream )
{
	c_tracking :: initialize();
	c_pupil_tracking :: initialize();
	
	if ( setup (	image_width,
					image_height,
					t_0,
					sqrt_q_0,
					f,
					sqrt_q,
					err_sigma,
					_err_stream ) )
	{
		throw (invalid_argument("Bad args."));
	} 
}

int c_pupil_tracking :: setup (	unsigned int image_width,
								unsigned int image_height,
								const gsl_vector * t_0,
								const gsl_matrix * sqrt_q_0,
								const gsl_matrix * f,
								const gsl_matrix * sqrt_q,
								const double & err_sigma,
								ostream * _err_stream )
{
	unsigned int size_x;
	c_tracking :: free();
	c_tracking :: initialize();
	c_pupil_tracking :: initialize();
	err_stream = _err_stream;
	if ( !image_width || !image_height || !t_0 )
	{
		if ( err_stream )
			*err_stream << "Error : Bad argument in int c_pupil_tracking :: setup( unsigned int , unsigned int , const gsl_vector * , const gsl_matrix *, const gsl_matrix *, const gsl_matrix *, const double &);" << endl;
		return 1;
	}
	_image_width = image_width;
	_image_height = image_height;
	
	
	size_x = t_0->size - 3;
	
	int err = c_tracking :: setup ( t_0, sqrt_q_0, f, sqrt_q, size_x, err_sigma, err_stream );
	if ( err )
		return err;


	
	return 0;
}

int c_pupil_tracking :: setup (	api_parameters & params,
								unsigned int image_width,
								unsigned int image_height,
								ostream * _err_stream,								
								const char * n_space,
								const char * t_0_name,
								const char * sqrt_q_0_name,
								const char * f_name,
								const char * sqrt_q_name,
								const char * n_name)
{
	int q = 0;
	double n;
	gsl_vector t_0;
	gsl_matrix sqrt_q_0,
			   f,
			   sqrt_q;
	stringstream oss;
	err_stream = _err_stream;
	//T0
	oss << n_space << "::" << t_0_name;
	if ( api_get_vector(	params,
							oss.str().c_str(),
							&t_0,
							err_stream	) )
		q = 1;
	oss.str("");
	//sqrt_Q0
	oss << n_space << "::" << sqrt_q_0_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&sqrt_q_0,
							err_stream	) )
		q = 1;
	oss.str("");
	//F
	oss << n_space << "::" << f_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&f,
							err_stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << sqrt_q_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&sqrt_q,
							err_stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::"  << n_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&n,
							err_stream	) )
		q = 1;
	oss.str("");
	if ( q )
		return 1;
	
	return setup (	image_width,
					image_height,
					&t_0,
					&sqrt_q_0,
					&f,
					&sqrt_q,
					n,
					err_stream ) ;
}

void c_pupil_tracking :: predict_next_position( double current_x,
												double current_y,
												double current_r,
												unsigned int frame_id  )
{
	//Prediction à l'état n-1
	double tmp[3];
	tmp[0] = current_x - _image_width / 2;
	tmp[1] = current_y - _image_height / 2;
	tmp[2] = current_r;
	
	//Prédiction
	c_tracking :: predict_next_position( tmp, frame_id );
	
	//Calcul du roi
	double tmp_x_b, tmp_y_b,
		   tmp_x_e, tmp_y_e;

	tmp_x_b = _x_min->data[0] - _x_max->data[2] + _image_width / 2;
	tmp_y_b = _x_min->data[1] - _x_max->data[2] + _image_height / 2;
	
	tmp_x_e = _x_max->data[0] + _x_max->data[2] + _image_width / 2;
	tmp_y_e = _x_max->data[1] + _x_max->data[2] + _image_height / 2;
	
	if ( tmp_x_b < 0 )
		tmp_x_b = 0;
	
	if ( tmp_x_e < tmp_x_b )
		tmp_x_e = tmp_x_b;

	
	if ( tmp_y_b < 0 )
		tmp_y_b = 0;
	
	if ( tmp_y_e < tmp_y_b )
		tmp_y_e = tmp_y_b;


	if ( tmp_x_e > _image_width )
		tmp_x_e = _image_width;
		
	if ( tmp_x_b > tmp_x_e )
		tmp_x_b = tmp_x_e;
		

	if ( tmp_y_e > _image_height)
		tmp_y_e = _image_height;
		
	if ( tmp_y_b > tmp_y_e )
		tmp_y_b = tmp_y_e;

	_width_roi = tmp_x_e - tmp_x_b;
	_height_roi = tmp_y_e - tmp_y_b;

	_x_roi = tmp_x_b;
	_y_roi = tmp_y_b;
	
}

void c_pupil_tracking :: reset_tracking ( )
{
	_width_roi = _image_width;
	_height_roi =_image_height;

	_x_roi = 0;
	_y_roi = 0;
	c_tracking::reset_tracking();
	
}
void c_pupil_tracking :: reset( )
{
	c_tracking:: reset( );
	_width_roi = _image_width;
	_height_roi =_image_height;

	_x_roi = 0;
	_y_roi = 0;
	
	
}


void c_pupil_tracking :: initialize ( )
{
	_image_width = 0;
	_image_height = 0;
	_x_roi = 0;
	_y_roi = 0;
	_width_roi = 0;
	_height_roi = 0;
}


