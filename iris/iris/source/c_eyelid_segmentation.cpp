#include "c_eyelid_segmentation.hpp"
#include <sstream>
using namespace std;
double x_par( 	const double & x,
				const double * params )
{
	return x * cos(params[5]) + ( params[2] * pow( ( x - params[3] ), 2 ) + params[4] ) * sin ( params[5] ) + params[0];
}

double y_par( 	const double & x,
				const double * params )
{
	return -x * sin(params[5]) + ( params[2] * pow( ( x - params[3] ), 2 ) + params[4] ) * cos ( params[5] ) + params[1];
}

double dx_par( 	const double & x,
				const double * params )
{
	return abs( cos(params[5]) + ( params[2] *  2 * ( x - params[3] ) ) * sin ( params[5] ) );
}

double dy_par( 	const double & x,
				const double * params )
{
	return abs(- sin(params[5]) + ( params[2] *  2 * ( x - params[3] ) ) * cos ( params[5] ) );
}

function_2d parabol2(  x_par, y_par);
function_2d d_parabol2(  x_par, y_par);



c_eyelid_segmentation :: c_eyelid_segmentation ( void )
{
	initialize();
}	

c_eyelid_segmentation :: c_eyelid_segmentation ( 	unsigned int width,
													unsigned int height,
													unsigned int nb_t,
													unsigned int nb_a,
													unsigned int nb_c,
													unsigned int nb_theta,
													double theta_min,
													double theta_max,
													unsigned int kernel_size,
													ostream * stream)
{
	initialize();
	
	setup(  width,
			height,
			nb_t,
			nb_a,
			nb_c,
			nb_theta,
			theta_min,
			theta_max,
			kernel_size,
			stream );
}

int c_eyelid_segmentation :: setup ( 	unsigned int width,
										unsigned int height,
										unsigned int nb_t,
										unsigned int nb_a,
										unsigned int nb_c,
										unsigned int nb_theta,
										double theta_min,
										double theta_max,
										unsigned int kernel_size,
										ostream * stream )
{
	free();
	initialize();
	err_stream = stream;
	
	if ( width == 0 || height == 0 || nb_a == 0 || nb_c == 0 || nb_theta == 0 || nb_t == 0 || kernel_size == 0)
	{
		if ( err_stream )
			*err_stream << "Error : ( width == 0 || height == 0 || nb_a == 0 || nb_c == 0 || nb_theta == 0 ) in  c_eyelid_segmentation :: setup !" << endl;
		return 1;
	}
	
	min_values = new double[6];
	max_values = new double[6];
	nbs_values = new unsigned int[6];
	_nb_t = nb_t;
	_kernel_size = kernel_size;
	
	_width = width;
	_height = height;
	min_values[5] = theta_min * M_PI / 180.0;
	max_values[5] = theta_max * M_PI / 180.0;
	min_values[3] = 0;
	max_values[3] = 0;

	
	nbs_values[0] = 1;
	nbs_values[1] = 1;
	nbs_values[2] = nb_a;
	nbs_values[3] = 1;
	nbs_values[4] = nb_c;
	nbs_values[5] = nb_theta;
	 
	alloc();
	
	return 0;
	
	
	
}

int c_eyelid_segmentation :: setup ( 	unsigned int width,
										unsigned int height,
										api_parameters & params,
										ostream * stream,
										const char * n_space,
										const char * nb_t_n,
										const char * nb_a_n,
										const char * nb_c_n,
										const char * nb_theta_n,
										const char * theta_min_n,
										const char * theta_max_n,
										const char * kernel_size_n)
{
	int q = 0;
	
	unsigned int nb_t, nb_a, nb_c, nb_theta, kernel_size;
	
	double theta_min, theta_max;
	
	stringstream oss;
	
	oss << n_space << "::" << nb_t_n;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_t,
									stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << kernel_size_n;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&kernel_size,
									stream	) )
		q = 1;
	oss.str("");
	
	
	
	oss << n_space << "::" << nb_a_n;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_a,
									stream	) )
		q = 1;
	oss.str("");
		
		
	oss << n_space << "::" << nb_c_n;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_c,
									stream	) )
		q = 1;
	oss.str("");
		
	oss << n_space << "::" << nb_theta_n;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&nb_theta,
									stream	) )
		q = 1;
	oss.str("");
		
		
	oss << n_space << "::" << theta_min_n;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&theta_min,
							stream	) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << theta_max_n;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&theta_max,
							stream	) )
		q = 1;
	oss.str("");
	
	if ( q )
		return 1;
		
	return ( setup(  	width,
						height,
						nb_t,
						nb_a,
						nb_c,
						nb_theta,
						theta_min,
						theta_max,
						kernel_size,
						stream ) );
}


int c_eyelid_segmentation :: compute( 	const IplImage * image,
										const IplImage * mask,
										const double & x_iris,
										const double & y_iris,
										const double & r_iris,
										const double & x_pupil,
										const double & y_pupil,
										const double & a_pupil,
										const double & b_pupil )
{
	
	memset( _mask->imageData, 255, _mask->widthStep * _mask->height );
	double params[6];
	//2 premières paraboles!
		min_values[0] = x_iris;
		max_values[0] = x_iris;
		
		min_values[1] = y_iris;
		max_values[1] = y_iris;
		
		
		min_values[2] = -1.0 / ( 2 * r_iris );
		max_values[2] = 0;
		
		
		min_values[4] = MIN( a_pupil, b_pupil ) - y_pupil + y_iris + _kernel_size;
		max_values[4] = r_iris;

		op->compute(	params,
						image,
						mask,
						min_values,
						max_values,
						nbs_values,
						6,
						-1.5 * r_iris,
						1.5 * r_iris,
						_nb_t,
						parabol2,
						d_parabol2,
						abs_more_than );
						
		_a_lower = params[2];
		_c_lower = params[4];
		_theta_lower = params[5];
		
		update_mask( params, r_iris, false );
	 
	
		min_values[0] = x_iris;
		max_values[0] = x_iris;
		
		min_values[1] = y_iris;
		max_values[1] = y_iris;
		
		
		max_values[2] = 1.0 / ( 2 * r_iris );
		min_values[2] = 0;
		
		
		max_values[4] = - MIN( a_pupil, b_pupil ) + y_pupil - y_iris - _kernel_size;
		min_values[4] = - r_iris;

		op->compute(	params,
						image,
						mask,
						min_values,
						max_values,
						nbs_values,
						6,
						-1.5 * r_iris,
						1.5 * r_iris,
						_nb_t,
						parabol2,
						d_parabol2,
						abs_more_than );
	
		_a_upper = params[2];
		_c_upper = params[4];
		_theta_upper = params[5];
	
		update_mask( params, r_iris, true );
	
		//Erosion du masque
		cvMorphologyEx( _mask,
						_mask,
						NULL,
						structuring_element_1,
						CV_MOP_ERODE,
						1);	
	
	return 0;
	
	
	
	
}

c_eyelid_segmentation :: ~c_eyelid_segmentation()
{
	free();
	initialize();
}



void c_eyelid_segmentation :: initialize()
{
	op = 0;
	_mask = 0;
	_mask_0 = 0;
	_width = 0;
	_height = 0;  	 
	err_stream = 0;
	min_values = 0;
	max_values = 0;
	nbs_values = 0;
	_nb_t = 0;
	_kernel_size = 0;
	structuring_element_1 = 0;
	_a_upper = 0;
	_c_upper = 0;
	_theta_upper = 0;
	
	_a_lower = 0;
	_c_lower = 0;
	_theta_lower = 0;
	
	
}

void c_eyelid_segmentation ::free()
{
	if ( op )
		delete op;
		
	if ( _mask )
		cvReleaseImage( &_mask );
		
	if ( _mask_0 )
		cvReleaseImage( &_mask_0 );
		
	if ( min_values )
		delete[] min_values;
	
	if ( max_values )
		delete[] max_values;
		
	if ( nbs_values )
		delete[] nbs_values;
	if ( structuring_element_1 )
		cvReleaseStructuringElement(&structuring_element_1);
}

void c_eyelid_segmentation :: alloc()
{

	op = new c_integro_differential_2d( _width, _height, 2 * _kernel_size + 1 );
	_mask = cvCreateImage( cvSize( _width, _height), IPL_DEPTH_8U, 1 );
	_mask_0 = cvCreateImage( cvSize( _width, _height), IPL_DEPTH_8U, 1 );
	memset( _mask_0->imageData, 255, 
			_mask_0->widthStep * _mask_0->height );
	structuring_element_1 = cvCreateStructuringElementEx( 2 * _kernel_size + 1,
														  2 * _kernel_size + 1,
														  _kernel_size,
														  _kernel_size,
														  CV_SHAPE_ELLIPSE,
														  NULL);
}

void c_eyelid_segmentation :: update_mask( const double * params, const double & r, bool up )
{
	double t_min = -r,
		   t_max = r,
		   dt = (t_max - t_min) / ( 50 - 1 );
	
	
	CvPoint pts[4];
	
	
	for ( unsigned int i = 1; i < 50; ++ i )
	{
		
		double 	t_0 = t_min + dt * (i - 1),
				t_1 = t_min + dt * i;
		
		double 	x_0, y_0,
				x_1, y_1;
	
		parabol2.get_value( x_0, y_0, t_0, params );
		parabol2.get_value( x_1, y_1, t_1, params );

			pts[0].x = x_0; 
			pts[0].y = y_0;

			pts[1].x = x_1;
			pts[1].y = y_1;

			pts[2].x = x_1; 
			pts[3].x = x_0; 
			

		if ( up )
		{ 
			pts[2].y = 0;
			pts[3].y = 0;
			
			pts[0].y += _kernel_size;
			pts[1].y += _kernel_size;	
			
		}
		else
		{
			pts[2].y = _height;
			pts[3].y = _height;
			
			pts[0].y -= _kernel_size;
			pts[1].y -= _kernel_size;	
		}

		cvFillConvexPoly( _mask, pts, 4, CV_RGB(0,0,0) );
	}
}
