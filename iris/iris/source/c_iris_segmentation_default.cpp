#include "c_iris_segmentation.hpp"
#include "iris_default.hpp"
int c_iris_segmentation :: default_setup( 	unsigned int width,
												unsigned int height,
												ostream * _err_stream )
{
	unsigned int nb_directions;
	IRIS__NB_DIRECTIONS( nb_directions);
	unsigned int kernel_size;
	IRIS__KERNEL_SIZE( kernel_size );
	
	unsigned int nb_samples;
	IRIS__NB_SAMPLES( nb_samples );
	
	unsigned int dx;
	IRIS__DX( dx );
	
	unsigned int nb_x;
	IRIS__NB_X( nb_x );
	
	unsigned int dy;
	IRIS__DY( dy );
	
	unsigned int nb_y;
	IRIS__NB_Y( nb_y );
	
	unsigned int r_min;
	IRIS__R_MIN( r_min );
	
	unsigned int r_max;
	IRIS__R_MAX( r_max );
	
	double r_ratio;
	IRIS__R_RATIO( r_ratio );
	
	unsigned int nb_r_ellipse;
	IRIS__NB_SAMPLES_ELLIPSE( nb_r_ellipse );
	
	unsigned int dx_ellipse;
	IRIS__DX_ELLIPSE(dx_ellipse);
	
	unsigned int nb_x_ellipse;
	IRIS__NB_X_ELLIPSE(nb_x_ellipse);
	
	unsigned int dy_ellipse;
	IRIS__DY_ELLIPSE(dy_ellipse);
	
	unsigned int nb_y_ellipse;
	IRIS__NB_Y_ELLIPSE(nb_y_ellipse);
	
	unsigned int nb_thetas;
	IRIS__NB_THETAS(nb_thetas);
	
	return ( setup (	width,
					height,
					nb_directions,
					kernel_size,
					nb_samples,
					dx,
					nb_x,
					dy,
					nb_y,
					r_min,
					r_max,
					r_ratio,
					nb_r_ellipse,
					dx_ellipse,
					nb_x_ellipse,
					dy_ellipse,
					nb_y_ellipse,
					nb_thetas,
					_err_stream	) );
}
