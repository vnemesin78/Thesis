#include "c_hough_ellipse.hpp"
#include <opencv/highgui.h>

c_hough_ellipse :: c_hough_ellipse()
{
	initialize();
}

c_hough_ellipse :: c_hough_ellipse(	unsigned int nb_x,
										unsigned int nb_y,
										unsigned int nb_thetas,
										unsigned int r_min,
										unsigned int r_max,
										unsigned int thickness )
{
	initialize();
	setup( 	nb_x,
			nb_y,
			nb_thetas,
			r_min,
			r_max,
			thickness );
}

void c_hough_ellipse :: setup ( void )
{
	free();
	initialize();
}

int c_hough_ellipse :: setup (	unsigned int nb_x,
									unsigned int nb_y,
									unsigned int nb_thetas,
									unsigned int r_min,
									unsigned int r_max,
									unsigned int thickness )
{
	free();
	initialize();
	
	if ( 	!nb_x			||
			!nb_y			||
			!nb_thetas		||
			!thickness		||
			r_min >= r_max	)
	{
		*err_stream << "Error in int c_hough_ellipse :: setup : Invalid argument(s)!" << endl;
		return 1;
	}
	
	_nb_x = nb_x;
	_nb_y = nb_y;
	_nb_thetas = nb_thetas;
	_r_min = r_min;
	_r_max = r_max;
	_nb_r = _r_max - _r_min;
	_thickness = thickness;
	//Super alloc mémoire!
	_hough_space = new unsigned int [ _nb_x * _nb_y * _nb_r * _nb_r * _nb_thetas];
	memset ( _hough_space, 0, sizeof(int) * _nb_x * _nb_y * _nb_r * _nb_r * _nb_thetas );
	//Génération des ellipses
	generate_templates();
	return 0;
}

int c_hough_ellipse :: compute_hough_transform( 	const unsigned char * img_data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													int x_start,
													int x_end,
													int y_start,
													int y_end,											
													unsigned int r_min,
													unsigned int r_max )
{
	unsigned int	nb_x = x_end - x_start,
					nb_y = y_end - y_start;
	
	if ( 	nb_x == 0 		||
			nb_y == 0		||
			img_data == 0 	||
			r_min >= r_max	||
			width == 0		||
			height == 0		)
	{
		*err_stream << "Error in int c_hough_ellipse :: compute_hough_transform : Invalid argument(s) " << endl;
		return 1;
		
	}
	
	if ( width_step == 0 )
		width_step = width;
	
	
	//Check dim.	
	if ( 	nb_x > _nb_x 		|| 
			nb_y > _nb_y 		|| 
			r_min < _r_min 	||
			r_max > _r_max		)
	{
		*err_stream << "Warning : c_hough_ellipse object re-allocation" << endl;
		setup( 	nb_x, 
				nb_y,
				_nb_thetas,
				r_min,
				r_max );
	}
	
	// 0
	memset ( 	_hough_space, 
				0, 
				sizeof(int) * _nb_x * _nb_y * _nb_r * _nb_r * _nb_thetas );
	
	//Computation
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( img_data[ i * width_step + j ] )
			{
				add_ellipses ( i, 
							   j,
							   x_start,
							   x_end,
							   y_start,
							   y_end, 
							   r_min, 
							   r_max );
			}
		}
	}	
	return 0;
}

int c_hough_ellipse :: compute_hough_transform( 	const unsigned int * contour,
													unsigned int nb_pts_contour,
													int x_start,
													int x_end,
													int y_start,
													int y_end,											
													unsigned int r_min,
													unsigned int r_max )
{
	unsigned int	nb_x = x_end - x_start,
					nb_y = y_end - y_start;
	
	if ( 	nb_x == 0 			||
			nb_y == 0			||
			contour == 0 		||
			! nb_pts_contour	||
			r_min >= r_max		)
	{
		*err_stream << "Error in int c_hough_ellipse :: compute_hough_transform : Invalid argument(s) " << endl;
		return 1;
		
	}
	
	//Check dim.	
	if ( 	nb_x > _nb_x 		|| 	
			nb_y > _nb_y 		|| 
			r_min < _r_min 	||
			r_max > _r_max		)
	{
		*err_stream << "Warning : c_hough_ellipse object re-allocation" << endl;
		setup( 	nb_x, 
				nb_y,
				_nb_thetas,
				r_min,
				r_max );
	}
	// 0
	memset ( 	_hough_space, 
				0, 
				sizeof(int) * _nb_x * _nb_y * _nb_r * _nb_r * _nb_thetas );
	
	for ( unsigned int i = 0; i < nb_pts_contour; ++ i )
	{
		add_ellipses ( contour[2 * i + 1], 
					   contour[2 * i],
					   x_start,
					   x_end,
					   y_start,
					   y_end, 
					   r_min, 
					   r_max );		
	}
	
	
	
	
	
	
	return 0;
}

unsigned int c_hough_ellipse :: search_best_ellipse( 	int & x,
															int & y,
															unsigned int & a,
															unsigned int & b,
															double & theta,
															int x_start,
															int x_end,
															int y_start,
															int y_end,											
															unsigned int r_min,
															unsigned int r_max )
{
	unsigned int	dx = x_end - x_start,
					dy = y_end - y_start,
					dr = r_max - r_min;
	x = 0;
	y = 0;
	a = 0;
	b = 0;
	theta = - M_PI / 2;
	unsigned int max = 0;
	for ( unsigned int i = r_min; i < r_max; ++ i )
	{
		for ( unsigned int j = r_min; j < r_max; ++ j )
		{
			for ( unsigned int k = 0; k < _nb_thetas; ++ k )
			{
				for ( unsigned int l = 0; l < dy; ++ l )
				{
					for ( unsigned int m = 0; m < dx; ++ m )
					{
						unsigned int v = _hough_space[ ( ( ( ( i - r_min) * dr + (j - r_min) ) * _nb_thetas + k ) * dy + l ) * dx + m ];
						if ( v >= max )
						{
							max = v;
							x = m + x_start;
							y = l + y_start;
							a = i;
							b = j;
							theta = ( k ) * M_PI / ( 2 * _nb_thetas );
						}				
					}					
				} 
			}
		}
	}
	return max;
}







void c_hough_ellipse :: add_ellipses ( 	unsigned int y,
											unsigned int x,
											int x_start,
											int x_end,
											int y_start,
											int y_end,											
											unsigned int r_min,
											unsigned int r_max )
{
	int dx = x_end - x_start,
		dy = y_end - y_start,
		dr = r_max - r_min;
	for ( unsigned int i = r_min; i < r_max; ++ i )
	{
		for ( unsigned int j = r_min; j < r_max; ++ j )
		{
			for ( unsigned int k = 0; k < _nb_thetas; ++ k )
			{
				unsigned int 	n = ( (i - _r_min ) * _nb_r + (j - _r_min ) ) * _nb_thetas + k,
								p = ( (i - r_min ) * dr + (j - r_min ) ) * _nb_thetas + k;
				
				for ( unsigned int l = 0; l < _nb_ellipses_points[n]; ++ l )
				{
					int 	_x = _template_ellipses[n][2 * l] + x,
							_y = _template_ellipses[n][2 * l + 1] + y;					
					
					if ( _x >= x_start && _x < x_end &&
						 _y >= y_start && _y < y_end )
					{
						_hough_space[ p * dx * dy + dx * ( _y - y_start ) + _x - x_start] ++; 
					}
				}
			}
		}
	}
}

void c_hough_ellipse :: generate_templates ( )
{
	

	unsigned int dim;
	unsigned char * temp_data;
	unsigned int width_step;
	
	IplImage * temp;					
	dim = _r_max * 2 + 1;
	temp =  cvCreateImage (	cvSize(dim, dim),
							8,
							1 );
	temp_data = (unsigned char*) temp->imageData;
	width_step = temp->widthStep;
	//Alloc mémoire
	_template_ellipses = new int *[ _nb_thetas * _nb_r * _nb_r ];
	_nb_ellipses_points = new unsigned int[ _nb_r * _nb_r * _nb_thetas ];
	
	for ( unsigned int i = _r_min; i < _r_max; ++ i )
	{
		for ( unsigned int j = _r_min; j < _r_max; ++ j )
		{
			for ( unsigned int k = 0; k < _nb_thetas; ++ k ) 
			{
				unsigned int n = ( ( i - _r_min ) * _nb_r + ( j - _r_min ) ) * _nb_thetas + k; 
				//Mise à zéro
				memset( temp_data,
						0,
						sizeof(char) * temp->widthStep * temp->height );
		
				cvEllipse( 	temp,
							cvPoint( _r_max, _r_max ),
							cvSize ( i, j ),
							(90.0 * k) / _nb_thetas,
							0,
							360,
							CV_RGB(255, 255, 255 ),
							_thickness ); 
				_nb_ellipses_points[ n ] = 0;
				for ( unsigned int y = 0; y < dim; ++ y )
				{
					for ( unsigned int x = 0; x < dim; ++ x )
					{
						if ( temp_data[ y * width_step + x ] )
							( _nb_ellipses_points[n] ) ++;
					}					
				}
			
				unsigned int nb_pts = 0;
				_template_ellipses[ n ] = new int[ 2 * _nb_ellipses_points[ n ] ];
				for ( unsigned int y = 0; y < dim; ++ y )
				{
					for ( unsigned int x = 0; x < dim; ++ x )
					{
						if ( temp_data[ y * width_step + x ] )
						{
							_template_ellipses[ n ][ 2 * nb_pts ] = x - _r_max;
							_template_ellipses[ n ][ 2 * nb_pts + 1 ] = y - _r_max;
							nb_pts ++;
						}
					}
				}
			}
		}		
	}
	cvReleaseImage( &temp );
}

void c_hough_ellipse :: free ( )
{
	if ( _hough_space )
		delete[] _hough_space;
		
	if ( _template_ellipses )
	{
		unsigned int nb = _nb_r * _nb_r * _nb_thetas;
		for ( unsigned int i = 0; i < nb; ++ i )
		{
			delete[] _template_ellipses[i];
		}
		delete[] _template_ellipses;
		
	}
	if ( _nb_ellipses_points )
		delete[] _nb_ellipses_points;
}

void c_hough_ellipse :: initialize()
{
	_hough_space = 0;
	_nb_x = 0;
	_nb_y = 0;
	_nb_r = 0;
	_nb_thetas = 0;
	_r_min = 0;
	_r_max = 0;
	err_stream = &cout;
	_template_ellipses = 0; //  4 PI ( r_max )³ x nb_thetas
	_nb_ellipses_points = 0;
}

c_hough_ellipse :: ~c_hough_ellipse()
{
	free();
	initialize();
}
