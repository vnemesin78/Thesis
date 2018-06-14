#include "c_hough.hpp"
#include <opencv/highgui.h>
#include <sstream>
#include <iostream>
#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
#include <stdexcept>
#include <cstring>
#include <cmath>
c_hough :: c_hough()
{
	initialize();
}

c_hough :: c_hough(	unsigned int width,
						unsigned int height,
						unsigned int r_max,
						unsigned int thickness ) // hough_space, width, height
{
	initialize();
	if ( setup ( 	width, 
					height, 
					r_max, 
					thickness ) )
	{
		throw ( invalid_argument("Arguments of c_hough !") );
	}
}

/**@fn
 * @param width : largeur de l'image
 * @param height : hauteur de l'image
 * @param r_max : rayon max.
 * @brief
 * Setup
 * 
 **/
int c_hough :: setup (	unsigned int width,
						unsigned int height,
						unsigned int r_max,
						unsigned int thickness ) // hough_space, width, height, r_step, _width_step
{
	free();
	initialize();
	
	if ( ! width || ! height || ! r_max || !thickness) 
	{
		*err_stream << "Error : bad arguments in int c_hough :: setup ( unsigned int width, unsigned int height, unsigned int r_max ); " << endl;
		return 1;
	}
	
	_width = width;
	_height = height;
	_r_max = r_max;
	_thickness = thickness;
	_width_step = _width + 2 * _r_max;
	_r_step =  ( _height + 2 * _r_max ) * _width_step;
	//Espace de hough
	_hough_space = new unsigned int[ _r_step * _r_max ];

	//Circle templates
	generate_templates();
	
	return 0;
}

/**@fn
 * @brief
 * Generate the circle images ( for r = 0 to r_max )
 * 
 * 
 **/
void c_hough :: generate_templates ( )  // r_max, templates, ROIs
{
	unsigned int dim;
	CvRect ROI;
	
	IplImage * temp;					
	dim = _r_max * 2 + 1;
	temp =  cvCreateImage (	cvSize(dim, dim),
							8,
							1 );
	_circles_points = new int * [_r_max];
	_nb_circles_points = new unsigned int[_r_max];
	
	//Circles
	for ( unsigned int i = 1; i <= _r_max; ++ i )
	{
		unsigned int j = i - 1,
					 width_step;
		unsigned char * temp_data;

		//Temp. data
		width_step = temp->widthStep;
		temp_data = (unsigned char *) temp->imageData;

		//0
		for (unsigned int x = 0; x < dim; ++ x )
		{
			memset( temp_data + x * width_step, 0, sizeof(char) * dim );
		}
		
		cvCircle ( 	temp,
					cvPoint( _r_max, _r_max ),
					i,
					CV_RGB(255, 255, 255 ),
					_thickness ); 

		//ROI computation
		ROI.x = _r_max - i;
		ROI.y = _r_max - i;
		ROI.width = 2 * i + 1;
		ROI.height = 2 * i + 1;
		
		//Nb points
		_nb_circles_points[j] = 0;
		for ( int y = 0; y < ROI.height; ++ y )
		{
			for ( int x = 0; x < ROI.width; ++ x )
			{
				if ( temp_data[ ( ROI.y  + y ) * width_step + x + ROI.x ] )
					( _nb_circles_points[j] ) ++;
			}
			
		}
		//Alloc
		_circles_points[j] = new int[ 2 * _nb_circles_points[j] ];
		
		unsigned int nb_pts = 0;
		for ( int y = 0; y < ROI.height; ++ y )
		{
			for ( int x = 0; x < ROI.width; ++ x )
			{
				if ( temp_data[ ( ROI.y  + y ) * width_step + x + ROI.x ] )
				{
					_circles_points[j][ 2 * nb_pts ] = x - _r_max + ROI.x;
					_circles_points[j][ 2 * nb_pts + 1 ] = y - _r_max + ROI.y;
					nb_pts ++;
				}
			}
		}
	
	
	}
	
	cvReleaseImage( &temp );
}

/**@fn 
 * @param img_data : image pixels
 * @param width : image width ( <= _width)
 * @param height : image height ( <= _height )
 * @param width_step : image widthStep
 * @param r_min : min radius ( >= 1 ) 
 * @param r_max : max radius ( <= _r_max )
 * @brief
 * Compute the hough transform.
 **/
int c_hough :: compute_hough_transform( 	const unsigned char * img_data,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											unsigned int r_min,
											unsigned int r_max )
{
	if ( ! img_data || ! height || ! width || ! r_min)
	{
		*err_stream << "Error : bad params in int  c_hough :: compute_hough_transform( const unsigned char * img_data, unsigned int width, unsigned int height, unsigned int width_step , unsigned int r_min, unsigned int r_max) " << endl;
		return -3;
	}
	
	if ( width_step == 0 )
		width_step = width;
	
	if ( r_max == 0 )
		r_max = _r_max;
	
	if ( r_max > _r_max || width > _width || height > _height )
	{
		*err_stream << "Error : too large width, height or r_max in int  c_hough :: compute_hough_transform( const unsigned char * img_data, unsigned int width, unsigned int height, unsigned int width_step , unsigned int r_min, unsigned int r_max " << endl;
		return -1;
	}
	
	// 0
	memset ( _hough_space, 0, sizeof(int) * _r_step * _r_max );
	
	//Computation
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( img_data[ i * width_step + j ] )
			{
				add_circles ( i, j, r_min, r_max );
			}
		}
	}
	return 0;
}

int c_hough :: compute_hough_transform( 	const unsigned int * contour,
											unsigned int nb_pts_contour,
											unsigned int width,
											unsigned int height,
											unsigned int r_min,
											unsigned int r_max )
{
	if ( 	! contour 	|| 
			! height 	|| 
			! width 	|| 
			! r_min     )
	{
		*err_stream << "Error : bad params in int  c_hough :: compute_hough_transform( const unsigned char * img_data, unsigned int width, unsigned int height, unsigned int width_step , unsigned int r_min, unsigned int r_max) " << endl;
		return -3;
	}
	
	if ( r_max == 0 )
		r_max = _r_max;
	
	if ( 	r_max > _r_max 	|| 
			width > _width 	|| 
			height > _height )
	{
		*err_stream << "Error : too large width, height or r_max in int  c_hough :: compute_hough_transform( const unsigned char * img_data, unsigned int width, unsigned int height, unsigned int width_step , unsigned int r_min, unsigned int r_max " << endl;
		return -1;
	}
	
	// 0
	memset ( _hough_space, 0, sizeof(int) * _r_step * _r_max );
	for ( unsigned int i = 0; i < nb_pts_contour; ++i) 
	{
		add_circles( 	contour[2 * i + 1],
						contour[2 * i],
						r_min,
						r_max );
	}
	
	
	
	return 0;
}


/**@fn
 * @param x, y : coordinate of pixel
 * @param r_min : min_radius.
 * @param r_max : max_radius.
 * @brief
 * Add a pixel on the hough transform.
 **/
void c_hough :: add_circles ( 	unsigned int y, 
								unsigned int x, 
								unsigned int r_min, 
								unsigned int r_max )
{
	
	
	for ( unsigned int j = r_min - 1; j< r_max; ++ j )
	{
		for ( unsigned int i = 0; i < _nb_circles_points[j]; ++ i )
		{

			_hough_space[ _r_step * j + _width_step * ( _r_max + y + _circles_points[j][ 2 * i + 1]) +  ( _r_max + x + _circles_points[j][ 2 * i]) ] ++;
		}
	}
}

void c_hough :: free ()
{
	if ( _hough_space )
		delete[] _hough_space;
		
	if ( _circles_points )
	{
		for ( unsigned int j = 0; j < _r_max; ++ j )
		{
			delete[] _circles_points[j];
		}
		delete[] _circles_points;
	}
	if ( _nb_circles_points )
		delete[] _nb_circles_points;
	
}


void c_hough :: initialize()
{
	err_stream = &cout;
	_hough_space = 0;
	_width = 0;
	_height = 0;
	_r_max = 0;
	_r_step = 0;
	_width_step = 0;
	_nb_circles_points = 0;
	_circles_points = 0;
	
	
	
}
c_hough :: ~c_hough()
{
	free();
	initialize();
}

unsigned int c_hough :: get_best_circle (	int & x,
											int & y,
											unsigned int & r,
											unsigned int r_min,
											unsigned int r_max,
											unsigned int width,
											unsigned int height )
{
	x = 0;
	y = 0;
	r = r_min;
	unsigned int max = _hough_space[ _r_step * ( r - 1 ) + _width_step * ( _r_max + y ) +  ( _r_max + x ) ];
	for ( unsigned int i = r_min; i <= r_max; ++ i )
	{
		for ( int j = -i; j < (long long) height + i; ++ j )
		{
			for ( int k = -i; k < (long long) width + i; ++ k )
			{
				unsigned int v =  _hough_space[ _r_step * ( i - 1 ) + _width_step * ( _r_max + j ) +  ( _r_max + k ) ] ;
				if ( v > max )
				{
					y = j;
					x = k;
					r = i;
					max = v;
				}
			}
		}
	}
	return max;
}
