#include "non_max_sup.hpp"

#include <iostream>
using namespace std;

int non_max_sup ( 	double * image_data,
					unsigned int image_width_step,
					const double * grad_data,
					unsigned int grad_width_step,
					const double * orient_data,
					unsigned int orient_width_step,
					unsigned int width,
					unsigned int height,
					double radius )
{
	double 	h_frac[ 181 ],
			v_frac[ 181 ],
			x_off [ 181 ],
			y_off [ 181 ];
	if ( radius < 1 )
		return 1;
	
	
	// Precalculate x and y offsets relative to centre pixel for each orientation angle 
	{
		double angle = 0;
		double step = M_PI / 180.0;
		for ( unsigned int i = 0; i <= 180; ++i )
		{
			//% Array of angles in 1 degree increments (but in radians).
			x_off[i] = radius * cos(angle);	//x and y offset of points at specified radius and angle
			y_off[i] = radius * sin(angle);	//from each reference position.
			h_frac[i] = x_off[i] - floor( x_off[i] ); //% Fractional offset of xoff relative to integer location
			v_frac[i] = y_off[i] - floor( y_off[i] ); //% Fractional offset of yoff relative to integer location
			angle += step;
		}
	}
	
	
	unsigned int i_radius = (unsigned int) ceil(radius);
	unsigned int row_max = (unsigned int) ceil(height - i_radius);
	unsigned int col_max = (unsigned int) ceil(width - i_radius);
	
	//0
	for ( unsigned int row = 0; row < height; ++ row)
	{
		memset ( image_data + row * image_width_step, 0, sizeof( double ) * width );
	}
	
	
	
	for ( unsigned int row = i_radius; row < row_max; ++ row)
	{

		for (unsigned col = i_radius; col < col_max; ++ col) 
		{
			double x,
				   y;

			unsigned int ori;
			
			double tl,
				   tr,
				   bl,
				   br;
				   
			int  fx,
				 cx,
				 fy,
				 cy;
			
			double 	upper_avg,
					lower_avg,
					v1;
			
			ori = (unsigned int) 180 * orient_data[row * orient_width_step + col] / M_PI + 0.5;   //% Index into precomputed arrays
			
			x = col + x_off[ori];     //% x, y location on one side of the point in question
			y = row - y_off[ori];

			fx = (int)floor(x);          //% Get integer pixel locations that surround location x,y
			cx = (int)ceil(x);
			fy = (int)floor(y);
			cy = (int)ceil(y);
			
			tl = grad_data[fy * grad_width_step + fx];   // % Value at top left integer pixel location.
			tr = grad_data[fy * grad_width_step + cx];  // % top right
			bl = grad_data[cy * grad_width_step + fx];   // % bottom left
			br = grad_data[cy * grad_width_step + cx];   // % bottom right

			upper_avg = tl + h_frac[ori] * (tr - tl);  //% Now use bilinear interpolation to
			lower_avg = bl + h_frac[ori] * (br - bl);  //% estimate value at x,y
			v1 = upper_avg + v_frac[ori] * (lower_avg - upper_avg);

			unsigned int n_pixel = row * grad_width_step + col;
			if ( grad_data[n_pixel] > v1) //% We need to check the value on the other side...
			{

				x = col - x_off[ori];    // % x, y location on the `other side' of the point in question
				y = row + y_off[ori];
	
				fx = (int) floor(x);
				cx = (int) ceil(x);
				fy = (int) floor(y);
				cy = (int) ceil(y);
				
				tl = grad_data[fy * grad_width_step + fx];   // % Value at top left integer pixel location.
				tr = grad_data[fy * grad_width_step + cx];  // % top right
				bl = grad_data[cy * grad_width_step + fx];   // % bottom left
				br = grad_data[cy * grad_width_step + cx];   // % bottom right

				upper_avg = tl + h_frac[ori] * (tr - tl);  //% Now use bilinear interpolation to
				lower_avg = bl + h_frac[ori] * (br - bl);  //% estimate value at x,y
				double v2 = upper_avg + v_frac[ori] * (lower_avg - upper_avg);

				if ( grad_data[n_pixel] > v2) //% We need to check the value on the other side...
				{
					image_data[row * image_width_step + col] = grad_data[n_pixel];
				}
			}
		}
	}
	return 0;
}


int non_max_sup ( 	IplImage * image_out,
					const IplImage * image_grad,
					const IplImage * image_orient,
					double radius )
{
	double * image_data;
	unsigned int image_width_step;
	const double * grad_data;
	unsigned int grad_width_step;
	const double * orient_data;
	unsigned int orient_width_step;
	
	
	
	unsigned int img_width;
	unsigned int img_height;
	
	unsigned int grad_width;
	unsigned int grad_height;
	
	
	unsigned int orient_width;
	unsigned int orient_height;
	
	
	//Check-up
	if ( image_out->depth != IPL_DEPTH_64F || image_grad->depth != IPL_DEPTH_64F || image_orient->depth != IPL_DEPTH_64F )
		return -1;
	
	//Dim.
	{
		image_width_step = image_out->widthStep / sizeof(double);
		if ( image_out->roi != NULL )
		{

			image_data = ( (double*) image_out->imageData ) + image_out->roi->xOffset + image_out->roi->yOffset * image_width_step;
			img_width = image_out->roi->width;
			img_height = image_out->roi->height;
		}
		else
		{
			image_data = ( (double*) image_out->imageData );
			img_width = image_out->width;
			img_height = image_out->height;
		}
	}
	
	{
		grad_width_step = image_grad->widthStep / sizeof(double);
		if ( image_grad->roi != NULL )
		{

			grad_data = ( (const double*) image_grad->imageData ) + image_grad->roi->xOffset + image_grad->roi->yOffset * grad_width_step;
			grad_width = image_grad->roi->width;
			grad_height = image_grad->roi->height;
		}
		else
		{
			grad_data  = ( (const double*) image_grad->imageData );
			grad_width = image_grad->width;
			grad_height = image_grad->height;
		}
	}
	
	{
		orient_width_step = image_orient->widthStep / sizeof(double);
		if ( image_orient->roi != NULL )
		{

			orient_data = ( (const double*) image_orient->imageData ) + image_orient->roi->xOffset + image_orient->roi->yOffset * orient_width_step;
			orient_width = image_orient->roi->width;
			orient_height = image_orient->roi->height;
		}
		else
		{
			orient_data = ( (const double*) image_orient->imageData );
			orient_width = image_orient->width;
			orient_height = image_orient->height;
		}
	}
	
	if ( img_width != grad_width || img_width != orient_width || img_height != grad_height || img_height != orient_height )
	{	
		cerr << img_width << "\t" << grad_width << "\t" << orient_width << endl;
		cerr << img_height << "\t" << grad_height << "\t" << orient_height << endl;
		return -1;
	}
		
	return non_max_sup ( 	image_data,
							image_width_step,
							grad_data,
							grad_width_step,
							orient_data,
							orient_width_step,
							img_width,
							img_height,
							radius );
	
}
