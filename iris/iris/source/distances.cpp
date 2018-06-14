#include "distances.hpp"

template <class type> int distance :: debug_convert( const type * code,
													 const type * mask,
													 const type * fb,
													 unsigned int width_step,
													 unsigned int width,
													 unsigned int height )
{
	unsigned int end = width_step - 1;
	IplImage * image = cvCreateImage( cvSize(width, height), IPL_DEPTH_8U, 3);
	if ( width % (8*sizeof(type)) == 0 )
		end = width_step;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		unsigned int n_pixel, n_pixel_img;
		for ( unsigned int j = 0; j < end; ++ j )
		{
			n_pixel = i * width_step + j;
			for( unsigned int k = 0; k < 8 * sizeof(type); ++ k )
			{
				n_pixel_img = i * image->widthStep + j * 3 * 8 * sizeof(type) + 3 * (8 * sizeof(type) - 1 - k);
				if ( (mask[n_pixel] >> k & 1 ) == 0 )
				{
					
						((unsigned char*) image->imageData)[n_pixel_img] = 0;
						((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
						((unsigned char*) image->imageData)[n_pixel_img + 2] = 255;	

				}
				else if ( ( fb[n_pixel] >> k & 1 ) == 0 )
				{
					((unsigned char*) image->imageData)[n_pixel_img] = 255;
					((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
					((unsigned char*) image->imageData)[n_pixel_img + 2] = 0;	
				}
				
				else
				{
					if ( (code[n_pixel] >> k & 1 ) == 0 )
					{
						((unsigned char*) image->imageData)[n_pixel_img] = 0;
						((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
						((unsigned char*) image->imageData)[n_pixel_img + 2] = 0;	
					}
					else
					{
						((unsigned char*) image->imageData)[n_pixel_img] = 255;
						((unsigned char*) image->imageData)[n_pixel_img + 1] = 255;
						((unsigned char*) image->imageData)[n_pixel_img + 2] = 255;
					}
					
				}
				
				
			}
			
		}
		n_pixel = i * width_step + width_step - 1;
		unsigned int r = width - (width_step-1) * 8 * sizeof(type);
		unsigned int r2d2 = 8 * sizeof(type) - r;
		for( unsigned int k = 0; k < width % (8*sizeof(type)); ++ k )
		{
			n_pixel_img = i * image->widthStep + (width_step - 1) * 3 * 8 * sizeof(type) + 3 * (width % (8*sizeof(type)) - 1 - k);
			if ( (mask[n_pixel] >> (k + r2d2)  & 1 ) == 0 )
			{
				
					((unsigned char*) image->imageData)[n_pixel_img] = 0;
					((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
					((unsigned char*) image->imageData)[n_pixel_img + 2] = 255;	

			}
			else if ( ( fb[n_pixel] >> (k + r2d2) & 1 ) == 0 )
			{
				((unsigned char*) image->imageData)[n_pixel_img] = 255;
				((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
				((unsigned char*) image->imageData)[n_pixel_img + 2] = 0;	
			}
			else
			{
				if ( (code[n_pixel] >> (k + r2d2) & 1 ) == 0 )
				{
					((unsigned char*) image->imageData)[n_pixel_img] = 0;
					((unsigned char*) image->imageData)[n_pixel_img + 1] = 0;
					((unsigned char*) image->imageData)[n_pixel_img + 2] = 0;	
				}
				else
				{
					((unsigned char*) image->imageData)[n_pixel_img] = 255;
					((unsigned char*) image->imageData)[n_pixel_img + 1] = 255;
					((unsigned char*) image->imageData)[n_pixel_img + 2] = 255;
				}
				
			}
			
			
		}
		
		
		
	}
	
	cvShowImage("Reconstruction", image);
	cvWaitKey(0);
	
	cvReleaseImage(&image);
	
	
	
}



template <class type> int distance :: convert( 	type * code_out,
												type * mask_out,
												type * fb_data_out,
												unsigned int & width_step,
												const unsigned char * code_data,
												const unsigned char * fragility_map_data,
												unsigned int width,
												unsigned int height,
												unsigned int code_width_step,
												unsigned int fragility_map_width_step,
												unsigned int fm_threshold )
{
	unsigned int 	width_step_m1, 
					r,
					r2d2;
	width_step = ( width + (8 * sizeof(type) - 1 ) ) / ( 8 * sizeof(type) );
	width_step_m1 = width_step - 1;
	r = width - width_step_m1 * 8 * sizeof(type);
	r2d2 = 8 * sizeof(type) - r;
	//0
	memset( code_out, 
			0, 
			sizeof(type) * width_step * height );
	memset( mask_out, 
			0, 
			sizeof(type) * width_step * height );
	memset( fb_data_out, 
			0, 
			sizeof(type) * width_step * height );
			
			
	for ( unsigned int i = 0; i < height; ++ i )
	{
		unsigned int n_pixel;
		for ( unsigned int j = 0; j < width_step_m1; ++ j )
		{
			n_pixel = i * width_step + j;
			for( unsigned int k = 0; k < 8 * sizeof(type); ++ k )
			{

				
				//Décalage des bits
				code_out[ n_pixel ] = ( code_out[ n_pixel ] << 1 );
				mask_out[ n_pixel ] = ( mask_out[ n_pixel ] << 1 );
				fb_data_out[ n_pixel ] = ( fb_data_out[ n_pixel ] << 1 );
				
				//Conversion des données
				if ( code_data[ i * code_width_step + j * 8 * sizeof(type) + k] )
					 ++ code_out[ n_pixel ];
				if ( fragility_map_data[ i * code_width_step + j * 8 * sizeof(type) + k] > 0 )
					 ++ mask_out[ n_pixel ];
				if ( fragility_map_data[ i * code_width_step + j * 8 * sizeof(type) + k] >= fm_threshold )
					 ++ fb_data_out[ n_pixel ]; 
					 
			}
		}
		n_pixel = i * width_step + width_step_m1;
		for( unsigned int k = 0; k < r; ++ k )
		{
			//Décalage des bits
			code_out[ n_pixel ] = ( code_out[ n_pixel ] << 1 );
			mask_out[ n_pixel ] = ( mask_out[ n_pixel ] << 1 );
			fb_data_out[ n_pixel ] = ( fb_data_out[ n_pixel ] << 1 );
			//Conversion des données
			if ( code_data[ i * code_width_step + width_step_m1 * 8 * sizeof(type) + k] )
				++ code_out[ n_pixel ];
			if ( fragility_map_data[ i * code_width_step + width_step_m1 * 8 * sizeof(type) + k] > 0 )
				++ mask_out[ n_pixel ];
			if ( fragility_map_data[ i * code_width_step + width_step_m1 * 8 * sizeof(type) + k] >= fm_threshold )
				++ fb_data_out[ n_pixel ];
				
		}
		
		code_out[ n_pixel ] = ( code_out[ n_pixel ] << r2d2 );
		mask_out[ n_pixel ] = ( mask_out[ n_pixel ] << r2d2 );
		fb_data_out[ n_pixel ] = ( fb_data_out[ n_pixel ] << r2d2 );
			
	}
	
	
}

											


template <class type> int distance :: rotate( 	type * code_out,
												type * mask_out,
												type * fb_data_out,
												const type * code_in,
												const type * mask_in,
												const type * fb_data_in,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												int theta )
{
	unsigned int 	width_step_m1, 
					r,
					r2d2;
	width_step_m1 = width_step - 1;
	r = width - width_step_m1 * 8 * sizeof(type);
	r2d2 = 8 * sizeof(type) - r;
	memset( code_out, 
			0, 
			sizeof(type) * width_step * height );
	memset( mask_out, 
			0, 
			sizeof(type) * width_step * height );
	memset( fb_data_out, 
			0, 
			sizeof(type) * width_step * height );	

	int right_delta,
		left_delta;
		
	if ( theta < 0 )
	{
		right_delta = theta + width;
		left_delta = theta;
	}
	else
	{
		right_delta = theta;
		left_delta = theta - width;
	}
	
	//Décalage à droite
	{
		int d_case = right_delta / (sizeof(type) * 8 );
		int r = right_delta % (sizeof(type) * 8 ) ;

		int r_inv = sizeof(type) * 8  - r;

		for ( unsigned int i = 0; i < height; ++ i )
		{
			unsigned int n_pixel_out,
						 n_pixel_in;
			for ( unsigned int j = width_step - d_case - 1; j != (unsigned int) (0); -- j )
			{
				n_pixel_out = i * width_step + j + d_case;
				n_pixel_in = i * width_step + j;
				code_out[ n_pixel_out ] += code_in[ n_pixel_in ] >> r;
				if ( r_inv != sizeof(type) * 8 )
					code_out[ n_pixel_out ] += code_in[ n_pixel_in - 1 ] << r_inv;
				
				mask_out[ n_pixel_out ] += mask_in[ n_pixel_in ] >> r;
				if ( r_inv != sizeof(type) * 8 )
					mask_out[ n_pixel_out ] += mask_in[ n_pixel_in - 1 ] << r_inv;
				
				fb_data_out[ n_pixel_out ] += fb_data_in[ n_pixel_in ] >> r;
				if ( r_inv != sizeof(type) * 8 )
					fb_data_out[ n_pixel_out ] += fb_data_in[ n_pixel_in - 1 ] << r_inv;
				
				
				
			}
			n_pixel_out = i * width_step + d_case;
			n_pixel_in = i * width_step;
			code_out[ n_pixel_out ] += code_in[ n_pixel_in ] >> r;
			mask_out[ n_pixel_out ] += mask_in[ n_pixel_in ] >> r;
			fb_data_out[ n_pixel_out ] += fb_data_in[ n_pixel_in ] >> r;
			
			
		}
	}

	//Décalage à gauche
	{
		int d_case = left_delta / (sizeof(type) * 8 ) + 1;
		int r = left_delta % (sizeof(type) * 8 ) ;
		int r_inv = sizeof(type) * 8  - r;
		for ( unsigned int i = 0; i < height; ++ i )
		{
			unsigned int n_pixel_out,
						 n_pixel_in;
			for ( unsigned int j = - d_case; j < width_step_m1; ++ j )
			{
				n_pixel_out = i * width_step + j + d_case;
				n_pixel_in = i * width_step + j;
				if ( r_inv != sizeof(type) * 8 )
					code_out[ n_pixel_out ] += code_in[ n_pixel_in ] << r_inv;
				code_out[ n_pixel_out ] += code_in[ n_pixel_in + 1 ] >> r;
				if ( r_inv != sizeof(type) * 8 )
					mask_out[ n_pixel_out ] += mask_in[ n_pixel_in ] << r_inv;
				mask_out[ n_pixel_out ] += mask_in[ n_pixel_in + 1 ] >> r;
				
				if ( r_inv != sizeof(type) * 8 )
					fb_data_out[ n_pixel_out ] += fb_data_in[ n_pixel_in ] << r_inv;
				fb_data_out[ n_pixel_out ] += fb_data_in[ n_pixel_in + 1] >> r;
				
			}
			n_pixel_out = i * width_step + width_step_m1 + d_case;
			n_pixel_in = i * width_step + width_step_m1;
			if ( r_inv != sizeof(type) * 8 )
				code_out[ n_pixel_out ] += code_in[ n_pixel_in ] << r_inv;
			if ( r_inv != sizeof(type) * 8 )
				mask_out[ n_pixel_out ] += mask_in[ n_pixel_in ] << r_inv;
			if ( r_inv != sizeof(type) * 8 )
				fb_data_out[  n_pixel_out ] += fb_data_in[ n_pixel_in ] << r_inv;
		}
	}
	//Derniere collonne
	{

		type j = (type) -1;
		j = j << r2d2;
		for ( unsigned int i = 0; i < height; ++ i )
		{
			unsigned int n_pixel_out = i * width_step + width_step_m1;
			code_out[ n_pixel_out ] = ( code_out[ n_pixel_out ] ) & j;
			mask_out[ n_pixel_out ] = (mask_out[ n_pixel_out ]) & j;
			fb_data_out[ n_pixel_out ] = fb_data_out[ n_pixel_out ] & j;	
				

		}
	}
}

template <class type> int distance :: Hamming_opt ( double & d,
													const type * code_1_data,
													const type * mask_1_data,
													const type * fb_data_1,
													const type * code_2_data,
													const type * mask_2_data,
													const type * fb_data_2,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													int theta,
													const void * params,
													type * _code_2_data,
													type * _mask_2_data,
													type * _fb_data_2 )
{


	d = 0;

	rotate( _code_2_data,
			_mask_2_data,
			_fb_data_2,
			code_2_data,
			mask_2_data,
			fb_data_2,
			width,
			height,
			width_step,
			theta );
	unsigned int n = 0;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width_step; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			//pixels valides et pixel différents
			type n1 = ( fb_data_1[n_pixel] ) & ( _fb_data_2[n_pixel] ), 
				 c1 = ( code_1_data[n_pixel] ^ _code_2_data[n_pixel] ) & n1;
			
			//Comptage des bits
			for ( unsigned int k = 0; k < 8 * sizeof(type); ++ k )
			{
				n += (n1 >> k) & 1;
				d += (c1 >> k) & 1;
			}
		}
	}
	if ( n == 0 )
	{
		d = 1;
		return 1;
	}
	d /= n;
	
	
	return 0;
	
	
}


template <class type> int  distance :: fragile_bit_distance_opt ( 	double & d,
																	const type * code_1_data,
																	const type * mask_1_data,
																	const type * fb_data_1,
																	const type * code_2_data,
																	const type * mask_2_data,
																	const type * fb_data_2,
																	unsigned int width,
																	unsigned int height,
																	unsigned int width_step,
																	int theta,
																	const void * params,
																	type * _code_2_data,
																	type * _mask_2_data,
																	type * _fb_data_2 )
{
	d = 0;

	rotate( _code_2_data,
			_mask_2_data,
			_fb_data_2,
			code_2_data,
			mask_2_data,
			fb_data_2,
			width,
			height,
			width_step,
			theta );
			
	unsigned int n = 0;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width_step; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			//pixels valides et pixel différents
			type n1 = mask_1_data[n_pixel] & _mask_2_data[n_pixel], 
				 c1 = ( fb_data_1[n_pixel] & _fb_data_2[n_pixel] ) & n1;
			
			//Comptage des bits
			for ( unsigned int k = 0; k < 8 * sizeof(type); ++ k )
			{
				n += (n1 >> k) & 1;
				d += (c1 >> k) & 1;
			}
		}
	}
	if ( n == 0 )
	{
		d = 1;
		return 1;
	}
	d /= n;
	d = 1 - d;
	
	return 0;
	
}
														

template <class type> int  distance :: Hamming_FBD_opt ( 	double & d,
															const type * code_1_data,
															const type * mask_1_data,
															const type * fb_data_1,
															const type * code_2_data,
															const type * mask_2_data,
															const type * fb_data_2,
															unsigned int width,
															unsigned int height,
															unsigned int width_step,
															int theta,
															const void * params,
															type * _code_2_data,
															type * _mask_2_data,
															type * _fb_data_2 )
{
	double alpha = *((double*) params);
	
	
	
	
	d = 0;
	double d1 = 0,d2 = 0;
	rotate( _code_2_data,
			_mask_2_data,
			_fb_data_2,
			code_2_data,
			mask_2_data,
			fb_data_2,
			width,
			height,
			width_step,
			theta );
	
	unsigned int n = 0;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width_step; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			//pixels valides et pixel différents
			type n1 = mask_1_data[n_pixel] & _mask_2_data[n_pixel], 
				 c1 = ( fb_data_1[n_pixel] & _fb_data_2[n_pixel] ) & n1,
				 c2 = ( code_1_data[n_pixel] ^ _code_2_data[n_pixel] ) & n1;
			//Comptage des bits
			for ( unsigned int k = 0; k < 8 * sizeof(type); ++ k )
			{
				n += (n1 >> k) & 1;
				d1 += (c1 >> k) & 1;
				d2 += (c2 >> k) & 1;
			}
		}
	}
	if ( n == 0 )
	{
		d = 1;
		return 1;
	}
	d1 /= n;
	d2 /= n;
	
	d = alpha * (1 - d1) + (1 - alpha) * d2;
	return 0;
	
	
}







template <class type> int distance :: convert(	type * code_out,
												type * mask_out,
												type * fb_data_out,
												unsigned int & width_step,
												const IplImage * image_1,
												const IplImage * mask_1,
												unsigned int fm_threshold )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;
	
	if ( !mask_1 || !image_1 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )
	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )	
	if (	image_1->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U  	)
		return  1;
		
	return distance :: convert( code_out,
								mask_out,
								fb_data_out,
								width_step,
								(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
								(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
								mask1_width,
								mask1_height,
								img1_width_step,
								mask1_width_step,
								fm_threshold );
	
}






int distance :: Hamming ( double & d,
						  const IplImage * image_1,
						  const IplImage * mask_1,
						  const IplImage * image_2,
						  const IplImage * mask_2,
						  int theta,
						  const void * params )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	img2_x_offset,
					img2_y_offset,
					img2_width,
					img2_height,
					img2_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;

	unsigned int 	mask2_x_offset,
					mask2_y_offset,
					mask2_width,
					mask2_height,
					mask2_width_step;
	
	
	if ( !mask_1 || !mask_2 || !image_1 || !image_2 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )

	 GET_IMAGE_DIM( image_2,
					img2_width,
					img2_height,
					img2_width_step,
					img2_x_offset,
					img2_y_offset )

	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )

	 GET_IMAGE_DIM( mask_2,
					mask2_width,
					mask2_height,
					mask2_width_step,
					mask2_x_offset,
					mask2_y_offset )
	if ( 	img1_width != img2_width 		||
			img1_width != mask1_width		||
			img1_width != mask2_width		||
			img1_height != img2_height 		||
			img1_height < mask1_height		||
			img1_height < mask2_height		)
	{
		return 1;
	}
	if (
			image_1->depth != IPL_DEPTH_8U 	||
			image_2->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U 	||
			mask_2->depth != IPL_DEPTH_8U 	)
		return  1;
		
	return distance :: Hamming (	d,
									(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
									(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
									(unsigned char*) image_2->imageData + img2_y_offset * img2_width_step + img2_x_offset,
									(unsigned char*) mask_2->imageData + mask2_y_offset * mask2_width_step + mask2_x_offset,
									img1_width,
									img1_height,
									img1_width_step,
									mask1_width_step,
									img2_width_step,
									mask2_width_step,
									theta,
									params );
}


int distance :: Hamming ( 	double & d,
							const unsigned char * code_1_data,
							const unsigned char * fragility_map_1_data,
							const unsigned char * code_2_data,
							const unsigned char * fragility_map_2_data,
							unsigned int width,
							unsigned int height,
							unsigned int code_1_width_step,
							unsigned int fragility_map_1_width_step,
							unsigned int code_2_width_step,
							unsigned int fragility_map_2_width_step,
							int theta,
							const void * params )
{
	d = 0;
	unsigned int n_pixel = 0;
	const Hamming_parameters * tmp = (const Hamming_parameters *) params;
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int 	n_pixel_code_1 = i * code_1_width_step + j,
							n_pixel_code_2 = i * code_2_width_step + ( ( ( (int) j ) + theta ) % width),
							n_pixel_map_1 = i * fragility_map_1_width_step + j,
							n_pixel_map_2 = i * fragility_map_2_width_step + ( ( ( (int) j ) + theta ) % width);

			//Test de la validté du pixel
			if ( 	( fragility_map_1_data[n_pixel_map_1] >= tmp->fragile_bit_threshold_1 ) 	&& 
					( fragility_map_2_data[n_pixel_map_2] >= tmp->fragile_bit_threshold_2 ) 	)
			{
				n_pixel ++;
				if ( code_1_data[n_pixel_code_1] != code_2_data[n_pixel_code_2] )
					d ++;
			}
		}
	}
	if ( n_pixel == 0 )
	{
		d = 1;
		return 1;
	}

	
	d /= n_pixel;

	return 0;
}


int distance :: fragile_bit_distance ( double & d,
									   const IplImage * image_1,
									   const IplImage * mask_1,
									   const IplImage * image_2,
									   const IplImage * mask_2,
									   int theta,
									   const void * params )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	img2_x_offset,
					img2_y_offset,
					img2_width,
					img2_height,
					img2_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;

	unsigned int 	mask2_x_offset,
					mask2_y_offset,
					mask2_width,
					mask2_height,
					mask2_width_step;
	
	
	if ( !mask_1 || !mask_2 || !image_1 || !image_2 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )

	 GET_IMAGE_DIM( image_2,
					img2_width,
					img2_height,
					img2_width_step,
					img2_x_offset,
					img2_y_offset )

	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )

	 GET_IMAGE_DIM( mask_2,
					mask2_width,
					mask2_height,
					mask2_width_step,
					mask2_x_offset,
					mask2_y_offset )
	if ( 	img1_width != img2_width 		||
			img1_width != mask1_width		||
			img1_width != mask2_width		||
			img1_height != img2_height 		||
			img1_height < mask1_height		||
			img1_height < mask2_height		)
	{
		return 1;
	}
	if (
			image_1->depth != IPL_DEPTH_8U 	||
			image_2->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U 	||
			mask_2->depth != IPL_DEPTH_8U 	)
		return  1;
		
	return distance :: fragile_bit_distance (	d,
												(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
												(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
												(unsigned char*) image_2->imageData + img2_y_offset * img2_width_step + img2_x_offset,
												(unsigned char*) mask_2->imageData + mask2_y_offset * mask2_width_step + mask2_x_offset,
												img1_width,
												img1_height,
												img1_width_step,
												mask1_width_step,
												img2_width_step,
												mask2_width_step,
												theta,
												params );
}

int distance :: fragile_bit_distance ( 	double & d,
										const unsigned char * code_1_data,
										const unsigned char * fragility_map_1_data,
										const unsigned char * code_2_data,
										const unsigned char * fragility_map_2_data,
										unsigned int width,
										unsigned int height,
										unsigned int code_1_width_step,
										unsigned int fragility_map_1_width_step,
										unsigned int code_2_width_step,
										unsigned int fragility_map_2_width_step,
										int theta,
										const void * params )
{
	d = 0;
	unsigned int n_pixel = 0;

	const Hamming_parameters * tmp = (const Hamming_parameters *) params;

	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int 	n_pixel_map_1 = i * fragility_map_1_width_step + j,
							n_pixel_map_2 = i * fragility_map_2_width_step + ( ( j + theta ) % width );

			//Test du pixel ( en dehors de la zone d'occlusion )
			if ( fragility_map_1_data[n_pixel_map_1] > 0 &&
				 fragility_map_2_data[n_pixel_map_2] > 0 )
			{
				n_pixel ++;
				if ( fragility_map_1_data[n_pixel_map_1] < tmp->fragile_bit_threshold_1 xor
					 fragility_map_2_data[n_pixel_map_2] < tmp->fragile_bit_threshold_2 )
					d ++; 
			}
		}
	}
	if ( n_pixel == 0 )
	{
		d = 1;
		return 1;
	}
	d /= n_pixel;
	return 0;
}

int distance :: Hamming_expectation ( double & d,
									  const IplImage * image_1,
									  const IplImage * mask_1,
									  const IplImage * image_2,
									  const IplImage * mask_2,
									  int theta,
									  const void * params )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	img2_x_offset,
					img2_y_offset,
					img2_width,
					img2_height,
					img2_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;

	unsigned int 	mask2_x_offset,
					mask2_y_offset,
					mask2_width,
					mask2_height,
					mask2_width_step;
	
	
	if ( !mask_1 || !mask_2 || !image_1 || !image_2 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )

	 GET_IMAGE_DIM( image_2,
					img2_width,
					img2_height,
					img2_width_step,
					img2_x_offset,
					img2_y_offset )

	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )

	 GET_IMAGE_DIM( mask_2,
					mask2_width,
					mask2_height,
					mask2_width_step,
					mask2_x_offset,
					mask2_y_offset )
	if ( 	img1_width != img2_width 		||
			img1_width != mask1_width		||
			img1_width != mask2_width		||
			img1_height != img2_height 		||
			img1_height < mask1_height		||
			img1_height < mask2_height		)
	{
		return 1;
	}
	if (
			image_1->depth != IPL_DEPTH_8U 	||
			image_2->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U 	||
			mask_2->depth != IPL_DEPTH_8U 	)
		return  1;
		
	return distance :: Hamming_expectation (	d,
												(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
												(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
												(unsigned char*) image_2->imageData + img2_y_offset * img2_width_step + img2_x_offset,
												(unsigned char*) mask_2->imageData + mask2_y_offset * mask2_width_step + mask2_x_offset,
												img1_width,
												img1_height,
												img1_width_step,
												mask1_width_step,
												img2_width_step,
												mask2_width_step,
												theta,
												params );
}


int distance :: Hamming_expectation ( 	double & d,
										const unsigned char * code_1_data,
										const unsigned char * fragility_map_1_data,
										const unsigned char * code_2_data,
										const unsigned char * fragility_map_2_data,
										unsigned int width,
										unsigned int height,
										unsigned int code_1_width_step,
										unsigned int fragility_map_1_width_step,
										unsigned int code_2_width_step,
										unsigned int fragility_map_2_width_step,
										int theta,
										const void * params )
{
	d = 0;
	unsigned int n_pixel = 0;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int 	n_pixel_code_1 = i * code_1_width_step + j,
							n_pixel_code_2 = i * code_2_width_step + ( ( ( (int) j ) + theta ) % width),
							n_pixel_map_1 = i * fragility_map_1_width_step + j,
							n_pixel_map_2 = i * fragility_map_2_width_step + ( ( ( (int) j ) + theta ) % width);

			//Test de la validté du pixel
			if ( fragility_map_1_data[n_pixel_map_1] > 0 &&
				 fragility_map_2_data[n_pixel_map_2] > 0 )
			{
				n_pixel ++;
				double p1 = 0.5 * (255.0 - fragility_map_1_data[n_pixel_map_1]) / 255.0,
					   q1 = 1 - p1,
					   p2 = 0.5 * (255.0 - fragility_map_2_data[n_pixel_map_2]) / 255.0,
					   q2 = 1 - p2;
				if ( code_1_data[n_pixel_code_1] == code_2_data[n_pixel_code_2] )
				{
					d += p1 * q2;
					d += q1 * p2;
				}
				else 
				{
					d += q1 * q2;
					d += p1 * p2;
				}				
			}
		}
	}
	if ( n_pixel == 0 )
	{
		d = 1;
		return 1;
	}

	d /= n_pixel;
	return 0;
}

int distance :: Hamming_FBD( 	double & d,
								const unsigned char * code_1_data,
								const unsigned char * fragility_map_1_data,
								const unsigned char * code_2_data,
								const unsigned char * fragility_map_2_data,
								unsigned int width,
								unsigned int height,
								unsigned int code_1_width_step,
								unsigned int fragility_map_1_width_step,
								unsigned int code_2_width_step,
								unsigned int fragility_map_2_width_step,
								int theta,
								const void * params )
{
	const Hamming_FBD_parameters * tmp = (const Hamming_FBD_parameters *) params;
	d = 0;
	double d_tmp = 0;
	if ( distance :: Hamming (	d_tmp,
								code_1_data,
								fragility_map_1_data,
								code_2_data,
								fragility_map_2_data,
								width,
								height,
								code_1_width_step,
								fragility_map_1_width_step,
								code_2_width_step,
								fragility_map_2_width_step,
								theta,
								(void*) & (tmp->p_Hamming) ) )
		return 1;
	d += tmp->alpha * d_tmp;
	if ( distance :: fragile_bit_distance(	d_tmp,
											code_1_data,
											fragility_map_1_data,
											code_2_data,
											fragility_map_2_data,
											width,
											height,
											code_1_width_step,
											fragility_map_1_width_step,
											code_2_width_step,
											fragility_map_2_width_step,
											theta,
											(void*) & (tmp->p_Hamming) ) )
		return 1;
	d += (1 - tmp->alpha) * d_tmp;
	return 0;
}

int distance :: Hamming_FBD ( double & d,
							  const IplImage * image_1,
							  const IplImage * mask_1,
							  const IplImage * image_2,
							  const IplImage * mask_2,
							  int theta,
							  const void * params )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	img2_x_offset,
					img2_y_offset,
					img2_width,
					img2_height,
					img2_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;

	unsigned int 	mask2_x_offset,
					mask2_y_offset,
					mask2_width,
					mask2_height,
					mask2_width_step;
	
	
	if ( !mask_1 || !mask_2 || !image_1 || !image_2 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )

	 GET_IMAGE_DIM( image_2,
					img2_width,
					img2_height,
					img2_width_step,
					img2_x_offset,
					img2_y_offset )

	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )

	 GET_IMAGE_DIM( mask_2,
					mask2_width,
					mask2_height,
					mask2_width_step,
					mask2_x_offset,
					mask2_y_offset )
	if ( 	img1_width != img2_width 		||
			img1_width != mask1_width		||
			img1_width != mask2_width		||
			img1_height != img2_height 		||
			img1_height < mask1_height		||
			img1_height < mask2_height		)
	{
		return 1;
	}
	if (
			image_1->depth != IPL_DEPTH_8U 	||
			image_2->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U 	||
			mask_2->depth != IPL_DEPTH_8U 	)
		return  1;
		
	return distance :: Hamming_FBD (	d,
										(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
										(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
										(unsigned char*) image_2->imageData + img2_y_offset * img2_width_step + img2_x_offset,
										(unsigned char*) mask_2->imageData + mask2_y_offset * mask2_width_step + mask2_x_offset,
										img1_width,
										img1_height,
										img1_width_step,
										mask1_width_step,
										img2_width_step,
										mask2_width_step,
										theta,
										params );
}

int distance :: registering (	double & d,
								int & theta,
								const IplImage * image_1,
								const IplImage * mask_1,
								const IplImage * image_2,
								const IplImage * mask_2,
								const void * params,
								distance :: function_prototype_bis dist,
								int theta_min,
								int theta_max )
{
	unsigned int 	img1_x_offset,
					img1_y_offset,
					img1_width,
					img1_height,
					img1_width_step;

	unsigned int 	img2_x_offset,
					img2_y_offset,
					img2_width,
					img2_height,
					img2_width_step;

	unsigned int 	mask1_x_offset,
					mask1_y_offset,
					mask1_width,
					mask1_height,
					mask1_width_step;

	unsigned int 	mask2_x_offset,
					mask2_y_offset,
					mask2_width,
					mask2_height,
					mask2_width_step;
	
	
	if ( !mask_1 || !mask_2 || !image_1 || !image_2 )
		return 1;


	 GET_IMAGE_DIM( image_1,
					img1_width,
					img1_height,
					img1_width_step,
					img1_x_offset,
					img1_y_offset )

	 GET_IMAGE_DIM( image_2,
					img2_width,
					img2_height,
					img2_width_step,
					img2_x_offset,
					img2_y_offset )

	 GET_IMAGE_DIM( mask_1,
					mask1_width,
					mask1_height,
					mask1_width_step,
					mask1_x_offset,
					mask1_y_offset )

	 GET_IMAGE_DIM( mask_2,
					mask2_width,
					mask2_height,
					mask2_width_step,
					mask2_x_offset,
					mask2_y_offset )
	if ( 	img1_width != img2_width 		||
			img1_width != mask1_width		||
			img1_width != mask2_width		||
			img1_height != img2_height 		||
			img1_height < mask1_height		||
			img1_height < mask2_height		)
	{
		return 1;
	}
	if (
			image_1->depth != IPL_DEPTH_8U 	||
			image_2->depth != IPL_DEPTH_8U 	||
			mask_1->depth != IPL_DEPTH_8U 	||
			mask_2->depth != IPL_DEPTH_8U 	)
		return  1;
	
	return distance :: registering( d,
									theta,
									(unsigned char*) image_1->imageData + img1_y_offset * img1_width_step + img1_x_offset,
									(unsigned char*) mask_1->imageData + mask1_y_offset * mask1_width_step + mask1_x_offset,
									(unsigned char*) image_2->imageData + img2_y_offset * img2_width_step + img2_x_offset,
									(unsigned char*) mask_2->imageData + mask2_y_offset * mask2_width_step + mask2_x_offset,
									img1_width,
									img1_height,
									img1_width_step,
									mask1_width_step,
									img2_width_step,
									mask2_width_step,
									params,
									dist,
									theta_min,
									theta_max
									);
	
	
	
	
}

int distance :: registering ( 	double & d,
								int & theta,
								const unsigned char * code_1_data,
								const unsigned char * fragility_map_1_data,
								const unsigned char * code_2_data,
								const unsigned char * fragility_map_2_data,
								unsigned int width,
								unsigned int height,
								unsigned int code_1_width_step,
								unsigned int fragility_map_1_width_step,
								unsigned int code_2_width_step,
								unsigned int fragility_map_2_width_step,
								const void * params,
								distance :: function_prototype_bis dist,
								int theta_min,
								int theta_max )
{
	d = 2;
	theta = theta_min;

	for( int i = theta_min; i < theta_max; ++ i )
	{

		double v = 2;
		if (	dist(	v,
						code_1_data,
						fragility_map_1_data,
						code_2_data,
						fragility_map_2_data,
						width,
						height,
						code_1_width_step,
						fragility_map_1_width_step,
						code_2_width_step,
						fragility_map_2_width_step,
						i,
						params ) == 0)
		{		
			if ( v < d )
			{
				theta = i;
				d = v;
			}
		}
	}

	if ( d == 2 )
		return 1;
	return 0;
}

int distance :: registering_64b(	double & d,
									int theta,
									const unsigned long long * code_1_data,
									const unsigned long long * mask_1_data,
									const unsigned long long * fb_1_data,
									const unsigned long long * code_2_data,
									const unsigned long long * mask_2_data,
									const unsigned long long * fb_2_data,
									unsigned int width,
									unsigned int height,
									unsigned int width_step,
									int theta_min,
									int theta_max,
									function_prototype_64b dist,
									const void * params,
									unsigned long long * _code_2_data,
									unsigned long long * _mask_2_data,
									unsigned long long * _fb_2_data )
{
	d = 2;
	theta = theta_min;
	for( int i = theta_min; i < theta_max; ++ i )
	{

		double v = 2;
		if (	dist(	v,
						code_1_data,
						mask_1_data,
						fb_1_data,
						code_2_data,
						mask_2_data,
						fb_2_data,
						width,
						height,
						width_step,
						i,
						params,
						_code_2_data,
						_mask_2_data,
						_fb_2_data ) == 0)
		{		
			if ( v < d )
			{
				theta = i;
				d = v;
			}
		}
	}

	if ( d == 2 )
		return 1;
	return 0;
}
DIST_ALL(unsigned long long)
DIST_ALL(unsigned char)
DIST_ALL(unsigned short int)
DIST_ALL(unsigned long int)

