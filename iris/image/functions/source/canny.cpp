#include "canny.hpp"
#include <iostream>
#include <opencv/highgui.h>
using namespace std;
int canny ( IplImage * grad,
			IplImage * ori,
			const IplImage * image_in,
			double sigma,
			double scaling,
			double vert,
			double horz )
{
	IplImage * tmp_img;
	int width,
		height;
	if (  image_in->roi )
	{
		tmp_img = cvCreateImage( cvSize( 	image_in->roi->width,
											image_in->roi->height ),
									IPL_DEPTH_64F,
									1 );
		width = image_in->roi->width / scaling + 0.5;
		height = image_in->roi->height / scaling + 0.5;
	}
	else
	{
		tmp_img = cvCreateImage( cvSize( 	image_in->width,
											image_in->height ),
									IPL_DEPTH_64F,
									1 );
		width = image_in->width / scaling + 0.5;
		height = image_in->height / scaling + 0.5;
	}
	


	//Lissage de l'image
	cvSmooth(	image_in, 
				tmp_img,
				CV_GAUSSIAN,
				2 * ( (int) ( ( sigma - 0.8 ) / 0.3 + 1 ) ) + 1, 
				2 * ( (int) ( ( sigma - 0.8 ) / 0.3 + 1 ) ) + 1 );/// @todo

	//Redimensionnement
	if ( scaling != 1 )
	{	
		IplImage * tmp_img_2 = cvCreateImage ( 	cvSize( width,
														height ),
												IPL_DEPTH_64F,
												1 );
		//Redimensionnement de l'image
		cvResize( tmp_img, tmp_img_2, 0 );

		cvReleaseImage(&tmp_img);
		tmp_img = tmp_img_2;
	}
	IplImage * tmp_img_2;
	tmp_img_2 = cvCreateImage ( cvSize ( tmp_img->width,
										 tmp_img->height ),
								IPL_DEPTH_32F,
								1 );
	for ( unsigned int i = 0; i < (unsigned int) tmp_img_2->height; ++ i )
	{
		for ( unsigned int j = 0; j < (unsigned int) tmp_img_2->width; ++ j )
		{
			( (float*) ( tmp_img_2->imageData + i * tmp_img_2->widthStep ) )[j] = ( (double*) ( tmp_img->imageData + i * tmp_img->widthStep ) )[j];
		}
	}
	
	IplImage * img_1, 
			 * img_2;
			 
	img_1 = cvCreateImage ( cvSize ( tmp_img->width,
									 tmp_img->height ),
							IPL_DEPTH_32F,
							1 );
							
	img_2 = cvCreateImage ( cvSize ( tmp_img->width,
									 tmp_img->height ),
							IPL_DEPTH_32F,
							1 );
	//Grad.
	cvSobel( 	tmp_img_2, 
				img_1, 
				1, 
				0, 
				3 );
	cvSobel( 	tmp_img_2, 
				img_2, 
				0, 
				1, 
				3 );
	//
	for ( unsigned int i = 0; i < (unsigned int) tmp_img_2->height; ++ i )
	{
		for ( unsigned int j = 0; j < (unsigned int) tmp_img_2->width; ++ j )
		{
			double x = ( (float*) img_1->imageData )[ i * img_1->widthStep / sizeof(float) + j] * horz,
				   y = ( (float*) img_2->imageData )[ i * img_2->widthStep / sizeof(float) + j] * vert;
			double angle;   
			( ( double*) grad->imageData )[ i * grad->widthStep / sizeof(double) + j] = sqrt ( x * x + y * y );
			angle = atan2( -y, x);
			if ( angle < 0 )
				angle += M_PI;
			( ( double*) ori->imageData )[ i * ori->widthStep / sizeof(double) + j] = angle;   
		}
	}
	
	cvReleaseImage ( &tmp_img );
	cvReleaseImage ( &tmp_img_2 );
	cvReleaseImage ( &img_1 );
	cvReleaseImage ( &img_2 );
	return 0;
	
}
			
