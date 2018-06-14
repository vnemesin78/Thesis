#include "display_functions.hpp"
#include <sstream>
#include <cmath>
using namespace std;
#define ALIGN_LEFT -1
#define ALIGN_CENTER 0
#define ALIGN_RIGHT 1
#define DANS_TON_CUL 1;
void display_text( 	IplImage * image,
					const char * text,
					const CvRect & rect,
					int h_align,
					int v_align,
					int font_size,
					int color )
{
	CvFont font;
	cvInitFont(  &font, 
				 CV_FONT_HERSHEY_SIMPLEX, 
				 font_size / 25.0, 
				 font_size / 25.0 );

	CvSize text_size;
	int baseline;
	cvGetTextSize( text, &font, &text_size, &baseline );


	CvPoint point;
	point.x = rect.x;
	point.y = rect.y;
	if ( h_align == ALIGN_LEFT )
		point.x = rect.x;
	else if ( h_align == ALIGN_RIGHT )
		point.x = rect.x + rect.width - text_size.width;
	else
		point.x = rect.x + (rect.width - text_size.width) / 2;
		
	if ( v_align == ALIGN_LEFT )
		point.y = rect.y + text_size.height;
	else if ( v_align == ALIGN_RIGHT )
		point.y = rect.y + rect.height;
	else
		point.y = rect.y + (rect.height + text_size.height) / 2;
	
	 cvPutText( image, 
				text, 
				point, 
				&font, 
				CV_RGB( (color / 65536) % 256, (color / 256) % 256, color % 256) );
}

void draw_parabol( IplImage * image,
				   double a,
				   double b,
				   double c,
				   double theta,
				   double x_iris,
				   double y_iris,
				   double r_iris,
				   double dx,
				   double dy,
				   int color,
				   int thickness )
{
	double x, y;
	for ( double t = - 1.5 * r_iris; t < 1.5 * r_iris; t = t + 1 )
	{
		x = x_iris + cos(theta) * t + sin(theta) * ( a * pow(( t ),2) + c );
		y = y_iris - sin(theta) * t + cos(theta) * ( a * pow(( t ),2) + c );
		
		x /= dx;
		y /= dy;
		
		cvCircle( image,
				  cvPoint(x,y),
				  thickness,
				  CV_RGB( (color / 65536) % 256, (color / 256) % 256, color % 256),
				  CV_FILLED );
	}
}

					
c_display_image_data :: c_display_image_data()
{
	initialize();
}

c_display_image_data :: c_display_image_data( 	unsigned int _img_width,
												unsigned int _img_height, 
												unsigned int _data_width,
												unsigned int _data_height,
												int _font_size,
												int _font_color )
{
	initialize();
	setup ( _img_width, 
			_img_height, 
			_data_width, 
			_data_height, 
			_font_size, 
			_font_color );
}

int c_display_image_data :: setup ( )
{
	free();
	initialize();
	return 0;
}

int c_display_image_data :: setup ( unsigned int _img_width,
									unsigned int _img_height, 
									unsigned int _data_width,
									unsigned int _data_height,
									int _font_size,
									int _font_color )
{
	setup();
	if ( 	! _img_width	||
			! _img_height	||
			! _data_width	||
			! _data_height	||
			! _font_size	)
		return DANS_TON_CUL;
	
	_seg_image = cvCreateImage( cvSize( _img_width, _img_height), IPL_DEPTH_8U, 3 );
	_data_image = cvCreateImage( cvSize( _data_width, _data_height), IPL_DEPTH_8U, 3 );
	
	font_color = _font_color;
	font_size = _font_size;
	
	return 0;
}
			
int c_display_image_data :: display( const image_data & data )
{
	y = 0;
	
	if ( ! data.img_ok )
		return DANS_TON_CUL;
	
	//Convertion de l'image
	{
		IplImage * tmp = 
			cvCreateImage( 	cvSize( data.image->width, data.image->height), 
							IPL_DEPTH_8U, 
							3);
		
		cvConvertImage(	data.image, 
						tmp, 
						CV_GRAY2BGR );
		cvResize( 	tmp, 
					_seg_image );
		cvReleaseImage(&tmp);
	
	}
	
	//Affichage des info
	memset( _data_image->imageData, 
			0, 
			_data_image->widthStep * _data_image->height );
	
	//Titre
	{
		display_text ( 	_data_image, 
						"Image data", 
						cvRect(0,y, _data_image->width, 1.5 * font_size + y ),
						ALIGN_CENTER,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;
	}
	
	
	//Image name
	{
		stringstream oss;
		oss << "Name : " << data.name; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
	}				
	
	
	//ID
	{
		stringstream oss;
		oss << "#" << data.frame_id; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
	}				
	

	//Score
	{
		stringstream oss;
		oss << "Q : " << ( (int) (100 * data.score) ) / 100.0 ; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
	}				

	return 0;
}

c_display_image_data :: ~c_display_image_data()
{
	free();
	initialize();
}

void c_display_image_data :: free()
{
	if ( _seg_image )
		cvReleaseImage( &_seg_image );
		
	if ( _data_image )
		cvReleaseImage (&_data_image );
	
}

void c_display_image_data :: initialize()
{
	_seg_image = 0;
	_data_image = 0;
	font_color = 0;
	font_size = 0;
	y = 0;
}

c_display_pupil_data :: c_display_pupil_data()
{
	initialize();
}

c_display_pupil_data :: c_display_pupil_data( 	unsigned int _img_width,
												unsigned int _img_height, 
												unsigned int _data_width,
												unsigned int _data_heigth,
												int _font_size,
												int _font_color,
												int _thickness,
												int _pupil_color )
{
	initialize();
	setup ( _img_width,
			_img_height, 
			_data_width,
			_data_heigth,
			_font_size,
			_font_color,
			_thickness,
			_pupil_color );
}

int c_display_pupil_data :: setup ( )
{
	free();
	initialize();
}

int c_display_pupil_data :: setup ( 	unsigned int _img_width,
										unsigned int _img_height, 
										unsigned int _data_width,
										unsigned int _data_height,
										int _font_size,
										int _font_color,
										int _thickness,
										int _pupil_color )
{
	free();
	initialize();
	
	d_image = new c_display_image_data;
	
	
	if ( d_image->setup(	_img_width,
							_img_height, 
							_data_width,
							_data_height,
							_font_size,
							_font_color ) )
		return DANS_TON_CUL;
		
	_pupil_img = cvCloneImage( d_image->seg_image() );	
	_data_image = d_image->data_image();
	font_color = _font_color;
	font_size = _font_size;
	pupil_color = _pupil_color;
	thickness = _thickness;
	
	return 0;
}


int c_display_pupil_data :: display( const pupil_data & data )
{
	double 	dx,
			dy;
				
			

	if ( ! data.seg_ok )
		return DANS_TON_CUL;
	
	if ( d_image->display( data._img_data ) )
		return DANS_TON_CUL;
		
	y = d_image->data_offset();
	cvCopyImage( d_image->seg_image(), _pupil_img);
	dx = data.smoothed_image->width / ( (double) _pupil_img->width );
	dy = data.smoothed_image->height / ( (double) _pupil_img->height );
	

	
	//Dessin du ROI
	{
		cvRectangle( 	_pupil_img, 
						cvPoint( 	data.roi.x * dx, 
									data.roi.y * dy), 
						cvPoint( 	(data.roi.x + data.roi.width) * dx, 
									(data.roi.y + data.roi.height) * dy ),
						CV_RGB( (pupil_color / 65536) % 256, (pupil_color / 256) % 256, pupil_color % 256),
						thickness ); 
	}				
	
	//Dessin de la pupille
	{
		cvEllipse( 	_pupil_img,
					cvPoint( 	data.x_pupil * dx,
								data.y_pupil * dy),
					cvSize( data.a_pupil * dx,
							data.b_pupil * dy ),
					180 / M_PI * data.theta_pupil,
					0,
					360,
					CV_RGB( (pupil_color / 65536) % 256, (pupil_color / 256) % 256, pupil_color % 256),
					thickness );
		
		cvLine ( 	_pupil_img, 
					cvPoint( 	(data.x_pupil - 1.5 * data.a_pupil * cos( data.theta_pupil )) * dx,
								(data.y_pupil - 1.5 * data.a_pupil * sin( data.theta_pupil )) * dy ),
					cvPoint( 	(data.x_pupil + 1.5 * data.a_pupil * cos( data.theta_pupil )) * dx,
								(data.y_pupil + 1.5 * data.a_pupil * sin( data.theta_pupil )) * dy ),		
								
					CV_RGB( (pupil_color / 65536) % 256, (pupil_color / 256) % 256, pupil_color % 256),
					thickness );
		cvLine ( 	_pupil_img, 
					cvPoint( 	(data.x_pupil - 1.5 * data.b_pupil * sin( data.theta_pupil )) * dx,
								(data.y_pupil + 1.5 * data.b_pupil * cos( data.theta_pupil )) * dy ),
					cvPoint( 	(data.x_pupil + 1.5 * data.b_pupil * sin( data.theta_pupil )) * dx,
								(data.y_pupil - 1.5 * data.b_pupil * cos( data.theta_pupil )) * dy ),		
					CV_RGB( (pupil_color / 65536) % 256, (pupil_color / 256) % 256, pupil_color % 256),
					thickness );						
								
	}				
	
	//Info
	//Titre
	{
		display_text ( 	_data_image, 
						"Pupil data", 
						cvRect(0,y, _data_image->width, 1.5 * font_size + y ),
						ALIGN_CENTER,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;

	}
	
	//ROI
	{
		stringstream oss;
		oss << "x_ROI : " << data.roi.x; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color )	;	
		y += 1.5 * font_size;	
		oss.str("");
		oss << "y_ROI : " << data.roi.y; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
		oss.str("");	
		oss << "w_ROI : " << data.roi.width; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;
		oss.str("");
		oss << "h_ROI : " << data.roi.height; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;
		oss.str("");
	}				
	
	//Pupille
	{
		stringstream oss;
		oss << "x_p : " << data.x_pupil; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
		oss.str("");
		oss << "y_p : " << data.y_pupil; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;		
		oss.str("");
		oss << "a_p : " << data.a_pupil; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;
		oss.str("");
		oss << "b_p : " << data.b_pupil; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;
		oss.str("");
		oss << "theta_p : " << data.theta_pupil; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;
		oss.str("");
		
	}		



	return 0;
}

c_display_pupil_data :: ~c_display_pupil_data()
{
	free();
	initialize();
}

void c_display_pupil_data :: free()
{
	if ( _pupil_img )
		cvReleaseImage( &_pupil_img );
		
	if ( d_image )
		delete d_image;
}

void c_display_pupil_data :: initialize()
{
	font_color = 0;
	font_size = 0;
	pupil_color = 0;
	thickness = 0;
	_pupil_img = 0;
	_data_image = 0;
	d_image = 0;
	y = 0;
}

c_display_iris_data :: c_display_iris_data()
{
	initialize();
}


c_display_iris_data :: c_display_iris_data( 	unsigned int _img_width,
												unsigned int _img_height, 
												unsigned int _data_width,
												unsigned int _data_heigth,
												int _font_size,
												int _font_color,
												int _thickness,
												int _pupil_color,
												int _color_upper_eyelid,
												int _color_lower_eyelid,
												int _color_iris,
												int _color_mask, //Nice try
												int _iris_code_width,
												int _iris_code_height,
												int _polar_img_width,
												int _polar_img_height )
{
	initialize();
	
	setup(	_img_width,
			_img_height, 
			_data_width,
			_data_heigth,
			_font_size,
			_font_color,
			_thickness,
			_pupil_color,
			_color_upper_eyelid,
			_color_lower_eyelid,
			_color_iris,
			_color_mask, //Nice try
			_iris_code_width,
			_iris_code_height,
			_polar_img_width,
			_polar_img_height );
			
	
	
	
	
}

int c_display_iris_data :: setup ( )
{
	free();
	initialize();
}

int c_display_iris_data :: setup ( 	unsigned int _img_width,
									unsigned int _img_height, 
									unsigned int _data_width,
									unsigned int _data_height,
									int _font_size,
									int _font_color,
									int _thickness,
									int _pupil_color,
									int _color_upper_eyelid,
									int _color_lower_eyelid,
									int _color_iris,
									int _color_mask,
									int _iris_code_width,
									int _iris_code_height,
									int _polar_img_width,
									int _polar_img_height )
{
	setup();
	
	if(	! _iris_code_width		|| 
		! _iris_code_height		|| 
		! _polar_img_width 		|| 
		! _polar_img_height 	)
		return DANS_TON_CUL;
		
	d_pupil = new c_display_pupil_data;
	if ( d_pupil->setup( 	_img_width,
							_img_height, 
							_data_width,
							_data_height,
							_font_size,
							_font_color,
							_thickness,
							_pupil_color ) )
		return DANS_TON_CUL;
		
	_iris_img = cvCloneImage( d_pupil->pupil_image() );	
	
	font_color = _font_color;
	font_size = _font_size;
	pupil_color = _pupil_color;
	thickness = _thickness;
	color_upper_eyelid = _color_upper_eyelid;
	color_lower_eyelid = _color_lower_eyelid;
	color_iris = _color_iris;
	color_mask = _color_mask;
	
	_data_image = d_pupil->data_image();
	
	_polar_img = cvCreateImage( cvSize( _polar_img_width, _polar_img_height ),
								IPL_DEPTH_8U,
								3 );
								
	_code_img = cvCreateImage( cvSize( _iris_code_width, _iris_code_height ),
								IPL_DEPTH_8U,
								3 );						
								
		
	return 0;
}

int c_display_iris_data :: display( const iris_data & data )
{

	
	if ( ! data.seg_ok )
		return DANS_TON_CUL;
		
	if ( d_pupil->display( data.p_data ) )
		return DANS_TON_CUL;
	y = d_pupil->data_offset();
	double 	dx,
			dy;

	dx = data.p_data.smoothed_image->width / ( (double) _iris_img->width );
	dy = data.p_data.smoothed_image->height / ( (double) _iris_img->height );
	
	cvCopyImage( d_pupil->pupil_image(), _iris_img );
	

	
	//Dessin de l'iris
	{
		cvEllipse( 	_iris_img,
					cvPoint( 	data.x_iris * dx,
								data.y_iris * dy),
					cvSize( data.a_iris * dx,
							data.b_iris * dy ),
					180 / M_PI * data.theta_iris,
					0,
					360,
					CV_RGB( (color_iris / 65536) % 256, (color_iris / 256) % 256, color_iris % 256),
					thickness );
		
		cvLine ( 	_iris_img,
					cvPoint( 	(data.x_iris - 1.5 * data.a_iris * cos( data.theta_iris )) * dx,
								(data.y_iris - 1.5 * data.a_iris * sin( data.theta_iris )) * dy ),
					cvPoint( 	(data.x_iris + 1.5 * data.a_iris * cos( data.theta_iris )) * dx,
								(data.y_iris + 1.5 * data.a_iris * sin( data.theta_iris )) * dy ),		
								
					CV_RGB( (color_iris / 65536) % 256, (color_iris / 256) % 256, color_iris % 256),
					thickness );
		cvLine ( 	_iris_img,
					cvPoint( 	(data.x_iris - 1.5 * data.b_iris * sin( data.theta_iris )) * dx,
								(data.y_iris + 1.5 * data.b_iris * cos( data.theta_iris )) * dy ),
					cvPoint( 	(data.x_iris + 1.5 * data.b_iris * sin( data.theta_iris )) * dx,
								(data.y_iris - 1.5 * data.b_iris * cos( data.theta_iris )) * dy ),		
					CV_RGB( (color_iris / 65536) % 256, (color_iris / 256) % 256, color_iris % 256),
					thickness );						
								
	}				

	//Info
	
	//Titre
	{
		display_text ( 	_data_image, 
						"Iris data", 
						cvRect(0,y, _data_image->width, 1.5 * font_size + y ),
						ALIGN_CENTER,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;
	}
	
	//Iris
	{
		stringstream oss;
		oss << "x_i : " << data.x_iris; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;	
		oss.str("");
		oss << "y_i : " << data.y_iris; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );		
		y += 1.5 * font_size;		
		oss.str("");
		oss << "a_i : " << data.a_iris; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );	
		y += 1.5 * font_size;
		oss.str("");
		oss << "b_i : " << data.b_iris; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );	
		y += 1.5 * font_size;
		oss.str("");
		oss << "theta_i : " << data.theta_iris; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;
		oss.str("");
	}		
	
	//Paupières
	{
		stringstream oss;
		oss << "a_upper : " << data.a_upper; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;	
		oss.str("");
		oss << "c_upper : " << data.c_upper; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;		
		oss.str("");
		oss << "theta_upper : " << data.theta_upper; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;
		oss.str("");
		oss << "a_lower : " << data.a_lower; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;	
		oss.str("");
		oss << "c_lower : " << data.c_lower; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;	
		oss.str("");	
		oss << "theta_lower : " << data.theta_lower; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		oss.str("");
		y += 1.5 * font_size;
	}		
	
	//Score
	{
		stringstream oss;
		oss << "score : " << data.nrj_ratio; 
		display_text ( 	_data_image, 
						oss.str().c_str(), 
						cvRect(font_size, y, _data_image->width - font_size, 1.5 * font_size + y ),
						ALIGN_LEFT,
						ALIGN_LEFT,
						font_size,
						font_color );
		y += 1.5 * font_size;		
	}
	
	
	//Dessin des paupières
	{
				
		draw_parabol( _iris_img, data.a_upper, data.x_iris, data.c_upper, data.theta_upper, data.x_iris, data.y_iris, data.new_r_iris, dx, dy, color_upper_eyelid, thickness );
		draw_parabol( _iris_img, data.a_lower, data.x_iris, data.c_lower, data.theta_lower, data.x_iris, data.y_iris, data.new_r_iris, dx, dy, color_lower_eyelid, thickness );
	}
	

	
	//Dessin de la transformée polaire
	{
		IplImage * tmp = 
			cvCreateImage( 	cvSize( data.polar_image->width, data.polar_image->height), 
							IPL_DEPTH_8U, 
							3);
		
		IplImage * tmp2 = 
			cvCreateImage( 	cvSize( data.polar_image->width, data.polar_image->height), 
							IPL_DEPTH_8U, 
							1);
		for ( int i = 0; i < data.polar_image->height; ++ i )
		{
			for ( int j = 0; j < data.polar_image->width; ++ j )
			{
				GET_IMAGE_PIXEL( tmp2, i, j, unsigned char ) = GET_IMAGE_PIXEL( data.polar_image, i, j, double );
				
			}
		}		
		cvConvertImage(	tmp2, 
						tmp, 
						CV_GRAY2BGR );
		cvReleaseImage(&tmp2);				
		cvResize( 	tmp, 
					_polar_img );
		cvReleaseImage(&tmp);	
		tmp = 
			cvCreateImage( 	cvSize( _polar_img->width, _polar_img->height), 
							IPL_DEPTH_8U, 
							1);
							
		cvResize( 	data.polar_mask, 
					tmp );
		
		for ( int i = 0; i < _polar_img->height; ++ i )
		{
			for ( int j = 0; j < _polar_img->width; ++ j )
			{
				if ( GET_IMAGE_PIXEL( tmp, i, j, unsigned char ) != 255 )
				{
					( (unsigned char *) ( _polar_img->imageData + _polar_img->widthStep * i ) )[ 3 * j ] = 0;
					( (unsigned char *) ( _polar_img->imageData + _polar_img->widthStep * i ) )[ 3 * j + 1 ] = 0;
					( (unsigned char *) ( _polar_img->imageData + _polar_img->widthStep * i ) )[ 3 * j + 2 ] = 255;
				}
				
			}
			
		}
		
		
		
		cvReleaseImage(&tmp);		
	}

	//Dessin de l'iris code
	{
		IplImage * tmp = 
			cvCreateImage( 	cvSize( data.code->width, data.code->height), 
							IPL_DEPTH_8U, 
							3);

		cvConvertImage(	data.code, 
						tmp, 
						CV_GRAY2BGR );

		
		cvResize( 	tmp, 
					_code_img );

		cvThreshold (  _code_img, _code_img, 128, 255, CV_THRESH_BINARY );
		cvReleaseImage(&tmp);	
		tmp = 
			cvCreateImage( 	cvSize( _code_img->width, _code_img->height / 2), 
							IPL_DEPTH_8U, 
							1);
							
		cvResize( 	data.code_mask, 
					tmp );
		
		for ( int i = 0; i < _code_img->height / 2; ++ i )
		{
			for ( int j = 0; j < _code_img->width; ++ j )
			{
				if ( GET_IMAGE_PIXEL( tmp, i, j, unsigned char ) != 255 )
				{
					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * i ) )[ 3 * j ] = 0;
					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * i ) )[ 3 * j + 1 ] = 0;
					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * i ) )[ 3 * j + 2 ] = 255;

					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * ( i + _code_img->height / 2 )  ) )[ 3 * j ] = 0;
					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * ( i + _code_img->height / 2 )  ) )[ 3 * j + 1 ] = 0;
					( (unsigned char *) ( _code_img->imageData + _code_img->widthStep * ( i + _code_img->height / 2 )  ) )[ 3 * j + 2 ] = 255;

				}
				
			}
			
		}
		cvReleaseImage(&tmp);		
	}

	return 0;
}

c_display_iris_data :: ~c_display_iris_data()
{
	free();
	initialize();
}

void c_display_iris_data :: free()
{
	if ( _iris_img )
		cvReleaseImage ( &_iris_img );
		
	if ( _polar_img )
		cvReleaseImage( &_polar_img );
	
	if ( _code_img )
		cvReleaseImage( &_code_img );
		
	if ( d_pupil )
		delete d_pupil;
}

void c_display_iris_data :: initialize()
{
	font_color = 0;
	font_size = 0;
	pupil_color = 0;
	thickness = 0;  
	color_upper_eyelid = 0;
	color_lower_eyelid = 0;
	color_iris = 0;
	color_mask = 0;
			
	_iris_img = 0;
	_polar_img = 0;
	_code_img = 0;
					 
	_data_image = 0;
			
	d_pupil = 0;
	y = 0;
}
