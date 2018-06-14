#include "c_preprocessing.hpp"
/**@fn
 * @brief
 * Constructeur
 **/
c_preprocessing :: c_preprocessing()
{
	initialize();
}

/**@fn 
* @param width : largeur maximale de travail
* @param height: hauteur maximale de travail
* @param opening : Ouverture dans le masque de la pupille
* @param erosion : Erosion dans le masque de l'iris
* @param median_filter : taille du filtre median pour l'image de la pupille
* @param err_stream : flux d'erreur ( NULL pour désactiver )
* @brief
* Constructeur
**/
c_preprocessing :: c_preprocessing( 	unsigned int width,
										unsigned int height,
										unsigned int closing_1,
										unsigned int opening_1,
										unsigned int closing_2,
										unsigned int opening_2,
										unsigned int median_filter,
										ostream * _err_stream )
{
	initialize();
	if ( setup( width, 
				height,
				closing_1,
				opening_1,
				closing_2,
				opening_2,
				median_filter,
				_err_stream ) )
		throw ( invalid_argument("Arguments of c_preprocessing!") );
}


/**@fn 
 * @param width : largeur maximale de travail
 * @param height: hauteur maximale de travail
 * @param opening : Ouverture dans le masque de la pupille
 * @param erosion : Erosion dans le masque de l'iris
 * @param median_filter : taille du filtre median pour l'image de la pupille
 * @param err_stream : flux d'erreur ( NULL pour désactiver )
 * @brief
 * Setup
 **/
int c_preprocessing :: setup( 	unsigned int width,
								unsigned int height,
								unsigned int closing_1,
								unsigned int opening_1,
								unsigned int closing_2,
								unsigned int opening_2,
								unsigned int median_filter,
								ostream * _err_stream )
{
	//Reset
	free();
	initialize();
	
	err_stream = _err_stream;
	

	//Arg
	_width_max = width;
	_height_max = height;
	_opening_1 = opening_1;
	_closing_1 = closing_1;
	_closing_2 = closing_2;
	_opening_2 = opening_2;
	_median_filter = median_filter;
	
	//Check des arg.
	if ( ! ( width && height ) )
	{
		if (err_stream)
		{
			*err_stream << "Error :\tInvalid image dimension in c_preprocessing :: setup" << endl;
			*err_stream << "\t\t\tWidth = " << width << endl;
			*err_stream << "\t\t\tHeight = " << height << endl;
		}
		return 1;
	}	
	
	alloc();
	
	return 0;
}

/**@fn 
 * @param width : largeur maximale de travail
 * @param height: hauteur maximale de travail
 * @param opening : Ouverture dans le masque de la pupille
 * @param erosion : Erosion dans le masque de l'iris
 * @param median_filter : taille du filtre median pour l'image de la pupille
 * @param err_stream : flux d'erreur ( NULL pour désactiver )
 * @brief
 * Constructeur
 **/
int c_preprocessing :: setup ( 	api_parameters & params,
								unsigned int width,
								unsigned int height,
								ostream * _err_stream,
								const char * n_space,
								const char * closing_1_name,
								const char * opening_1_name,
								const char * closing_2_name,
								const char * opening_2_name,
								const char * median_name
								)
{
	stringstream oss;
	int q = 0;
	//Reset
	free();
	initialize();
	
	//Flux d'errerus
	err_stream = _err_stream;
	
	unsigned int opening_1, closing_1, opening_2, closing_2;
	oss << n_space << "::" << opening_1_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&opening_1,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << closing_1_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&closing_1,
									err_stream ) )
		q = 1;
	oss.str("");
	oss << n_space << "::" << closing_2_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&closing_2,
									err_stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << opening_2_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&opening_2,
									err_stream ) )
		q = 1;
	oss.str("");
	
	
	unsigned int median_filter;
	oss << n_space << "::" << median_name;
	if (api_get_positive_integer(	params, 
									oss.str().c_str(), 
									&median_filter,
									err_stream ) )
		q = 1;
	oss.str("");
	
	if (q)
		return 1;
		
	return ( 	setup( 	width, 
						height,
						closing_1,
						opening_1,
						closing_2,
						opening_2,
						median_filter,
						_err_stream ) );
}

int c_preprocessing :: process ( const IplImage * image )
{		
	if ( image == NULL )
		return 1;
		
	//Vérification du format de l'image
	if ( image->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
			*err_stream << "Error in int c_preprocessing :: process ( const IplImage * image ) : Image format is not  IPL_DEPTH_8U !" << endl;
		return 1;
	}
	
	if ( 	GET_IMAGE_WIDTH(image) > _width_max 	||
			GET_IMAGE_HEIGHT(image) > _height_max	)
	{
			if (	setup( 	GET_IMAGE_WIDTH(image), 
							GET_IMAGE_HEIGHT(image),
							_closing_1,
							_opening_1,
							_closing_2,
							_opening_2,
							_median_filter ) )
				return 1;
	}
	
	
	//ROI
	cvSetImageROI( 	_smoothed_image,
					cvRect( 0, 
							0,
							GET_IMAGE_WIDTH(image),
							GET_IMAGE_HEIGHT(image) ) );
	
	
	//Filtre median
	if ( _median_filter )
		cvSmooth( image, 
				  _smoothed_image,
				  CV_MEDIAN,
				  _median_filter,
				  _median_filter );
				  
				  
	else
		cvCopyImage( 	image,
						_smoothed_image );

	//Suppression des cils ( Fermeture )
	if ( structuring_element_4 )
		cvMorphologyEx( _smoothed_image,
						_smoothed_image,
						NULL,
						structuring_element_4,
						CV_MOP_CLOSE,
						1);	
	

	//Suppression des spots (Ouverture)
	if ( structuring_element_1 )
		cvMorphologyEx( _smoothed_image,
						_smoothed_image,
						NULL,
						structuring_element_1,
						CV_MOP_OPEN,
						1);	
	
	//Suppression des ombres ( Fermeture)
	if ( structuring_element_2 )
		cvMorphologyEx( _smoothed_image,
						_smoothed_image,
						NULL,
						structuring_element_2,
						CV_MOP_CLOSE,
						1);	
						
	//Correction de la pupille ( Ouverture )
	if ( structuring_element_3 )
		cvMorphologyEx( _smoothed_image,
						_smoothed_image,
						NULL,
						structuring_element_3,
						CV_MOP_OPEN,
						1);	
	

	return 0;
}


c_preprocessing :: ~c_preprocessing()
{
	free();
	initialize();
}

void c_preprocessing :: free()
{
	if ( _smoothed_image )
		cvReleaseImage( &_smoothed_image);
	if ( structuring_element_1)
		cvReleaseStructuringElement( & structuring_element_1 );
	if ( structuring_element_2)
		cvReleaseStructuringElement( & structuring_element_2 );
	if ( structuring_element_3)
		cvReleaseStructuringElement( & structuring_element_3 );
	if ( structuring_element_4)
		cvReleaseStructuringElement( & structuring_element_4 );
}

void c_preprocessing :: initialize()
{
	err_stream = 0;
	_width_max = 0;
	_height_max = 0;
	_opening_1 = 0;
	_closing_1 = 0;
	_closing_2 = 0;
	_opening_2 = 0;
	_median_filter = 0;
	_width = 0;
	_height = 0;
	structuring_element_1 = 0;
	structuring_element_2 = 0;
	structuring_element_3 = 0;
	structuring_element_4 = 0;
	_smoothed_image = 0;
	
}

void c_preprocessing :: alloc()
{
	 _smoothed_image = cvCreateImage (	cvSize ( _width_max, _height_max),
										IPL_DEPTH_8U,
										1 );						
	if ( _opening_1 )
		structuring_element_1 = cvCreateStructuringElementEx( _opening_1 * 2 + 1,
															  _opening_1 * 2 + 1,
															  _opening_1,
															  _opening_1,
															  CV_SHAPE_ELLIPSE,
															  NULL);
	if ( _closing_2 )
		structuring_element_2 = cvCreateStructuringElementEx( _closing_2 * 2 + 1,
															  _closing_2 * 2 + 1,
															  _closing_2,
															  _closing_2,
															  CV_SHAPE_ELLIPSE,
															  NULL);	
															  
	if ( _opening_2 )
		structuring_element_3 = cvCreateStructuringElementEx( _opening_2 * 2 + 1,
															  _opening_2 * 2 + 1,
															  _opening_2,
															  _opening_2,
															  CV_SHAPE_ELLIPSE,
															  NULL);
	if ( _closing_1 )
		structuring_element_4 = cvCreateStructuringElementEx( _closing_1 * 2 + 1,
															  _closing_1 * 2 + 1,
															  _closing_1,
															  _closing_1,
															  CV_SHAPE_ELLIPSE,
															  NULL);	
														
															  
}


