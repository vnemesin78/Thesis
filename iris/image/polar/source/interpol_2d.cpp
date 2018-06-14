#include "interpol_2d.hpp"


/**@fn 
 * @param[out] p_value : valeur du pixel
 * @param[out] p_mask : masque de validité du pixel
 * @param[in] x : abscisse du pixel
 * @param[in] y : ordonnée du pixel
 * @param[in] img_data : données de l'image
 * @param[in] img_mask : masque de l'image
 * @param[in] width : largeur de l'image
 * @param[in] height : hauteur de l'image
 * @brief
 * Cette fonction effectue un filtrage billinéaire.
 * 
 */
template <class type> void interpol_2d ( 	double & p_value,
											unsigned char & p_mask,
											const double & x,
											const double & y,
											const type * img_data,
											const unsigned char * img_mask,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											unsigned int mask_width_step )
											
{
	type pixel_00 = 0,
		 pixel_01 = 0,
		 pixel_10 = 0,
		 pixel_11 = 0;
		 
	double coef_00 = 0,
		   coef_01 = 0,
		   coef_10 = 0,
		   coef_11 = 0,
		   sum_coef = 0;

	double dx,
		   dy;

	long long x_0,
			  x_1,
			  y_0,
			  y_1;
			  
	x_0 = (long long) x;
	x_1 = (long long) x + 1;

	y_0 = (long long) y;
	y_1 = (long long) y + 1;

	dx = x - x_0;
	dy = y - y_0;

	//Vérification des pixels
	if ( check_pixel( x_0, y_0, img_mask, width, height, mask_width_step ) )
	{
		pixel_00 = img_data[ x_0 + y_0 * width_step ];		
		coef_00 = ( 1 - dx ) * ( 1 - dy );		
		sum_coef += coef_00;
	}
	if ( check_pixel( x_1, y_0, img_mask, width, height, mask_width_step ) )
	{
		pixel_01 = img_data[ x_1 + y_0 * width_step ];		
		coef_01 = dx * ( 1 - dy );	
		sum_coef += coef_01;
	}

	if ( check_pixel( x_0, y_1, img_mask, width, height, mask_width_step ) )
	{
		pixel_10 = img_data[ x_0 + y_1 * width_step ];		
		coef_10 = ( 1 - dx ) * dy;	
		sum_coef += coef_10;
	}

	if ( check_pixel( x_1, y_1, img_mask, width, height, mask_width_step ) )
	{
		pixel_11 = img_data[ x_1 + y_1 * width_step ];		
		coef_11 = dx * dy;	
		sum_coef += coef_11;
	}

	//Validité de la transformée polaire
	if (sum_coef == 0)
	{
		p_value = 0x00000000;
		p_mask = 0x00000000;
	}
	else
	{
		p_mask = 255;
		p_value = (coef_00 * pixel_00 + coef_01 * pixel_01 + coef_10 * pixel_10 + coef_11 * pixel_11) / sum_coef;
	}

}
				   
INTERPOL_2D_DEF( unsigned char )
INTERPOL_2D_DEF( char )

INTERPOL_2D_DEF( unsigned short int )
INTERPOL_2D_DEF( short int )

INTERPOL_2D_DEF( unsigned long int )
INTERPOL_2D_DEF( long int )

INTERPOL_2D_DEF( float )
INTERPOL_2D_DEF( double )

















