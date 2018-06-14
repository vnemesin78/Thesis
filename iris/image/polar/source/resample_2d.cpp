#include "resample_2d.hpp"
#include <iostream>
using namespace std;
template <class type> void resample_2d(	type & p_value,
											unsigned char & p_mask,
											const double & x,
											const double & y,
											const double & dx,
											const double & dy,
											const type * img_data,
											const unsigned char * img_mask,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											unsigned int mask_width_step )
{
	double sum_coef = 0,
		   coef = 0,
		   value = 0,
		   x_end,
		   y_end;
	unsigned int i, 
				 j, 
				 i_end, 
				 j_end;

	//Mise à zéro de la valeur
	value = 0;
	
	//Limites
	x_end = x + dx;
	y_end = y + dy;
	i_end = (int) (y_end);
	j_end = (int) (x_end);
	
	//Première ligne
	{
		j = (int) x;
		i = (int) y;
		//Première colonne
		if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
		{
			coef = ( 1 - (x - j) ) * ( 1 - (y - i) );
			sum_coef += coef;
			
			value += coef * img_data[ i * width_step + j ];
		}
		//Autres colonnes
		for ( ++ j; j < j_end; ++ j )
		{
			if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
			{
				coef = ( 1 - (y - i) );
				sum_coef += coef;
				value += coef * img_data[ i * width_step + j ];
			}			
		}
		//Dernière colonne
		if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
		{

			coef = ( x_end - j ) * ( 1 - (y - i) );
			sum_coef += coef;
			
			value += coef * img_data[ i * width_step + j ];
		}
		
		
	}
	//Autres lignes
	{

		for ( ++ i; i < i_end; ++ i )
		{

			j = (int) x;
			//Première colonne
			if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
			{
				coef = ( 1 - (x - j) );
				sum_coef += coef;
				
				value += coef * img_data[ i * width_step + j ];
			}
			//Autres colonnes
			for ( ++ j; j < j_end; ++ j )
			{
				if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
				{

					sum_coef += 1;
					value += img_data[ i * width_step + j ];
				}			
			}
			//Dernière colonne
			if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
			{
				coef = ( x_end - j );
				
				sum_coef += coef;
				value += coef * img_data[ i * width_step + j ];
			}
		}
	}
	//Dernière ligne
	{

		j = (int) x;
		//Première colonne
		if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
		{
			
			coef = ( y_end - i ) * ( 1 - (x - j) );
			sum_coef += coef;
				
			value += coef * img_data[ i * width_step + j ];
		}
	
		//Autre colonnes
		for ( ++ j; j < j_end; ++ j )
		{
			if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
			{

				coef = ( y_end - i );
				
				sum_coef += coef;
				value += coef * img_data[ i * width_step + j ];
			}			
		}
		//Dernière colonne
		if ( check_pixel( j, i, img_mask, width, height, mask_width_step ) )
		{

			coef = ( y_end - i ) * ( x_end - j );
			sum_coef += coef;
			value += coef * img_data[ i * width_step + j ];
		}
		
	
	
	}
	
	//Vérification
	if ( sum_coef == 0 ) //Pixel invalide!
	{
		p_mask = 0x00000000; // Porc!
		p_value = 0x00000000; // Voleur de poules!
	}
	else
	{

		p_mask = 255; //Saligot!
		p_value = (type) ( value / sum_coef );
	}
	return;
}

RESAMPLE_2D_DEF(unsigned char)
RESAMPLE_2D_DEF(char)

RESAMPLE_2D_DEF(unsigned short int)
RESAMPLE_2D_DEF(short int)

RESAMPLE_2D_DEF(unsigned long int)
RESAMPLE_2D_DEF(long int)

RESAMPLE_2D_DEF(double)
RESAMPLE_2D_DEF(float)

