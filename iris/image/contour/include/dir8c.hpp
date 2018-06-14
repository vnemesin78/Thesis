/**\file dir8c.hpp
 * \author Valérian Némesin
 * \date 27/12/2010
 * \brief Ce fichier contient les définitions des directions de recherche pour le contour 8-connexe
 */
#ifndef _DIR8C_HPP_
#define _DIR8C_HPP_
/**\var xDir8c
 * \brief La variable contient les abscisses des directions de recherche pour la recherche du contour 8-connexe.
**/
const int xDir8c[8] = {1, 
					   1, 
					   0,
					   -1,
					   -1,
					   -1, 
					   0, 
					   1};
				 
				 
/**\var yDir8c
 * \brief La variable contient les ordonnées des directions de recherche pour la recherche du contour 4-connexe.
**/			 
const int yDir8c[8] = {0, 
					   1, 
					   1, 
					   1, 
					   0,
					   -1,
					   -1,
					   -1};

#endif
