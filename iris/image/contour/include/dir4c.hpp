/**\file dir4c.hpp
 * \author Valérian Némesin
 * \date 27/12/2010
 * \brief Ce fichier contient les définitions des directions de recherche pour le contour 4-connexe
 */
#ifndef XDIR4C
#define XDIR4C
/**\var xDir4c
 * \brief La variable contient les abscisses des directions de recherche pour la recherche du contour 4-connexe.
**/
const int xDir4c[4] = {1,
					   0,
					   -1,
					   0};
				 
				 
/**\var yDir4c
 * \brief La variable contient les ordonnées des directions de recherche pour la recherche du contour 4-connexe.
**/			 
const int yDir4c[4] = {0,
					   1,
					   0,
					   -1};
/**\var x4c
 * \brief Position du contour en fonction de la frontière
 */
const int x4c[4] = {1, 
					1, 
					0, 
					0};
/**\var y4c
 * \brief Position du contour en fonction de la frontière
 */
const int y4c[4] = {1, 
					0, 
					1, 
					0};
#endif
