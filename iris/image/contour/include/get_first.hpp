/**\file get_first.h
 * \author Valérian Némesin
 * \brief Ce fichier contient la définition de la fonction qui recherche le premier point du masque
 */
#ifndef _GET_FIRST_HPP_
#define _GET_FIRST_HPP_
/**\fn template <class T> int get_first(unsigned int & x,
							     unsigned int & y,
							     const unsigned T * pixels,
							     unsigned int width,
							     unsigned int height,
							     unsigned int width_step,
							     const T & value);
 * \param[out] x : absisse du premier point du masque
 * \param[out] y : ordonnée du premier point du masque
 * \param[in] pixels : pixels du masque
 * \param[in] width : largeur du masque
 * \param[in] height : hauteur du masque
 * \param[in] width_step : taille de la ligne réelle
 * \param[in] value : valeur des pixels n'appartemant pas à l'arrière plan
 * \return 
 * 1 si le masque est vide
 * 0 sinon
 * \brief Cette fonction recherche le premier point du masque en partant du coin en haut à gauche.
 **/
template <class T> int get_first(unsigned int & x,
							     unsigned int & y,
							     const T * pixels,
							     unsigned int width,
							     unsigned int height,
							     unsigned int width_step,
							     const T & value = 0);

#define GET_FIRST(type) template int get_first(unsigned int & x,\
											   unsigned int & y,\
											   const type * pixels,\
											   unsigned int width,\
											   unsigned int height,\
											   unsigned int width_step,\
											   const type & value = 0);

#endif
