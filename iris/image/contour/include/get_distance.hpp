/**\file get_distance.hpp
 * \author Valérian Némesin
 * \brief Le fichier contient la définition de la fonction inline qui calcule la distance entre deux points.
 */
#ifndef _GET_DISTANCE_HPP_
#define _GET_DISTANCE_HPP_

/**\fn inline double get_distance(unsigned int x0,
								  unsigned int y0,
								  unsigned int x1,
								  unsigned int y1)
 * \param[in] x0 : abscisse du premier point
 * \param[in] y0 : ordonnée du premier point
 * \param[in] x1 : abscisse du second point
 * \param[in] y1 : ordonnée du second point*
 * @brief
 * Calcule la distance euclidienne entre 2 points.
 */
template<class type> double get_distance(  const type & x0,
										   const type & y0,
										   const type & x1,
										   const type & y1);

#endif
