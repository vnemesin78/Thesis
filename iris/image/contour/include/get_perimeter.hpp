/**\file get_perimeter.hpp
 * \author Valérian Némesin
 * \brief Ce fichier contient le prototype de la fonction qui calcule le périmètre d'un contour affine par morceau et continu.
 */
#ifndef _GET_PERIMETER_HPP_
#define _GET_PERIMETER_HPP_
/**\fn double get_perimeter(const unsigned int * contour,
							unsigned int nb_points_contour);
 * \param[in] contour : contour
 * \param[in] nb_points_contour : nombre de points du contour
 * \return Périmètre (double)
 * \brief Calcule le périmètre du contour.
 */
template <class type> double get_perimeter(	const type * contour,
											unsigned int nb_points_contour);
#endif
