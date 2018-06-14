/**\file digitalize_euclide_contour.hpp
 * \author Valérian Némesin
 * \brief Fichier qui contient le prototype de la fonction qui effectue un échantionnage euclidien du contour
 */
#ifndef _DIGITALIZE_EUCLIDE_HPP_
#define _DIGITALIZE_EUCLIDE_HPP_

/**\fn int digitalize_euclide_contour(double * d_contour,
									  unsigned int nb_points_d_contour,
									  const unsigned int * contour,
									  unsigned int nb_points_contour);
 * \param[out] d_contour : contour échantillonné
 * \param[in] nb_points_d_contour : nombre de points du contour échantillonné
 * \param[in] contour : contour à échantillonner
 * \param[in] nb_points_contour : nombre de points du contour à échantillonner
 * \return 
 * 0 si le ré-échantillonnage du contour a fonctionné
 * \brief La fonction échantillonne le contour. \n
 * Elle calcule des nouveaux points tels qu'ils soient équi-distant si on considère le contour comme étant affine par morceau et continu.
 */
template <class type> int digitalize_euclide_contour(  double * d_contour,
													   unsigned int nb_points_d_contour,
													   const type * contour,
													   unsigned int nb_points_contour);
#endif
