/**\file compute_contour_8c.hpp
 * \author Valérian Némesin
 * \date 27/12/2010
 * \brief Ce fichier contient le prototype de la fonction qui détermnine l'ensemble des pixels de la frontière 8-connexe.
 */
#ifndef COMPUTECONTOUR8C_H
#define COMPUTECONTOUR8C_H
/**\fn template <class T> int compute_contour_8c(unsigned int * contour,
												 unsigned int * pNbPointsContour,
												 unsigned int maxPointsContour,
												 const T * pixels,
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step,
												 const T & value)
 * \param[out] contour : contour calculé (Doit être préalloué)
 * \param[out] pNbPointsContour : pointeur vers le nombre de points du contour.
 * \param[in] maxPointsContour : nombre de points maximal du contour.
 * \param[in] pixels : pixels de l'image
 * \param[in] width : largeur de l'image
 * \param[in] height : hauteur de l'image
 * @param[in] width_step : largeur de la ligne
 * @param value : valeur du masque
 * \return La fonction renvoie 0 si le calcul du contour fonctionne.
 * \brief La fonction liste l'ensemble des points de la frontière 8 - connexe du masque. \n
    Explication: \n
     - On recherche le premier point du contour. \n
     - On utilise des directions dans cet ordre GAUCHE, BAS, HAUT, DROITE. \n
     - On regarde si le pixel est blanc dans l'une des directions. \n
     - Si il est blanc, ce pixel est un point du contour. \n
     - On décale de 5 cases la direction, ce qui permet de trouver un seul contour. \n
     - On s'arrète quand on tombe sur le second point et que la première direction de recherche est identique à la première. \n
 */
template <class T> int compute_contour_8c(unsigned int * contour,
										  unsigned int * pNbPointsContour,
										  unsigned int maxPointsContour,
										  const T * pixels,
									      unsigned int width,
										  unsigned int height,
										  unsigned int width_step,
										  const T & value = 0);

#define COMPUTE_CONTOUR_8C(type) template int compute_contour_8c(unsigned int * contour,\
																 unsigned int * pNbPointsContour,\
																 unsigned int maxPointsContour,\
																 const type * pixels,\
																 unsigned int width,\
																 unsigned int height,\
																 unsigned int width_step,\
																 const type & value = 0);
#endif
