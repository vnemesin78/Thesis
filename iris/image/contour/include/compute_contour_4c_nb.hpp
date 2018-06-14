/**\file compute_contour_4c.hpp
 * \author Valérian Némesin
 * \brief Ce fichier contient le prototype de la fonction qui détermnine l'ensemble des pixels de la frontière 4-connexe.
 */
#ifndef COMPUTECONTOUR4C_NB_HPP
#define COMPUTECONTOUR4C_NB_HPP
/**\fn template <class T> int compute_contour_4c(unsigned int * contour,
												 unsigned int * pNbPointsContour,
												 unsigned int maxPointsContour,
												 const T * pixels,
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step,
												 const T & value);		
 * \param[out] contour : contour calculé (Doit être préalloué)
 * \param[out] pNbPointsContour : pointeur vers le nombre de points du contour.
 * \param[out] pNbPointsContour : pointeur vers le nombre de points du contour.
 * \param[in] maxPointsContour : nombre de points maximal du contour.
 * \param[in] pixels : pixels de l'image
 * \param[in] width : largeur de l'image
 * \param[in] height : hauteur de l'image
 * \param[in] height : hauteur de l'image
 * @param[in] width_step : largeur de la ligne
 * \return La fonction renvoie 0 si le calcul du contour fonctionne.
 * \brief La fonction liste l'ensemble des points de la frontière 4 - connexe du masque. \n
 Explications: \n
 \t A) La fonction détermine le premier point blanc à partir du coin en haut à gauche du masque. \n
 \t B) Tant que la fonction ne retombe pas sur le second point dans la même direction, elle effectue cet ordre: \n
 \t \t 1) La fonction regarde dans les directions dans cet ordre : gauche, bas, droite, haut. \n
 \t \t 2) Si le pixel observé est blanc, alors c'est un pixel du contour, donc la fonction ajoute le point à la liste et se repositionne sur ce point. \n
 \t \t 3) La fonction décale la position de 3 cases pour mettre la posibilité de demi-tour en dernière. \n
 * \test Analyse d'un masque simple: \n
 \t Soit le masque à analyser: \n
 \t  xxx \n
 \t  x x \n
 \t xxxxx \n
 \t La fonction nous renvoie la liste de points: \n
 \t \t (0,1) (0,2) (0,3) (1,3) (2,3) (2,4) (2,3) (2,2) (2,1) (2,0) (2,1) (1,1) (0,1) (0,2)fin \n
 */
template <class T> int compute_contour_4c_nb(unsigned int * contour,
										  unsigned int * pNbPointsContour,
										  unsigned int maxPointsContour,
										  const T * pixels,
									      unsigned int width,
										  unsigned int height,
										  unsigned int width_step,
										  const T & value = 0);	
#define COMPUTE_CONTOUR_4C_NB(type) template int compute_contour_4c_nb(  unsigned int * contour,\
																		 unsigned int * pNbPointsContour,\
																		 unsigned int maxPointsContour,\
																		 const type * pixels,\
																		 unsigned int width,\
																		 unsigned int height,\
																		 unsigned int width_step,\
																		 const type & value); 
#endif
