/**\file is_mask.hpp
 * \author Valérian Némesin
 */
#ifndef _IS_MASK_HPP_
#define _IS_MASK_HPP_
/**\fn template <class U> int is_mask(const unsigned int & x,
									 const unsigned int & y,
									 const T * pixels,
									 unsigned int width,
									 unsigned int height,
									 unsigned int width_step,
									 const T & value)
 * \param[in] x : abscisse
 * \param[in] y : ordonnée
 * \param[in] pixels : pixels de l'image
 * \param[in] width : largeur de l'image
 * \param[in] height : hauteur de l'image
 * \param[in] width_step : largeur réelle d'une ligne
 * @param[in] value : valeur des pixels n'appartenant pas à l'arrière plan
 * \return 
 * 1 si le point appartient au masque \n
 * 0 sinon
 * \brief La fonction teste l'appartenance d'un point au masque.\n
 * Attention : Un point en dehors de l'image est considéré comme n'appartenant pas au masque!
 */
template <class T> int is_mask(const unsigned int & x,
				               const unsigned int & y,
							   const T * pixels,
							   unsigned int width,
							   unsigned int height,
							   unsigned int width_step,
							   const T & value = 0);

#define IS_MASK(type) template int is_mask(const unsigned int & x,\
										  const unsigned int & y,\
										  const type * pixels,\
										  unsigned int width,\
										  unsigned int height,\
										  unsigned int width_step,\
										  const type & value = 0);
#endif
