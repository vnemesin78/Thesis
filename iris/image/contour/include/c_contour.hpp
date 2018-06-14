/**\file c_contour.hpp
 * \author Valérian Némesin
 * \brief Ce fichier contient le prototype de la classe qui permet de calculer le contour.
 */
#ifndef C_CONTOUR_HPP
#define C_CONTOUR_HPP


/**\class c_contour
 * \brief Cette classe permet de calculer un contour à partir d'un masque binaire. \n
 * Cette classe préalloue la mémoire pour pouvoir stocker un contour de taille maximale. \n
 **/
class c_contour
{
	public:
		/**@fn c_contour :: c_contour(int max_nb_points)
		 * @param max_nb_points : nombre de points maximum du contour
		 * @brief
		 * Constructeur
		 */
		c_contour(unsigned int max_nb_points);
		
		/**@fn void c_contour :: setup(int max_nb_points);
		 * @param max_nb_points : nombre de points maximum du contour
		 * @brief
		 * Setup
		 **/
		void setup(unsigned int max_nb_points);
	
		/**@fn template <class T> void c_contour :: compute_4c(const T * pixels,
												  unsigned int width,
												  unsigned int height,
												  unsigned int width_step,
												  const T & value = 0);
		 * @param pixels : pixels
		 * @param width : largeur
		 * @param height : hauteur
		 * @param width : taille réelle de la ligne
		 * @param value : valeur des pixels du masque
		 * @brief
		 * Cette méthode calcule le contour 4 connexe.
		 */
		template <class T> void compute_4c(const T * pixels,
										   unsigned int width,
										   unsigned int height,
										   unsigned int width_step,
										   const T & value = 0);
	
		/**@fn template <class T> void c_contour :: compute_4c_nb(const T * pixels,
												  unsigned int width,
												  unsigned int height,
												  unsigned int width_step,
												  const T & value = 0);
		 * @param pixels : pixels
		 * @param width : largeur
		 * @param height : hauteur
		 * @param width : taille réelle de la ligne
		 * @param value : valeur des pixels du masque
		 * @brief
		 * Cette méthode calcule le contour 4 connexe.
		 */
		template <class T> void compute_4c_nb(const T * pixels,
										      unsigned int width,
										      unsigned int height,
										      unsigned int width_step,
										      const T & value = 0);
	
	
		/**@fn template <class T> void c_contour :: compute_8c(const T * pixels,
												  unsigned int width,
												  unsigned int height,
												  unsigned int width_step,
												  const T & value = 0);
		 * @param pixels : pixels
		 * @param width : largeur
		 * @param height : hauteur
		 * @param width : taille réelle de la ligne
		 * @param value : valeur des pixels du masque
		 * @brief
		 * Cette méthode calcule le contour 8 connexe.
		 */
		template <class T> void compute_8c(const T * pixels,
										   unsigned int width,
										   unsigned int height,
										   unsigned int width_step,
										   const T & value = 0);

		/**@fn int c_contour :: digitalize(unsigned int nb_points);
		 * @param nb_points : nombre de points du contour
		 * @brief Cette méthode rééchantillonne le contour.
		 **/
		int digitalize(unsigned int nb_points);
		
		/**@fn int c_contour :: perimeter(double & c1,
											  double & c2);
		 */
		int perimeter(double & c1,
					  double & c2);
		
		
		//Acceseurs
		/**@fn inline unsigned int c_contour :: nb_points_contour() const
		 * @return nombre de points du contour.
		 */
		inline unsigned int nb_points_contour() const
		{
			return _nb_points_contour;
		}
		/**@fn inline unsigned int c_contour :: max_nb_points_contour() const
		 * @return nombre maximal de points du contour
		 */
		inline unsigned int max_nb_points_contour() const
		{
			return _max_nb_points_contour;
		}
		
		/**@fn inline unsigned int c_contour :: nb_points_d_contour() const
		 * @return nombre de points du contour rééchantillonné.
		 */
		inline unsigned int nb_points_d_contour() const
		{
			return nb_points_contour_2;
		}
		
		/**@fn inline unsigned int c_contour :: max_nb_points_d_contour() const
		 * @return nombre maximal de points du contour rééchantillonné.
		 */
		inline unsigned int max_nb_points_d_contour() const
		{
			return max_nb_points_contour_2;
		}
		
		/**@fn inline const unsigned int * c_contour :: contour() const
		 * @return contour
		 */
		inline const unsigned int * contour() const
		{
			return _contour;
		}
		/**@fn inline const unsigned int * c_contour :: d_contour() const
		 * @return contour rééchantillonné
		 */
		inline const double * d_contour() const
		{
			return contour_2;
		}
		
		/**@fn  c_contour :: ~c_contour();
		 * @brief
		 * Destructeur
		 **/
		~c_contour();
	protected:
		/**@fn void c_contour :: initialize();
		 * @brief
		 * Mets tous les attributs de l'objet à 0.
		 */
		void initialize();
		/**@fn void c_contour :: free();
		 * @brief
		 * libère la mémoire occupée par les attributs de l'objet.
		 */
		void free();
	
		unsigned int _nb_points_contour;
		unsigned int _max_nb_points_contour;
		unsigned int nb_points_contour_2;
		unsigned int max_nb_points_contour_2;
		unsigned int * _contour;
		double * contour_2;
};
#define COMPUTE_4C_NB(type) template void c_contour :: compute_4c_nb(const type * pixels,\
																	 unsigned int width,\
																	 unsigned int height,\
																	 unsigned int width_step,\
																	 const type & value = 0); 
#define COMPUTE_4C(type) template void c_contour :: compute_4c(const type * pixels,\
															   unsigned int width,\
															   unsigned int height,\
															   unsigned int width_step,\
															   const type & value = 0); 
															   
#define COMPUTE_8C(type) template void c_contour :: compute_8c(const type * pixels,\
															   unsigned int width,\
															   unsigned int height,\
															   unsigned int width_step,\
															   const type & value = 0); 
#endif
