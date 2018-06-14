/**@file convex_boundary.hpp
 * 
 */
#ifndef _CONVEXE_BOUNDARY_HPP_
	#define _CONVEXE_BOUNDARY_HPP_
	#include "m_c_convexe.hpp"

	/**@class convex_boundary
	 * @brief
	 * Classe qui permet de calculer l'enveloppe convexe d'un ensemble de points
	 */
	template <class type> class convex_boundary
	{
		public:
			/**@fn convex_boundary(unsigned int max_nb_points);
			 * @param max_nb_points : nombre maximal de points
			 */
			convex_boundary(unsigned int max_nb_points);
			
			/**@fn void setup(unsigned int max_nb_points);
			 * @param max_nb_points : nombre maximal de points
			 */
			virtual void setup(unsigned int max_nb_points);
		
			/**@fn int  compute( const type * points,
								 unsigned int nb_points);
			 * @brief
			 * Cette méthode calcule l'enveloppe convexe de l'ensemble de points
			 */
			int compute( const type * points,
						 unsigned int nb_points,
						 const unsigned char * mask );
			/**@fn double surface() const;
			 * @returen
			 * Surface
			 */
			double surface() const;
			
			/**@fn const unsigned int * labels() const;
			 * Indices des points de l'enveloppe convexe
			 */
			const unsigned int * labels() const;
			/**@fn unsigned int nb_points() const;
			 * Nombre de points de l'enveloppe convexe
			 */
			unsigned int nb_points() const;
			/**@fn unsigned int max_nb_points() const;
			 * Nombre maximum de points supportés
			 */
			unsigned int max_nb_points() const;
			/**@fn ~convex_boundary()
			 * @brief
			 * Destructeru
			 */
			~convex_boundary();
		protected:
			/**@fn void ini_string();
			 * @brief
			 * initialise la chaine.			 
			 **/
			void ini_string();
			
			/**@fn void sort_string();
			 * @brief
			 * trie la chaine.
			 */
			void sort_string();
			
			/**@fn void compute_boundary();
			 * 
			 */
			void compute_boundary();
			
			/**@fn void check_point(m_c_convexe<type> * p);
			 * @brief
			 * Verif d'un point
			 */
			void check_point(m_c_convexe<type> * p);
			/**@fn void compute_label_array();
			 * @brief
			 * Cette méthode calcule le tableau d'indice.
			 */
			void compute_label_array();
			
			/**@fn void free();
			 * 
			 */
			void free();
			/**@fn void initialize();
			 * 
			 */
			void initialize();
			//Paramètres internes
			m_c_convexe<type> * string;
			unsigned int _max_nb_points;
			unsigned int _nb_points;
			unsigned int * _labels;
			double _surface;
			//Autres
			const type * points_;
			unsigned int nb_points_;
			unsigned int nb_points_r;
			const unsigned char * mask_;
	};
	
#endif
