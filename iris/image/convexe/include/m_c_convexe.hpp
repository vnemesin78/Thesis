/**@file m_c_convexe.hpp
 * @author Valérian Némesin
 */
#ifndef _M_C_CONVEXE_
	#define _M_C_CONVEXE_
	/**@class m_c_convexe
	 * @brief
	 * Cette classe gère un point d'un ensemble convexe
	 */
	template<class type> class m_c_convexe
	{
		public:
			/**@fn m_c_convexe()
			 * @brief
			 * Constructeur par défaut.
			 */
			m_c_convexe();
		
			/**@fn m_c_convexe( const m_c_convexe<type> * first,
			 * 					const type * point,
								unsigned int label);
			 * @param first : premier point du maillon
			 * @param point : point
			 * @param label : n° du point
			 */
			m_c_convexe(const m_c_convexe<type> * first,
						const type * point,
						unsigned int label);
		
			/**@fn void setup( const m_c_convexe<type> * first,
							   const type * point,
							   unsigned int label);
			 * @param point : point
			 * @param label : n° du point
			 */
			void setup( const m_c_convexe<type> * first,
						const type * point,
						unsigned int label);
			
			/**@fn m_c_convexe<type> * next(m_c_convexe<type> * _next = NULL);
			 * @param _next : maillon suivant
			 * @brief
			 * Cette méthode renvoie le maillon suivant ou le maillon suivant défini si _next = NULL
			 */
			m_c_convexe<type> * next(m_c_convexe<type> * _next = 0);
		
			/**@fn m_c_convexe<type> * previous(m_c_convexe<type> * _previous = NULL);
			 * @param _previous : maillon suivant
			 * @brief
			 * Cette méthode renvoie le maillon suivant ou le maillon suivant défini si _previous = NULL
			 */
			m_c_convexe<type> * previous(m_c_convexe<type> * _previous = 0);
			
			/**@fn const type * point() const;
			 * @return point du contour
			 */
			const type * point() const;
			
			/**@fn unsigned int label() const;
			 * @return
			 * n° du point
			 **/
			unsigned int label() const;
			
			/**@fn double cosinus() const
			 * @return 
			 * cosinus de l'angle avec le premier point
			 */
			double cosinus() const;
			
			/**@fn double radius() const
			 * @return 
			 * cosinus de l'angle avec le premier point
			 */
			double radius() const;
			
			/**@fn double scalar_product() const;
			 * @return
			 * Produit scalaire entre le point précédent et ce point.
			 */
			double scalar_product() const;
			
			/**@fn double determinant() const;
			 * @return
			 * |1	x_p	y_p	|
			 * |1	x	y	|
			 * |1	x_n	y_n	|
			 **/
			double determinant() const;
			
			/**@fn void erase();
			 * @brief
			 * Supprime ce point de la chaine.			 
			 * */
			void erase();
			
			/**@fn bool operator<(const m_c_convexe<type> & m) const;
			 * 
			 */
			bool operator<(const m_c_convexe<type> & m) const;
			
		protected:
			const type * _point;
			unsigned int _label;
			m_c_convexe<type> * _next;
			m_c_convexe<type> * _previous;
			const m_c_convexe<type> * _first;
	};
	
#endif
