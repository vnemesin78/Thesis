#ifndef _DISCRETE_CONVEXE_BOUNDARY_HPP_
	#define _DISCRETE_CONVEXE_BOUNDARY_HPP_
	#include "convexe_boundary.hpp"
	/**@class template<class type> class discrete_convex_boundary : public convex_boundary<type>
	 * 
	 **/
	template <class type> class discrete_convex_boundary : public convex_boundary<type>
	{
		public:
			/**@fn discrete_convex_boundary(unsigned int max_nb_points);
			 * @param max_nb_points : nombre max de points
			 */
			discrete_convex_boundary(unsigned int max_nb_points);
			
			/**@fn void setup(unsigned int max_nb_points);
			 * @param max_nb_points : nombre max de points
			 * 
			 **/
			virtual void setup(unsigned int max_nb_points);
			
			/**@fn int compute_segments(type dy);
			 * @brief
			 * Cette méthode calcule les points du contour convexe entre les points de l'enveloppe convexe.
			 */
			int compute_segments(type dy = 1);
			
			/**@fn ~discrete_convex_boundary();
			 * 
			 **/
			~discrete_convex_boundary();
			
			const unsigned int * get_left_contour() const;
			
			const unsigned int * get_right_contour() const;	
			
			unsigned int get_y_min() const;
			
			unsigned int get_y_max() const;
			
		protected:
			/**@fn void free();
			 * 
			 */
			void free();
			/**@fn void initialize();
			 * 
			 */
			void initialize();
			
			/**@fn void compute_extrema();
			 * 
			 */
			void compute_extrema();
			unsigned int * left_contour;
			unsigned int * right_contour;
			unsigned int y_min,
						 y_max;
	};
	
	
	
	
#endif
