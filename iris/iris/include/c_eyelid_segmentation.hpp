
#ifndef _C_EYELID_SEGMENTATION_HPP_
	#define  _C_EYELID_SEGMENTATION_HPP_
	#include "lib_image.hpp"
	#include <iostream>
	#include "lib_api.hpp"
	using namespace std;
	/**@class 
	 * @brief
	 * Classe pour la segmentation des paupières par l'approximation par une parabole
	 * 
	 */
	class c_eyelid_segmentation
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 **/
			c_eyelid_segmentation ( void );
		
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_eyelid_segmentation ( unsigned int width,
									unsigned int heigth,
									unsigned int nb_t,
									unsigned int nb_a,
									unsigned int nb_c,
									unsigned int nb_theta,
									double theta_min,
									double theta_max,
									unsigned int kernel_size,
									ostream * stream = NULL );
		
			/**@fn
			 * @brief
			 * Setup.
			 * 
			 */
			int setup ( unsigned int width,
						unsigned int heigth,
						unsigned int nb_t,
						unsigned int nb_a,
						unsigned int nb_c,
						unsigned int nb_theta,
						double theta_min,
						double theta_max,
						unsigned int kernel_size,
						ostream * stream = NULL );
		
		
			/**@fn
			 * @brief
			 * Setup.
			 * 
			 */
			int setup ( unsigned int width,
						unsigned int heigth,
						api_parameters & params,
						ostream * stream = NULL,
						const char * n_space = "eyelid",
						const char * nb_t_n = "nb_t",
						const char * nb_a_n = "nb_a",
						const char * nb_c_n = "nb_c",
						const char * nb_theta_n = "nb_theta",
						const char * theta_min_n = "theta_min",
						const char * theta_max_n = "theta_max",
						const char * kernel_size_n = "kernel_size" );
						
					
			/**@fn
			 * @brief
			 * Calcule les paraboles pour l'image en question.
			 * 
			 **/	
			int compute( const IplImage * image,
						 const IplImage * mask,
						 const double & x_iris,
						 const double & y_iris,
						 const double & r_iris,
						 const double & x_pupil,
						 const double & y_pupil,
						 const double & a_pupil,
						 const double & b_pupil );
		
		
			/**@fn
			 * @brief
			 * Retourne le masque calculé de la paupière
			 * 
			 **/
			inline const IplImage * mask() const
			{
				return _mask;
			}
		
			inline double a_upper() const
			{
				return _a_upper;
			}
		
			inline double c_upper() const
			{
				return _c_upper;
			}
		
			inline double theta_upper() const
			{
				return _theta_upper;
			}	
		
			inline double a_lower() const
			{
				return _a_lower;
			}
		
			inline double c_lower() const
			{
				return _c_lower;
			}
		
			inline double theta_lower() const
			{
				return _theta_lower;
			}	
		
		
		
		
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_eyelid_segmentation(); 
		protected:
			
			
			void initialize();
			
			void free();
			
			void alloc();
			
			/**@fn
			 * @brief
			 * Mise à jour du masque!
			 * 
			 */
			void update_mask( const double * params, const double & r_iris, bool up );
			
			
			c_integro_differential_2d * op;
			double 	* max_values,
					* min_values;
			unsigned int * nbs_values;
			
			IplImage * _mask,
					 * _mask_0;
			unsigned int _width,
						 _height;
			unsigned int _nb_t, _kernel_size;
			ostream * err_stream;
			
			
			double _a_upper,
				   _c_upper,
				   _theta_upper;
				   
			double _a_lower,
				   _c_lower,
				   _theta_lower;	   
				  
			
			
			
			IplConvKernel * structuring_element_1;
			
	};
	
	
	
	
	
#endif
