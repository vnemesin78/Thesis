
#ifndef _C_INTEGRO_DIFFERENTIAL_PAR_HPP_
	#define _C_INTEGRO_DIFFERENTIAL_PAR_HPP_
	
	#include <opencv/cv.h>
	#include "../../interpol_2d/lib_inter_2d.hpp"
	#include "../../polar/lib_polar.hpp"
	using namespace std;
	
	/**@typedef
	 * @brief
	 * typedef des fonctions de courbes
	 * 
	 **/
	typedef double (*curve_function) ( 	const double & x,
											const double * params );


	/**@typedef
	 * @brief
	 * typedef des fonctions de courbes
	 * 
	 **/
	typedef bool (*order_function) ( 	const double & x,
										const double & y );
	
	/**
	 * 
	 * 
	 **/
	class c_integro_differential_par
	{
		public:
		
		
			/**@fn
			 * @brief
			 * Constructor
			 * 
			 **/
			c_integro_differential_par ( );
		
			/**@fn
			 * @param width
			 * @param height
			 * @param depth
			 * - IPL_DEPTH_8U
			 * - IPL_DEPTH_32F
			 * - IPL_DEPTH_64F
			 * @param kernel_size : size of gaussian kernel
			 * @brief
			 * Constructor.
			 **/
			c_integro_differential_par ( 	unsigned int width,
											unsigned int height,
											unsigned int kernel_size );
		
			/**@fn
			 * @param width
			 * @param height
			 * @param depth
			 * - IPL_DEPTH_8U
			 * - IPL_DEPTH_32F
			 * - IPL_DEPTH_64F
			 * @param kernel_size : size of gaussian kernel
			 * @return
			 * 1 - Error
			 * @brief
			 * Setup
			 **/
			int setup (	unsigned int width,
						unsigned int height,
						unsigned int kernel_size );
		
		
			/**@fn
			 * @param params : paramètres de la courbe optimale
			 * @param image : image
			 * @param mask : masque
			 * @param min_values : valeurs min des params.
			 * @param max_values : valeurs max des params.
			 * @param nbs_values : nombre de valeurs testées pour chaques paramètres
			 * @param nb_params : nombre de paramètres
			 * @param curve : function de la courbe
			 * @param order : function d'ordre
			 * @return 
			 * - nan
			 * - Valeur de l'intégrale
			 * @brief
			 * Calcule la meilleure courbe pour l'image
			 *
			 **/
			double compute(	double * params,
							const IplImage * image,
							const IplImage * mask,
							const double * min_values,
							const double * max_values,
							unsigned int * nbs_values,
							unsigned int nb_params,
							curve_function curve,
							curve_function derivate,
							order_function order);
							
			/**@fn
			 * @brief
			 * Destructor
			 * 
			 **/
			~c_integro_differential_par ( );
			
		protected:
			/**@fn
			 * @brief
			 * Ini.
			 */
			void initialize();
			
			
			/**@fn 
			 * @brief
			 * Free.
			 * 
			 **/
			void free();
			
			/**@fn
			 * @brief
			 * Calcule l'intégrale d'une ligne
			 * 
			 **/
			double compute_integral( 	const IplImage * mask, 
										const double * curr_params, 
										curve_function curve,
										curve_function derivate );
			/**@fn
			 * @brief
			 * Trouve le maxima
			 * 
			 **/
			double search_maxima( 	double * params,
									double * curr_params,
									const double * min_values,
									const double * max_values,
									unsigned int * nbs_values,
									unsigned int i,
									unsigned int nb_params,
									const IplImage * mask,
									curve_function curve, 
									curve_function derivate,
									order_function order );
			
			unsigned int 	_width,
							_height,
							_kernel_size;
			
			//Params
			ostream * err_stream;
			
			
			//Image tmp
			c_interpol_2d interpol;
			IplImage * img_32f,
					 * h_grad,
					 * v_grad;

	};

#endif
