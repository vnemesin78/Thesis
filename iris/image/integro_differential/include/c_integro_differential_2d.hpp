
#ifndef _C_INTEGRO_DIFFERENTIAL_2D_HPP_
	#define _C_INTEGRO_DIFFERENTIAL_2D_HPP_
	#include <opencv/cv.h>
	#include "../../interpol_2d/lib_inter_2d.hpp"
	#include "c_integro_differential_par.hpp"
	using namespace std;
	
	/**@struct
	 * @brief
	 * Structure pour les courbes 2d
	 * 
	 **/
	struct function_2d
	{
		curve_function x_curve;
		curve_function y_curve;
		
		/**@fn
		 * @brief
		 * Constructeur
		 * @param x_curve : function x: t->x(t)
		 * @param y_curve : function y: t->y(t)
		 * @param t_min : t_min
		 * @param t_max : t_max
		 * 
		 **/
		function_2d ( 	curve_function x_curve = NULL,
						curve_function y_curve = NULL);
		
		/**@fn
		 * @brief
		 * Setup
		 * @param x_curve : function x: t->x(t)
		 * @param y_curve : function y: t->y(t)
		 * @param t_min : t_min
		 * @param t_max : t_max
		 * 
		 **/
		void setup ( 	curve_function x_curve = NULL,
						curve_function y_curve = NULL);
		
		/**@fn
		 * @brief
		 * Renvoie le point de la courbe pour t.
		 * 
		 **/
		void get_value (	double & x,
							double & y,
							const double & t,
							const double * params ) const;
							
		bool operator!() const;
	};
	
	
	
	/**
	 * 
	 * 
	 **/
	class c_integro_differential_2d
	{
		public:
		
		
			/**@fn
			 * @brief
			 * Constructor
			 * 
			 **/
			c_integro_differential_2d ( );
		
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
			c_integro_differential_2d (	unsigned int width,
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
								const double & t_min,
								const double & t_max,
								unsigned int nb_t,
								const function_2d & f,
								const function_2d & df,
								order_function order );
							
			/**@fn
			 * @brief
			 * Destructor
			 * 
			 **/
			~c_integro_differential_2d ( );
			
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
										const double & t_min,
										const double & t_max,
										unsigned int nb_t,
										const function_2d & f,
										const function_2d & df);
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
									const double & t_min,
									const double & t_max,
									unsigned int nb_t,
									const function_2d & f,
									const function_2d & df,
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
