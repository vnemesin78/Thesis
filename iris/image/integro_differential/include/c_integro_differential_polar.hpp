/**@fn
 * @brief
 * Ce fichier contient la classe qui permet d'utiliser l'op. intégro diff. sur une transformée polaire.
 */
#ifndef _C_INTEGRO_DIFF_POLAR_HPP_
	#define _C_INTEGRO_DIFF_POLAR_HPP_
	#include <opencv/cv.h>
	#include "../../interpol_2d/lib_inter_2d.hpp"
	#include "../../polar/lib_polar.hpp"
	#include "c_integro_differential_par.hpp"
	

	using namespace std;
	/**@struct
	 * @brief
	 * Paramètres du courbe polaire
	 * 
	 */
	struct polar_parameters
	{
		polar_parameters ( void );
		polar_parameters ( unsigned int nb_parameters );
		void setup ( unsigned int nb_params );
		polar_parameters & operator=( const polar_parameters & params );
		~polar_parameters ();
		
		double _x;
		double _y;
		double * _params;
		unsigned int _nb_params;
	};
	
	
	/**@class 
	 * @brief
	 * Cette classe permet de gérer un op. intégro différentiel à partir d'une transformée polaire d'une image.
	 * 
	 */
	class c_integro_differential_polar
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur.
			 */
			c_integro_differential_polar( void );
		
			/**@fn
			 * @param nb_directions : nombre de directions testés
			 * @param nb_samples : nombre d'échantillons
			 * @brief
			 * Constructeur.
			 * 
			 */
			c_integro_differential_polar ( unsigned int nb_directions,
											   unsigned int nb_samples );
		
			/**@fn
			 * @param nb_directions : nombre de directions testés
			 * @param nb_samples : nombre d'échantillons
			 * @brief
			 * Setup
			 * 
			 */
			int setup ( 	unsigned int nb_directions,
							unsigned int nb_samples );
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~c_integro_differential_polar( void );
	
	
			/**@fn
			 * @param[out] computed_params : paramètre optimaux
			 * @param image : image
			 * @param mask : masque
			 * @param min_values : valeurs min pour les paramètres
			 * @param max_values : valeurs max pour les paramètres
			 * @param nbs_values : nombre de valeurs testées pour chaque paramètres
			 * @param radius_function : function pour le rayon r(t) = ...
			 * @param radius_derivate : dérivée pour le rayon
			 * @param order : ordre
			 * @brief
			 * Trouve la courbe qui fitte le mieux par un op. intégro différentiel.
			 * @return
			 * - nan en cas de paramètre invalide
			 * - valeur du gradient trouvée
			 */
			double compute ( 	polar_parameters & computed_params,
								const IplImage * image,
								const IplImage * mask,
								const polar_parameters & min_values,
								const polar_parameters & max_values,
								unsigned int * nbs_values,
								curve_function radius_function,
								curve_function radius_derivate,
								order_function order,
								unsigned int kernel_size,
								const double & r_min,
								const double & r_max );
						
		protected:
			
			/**@fn
			 * @brief
			 * Lib. mémoire
			 */
			void free();
			
			/**@fn
			 * @brief
			 * Ini. mémoire
			 * 
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Alloc. mémoire.
			 * 
			 */
			void alloc();
		//Gradients calculés
			IplImage * _h_grad,
					 * _v_grad;
			IplImage * img_32f;
		//Erreurs
			ostream * err_stream;
		protected:
		//Nombre de directions
			unsigned int _nb_directions;
		public:
			inline unsigned int nb_directions() const
			{
				return _nb_directions;
			}
		protected:
		//Nombre d'échantillons
			unsigned int _nb_samples;
		public:
			inline unsigned int nb_samples() const
			{
				return _nb_samples;
			}			
			
		//Objets
		protected:
			c_polar * polar_obj;
			c_interpol_2d * inter_2d;
			
		protected:
			const polar_parameters * _min_params,
									  * _max_params;
			const unsigned int * _nbs_values;
			unsigned int _nb_params;
			
			//Functions
			curve_function _radius_function,
							 _radius_derivate;
			//Ordre
			order_function _order;	 
			
			double _r_min,
					_r_max;
					
			double _r_step,
					_theta_step;
					
					
					
		protected:	

			polar_parameters * _params;
			polar_parameters _step;
			
		protected:
			/**@fn
			 * @brief
			 * Calcul l'intégrale.
			 * 
			 */
			double compute_integral ( const polar_parameters & params );


			/**@fn
			 * @param i : degré de profondeur de la fonction récursive
			 * @brief
			 * Trouve le maxima de manière récursive
			 * 
			 **/
			double search_maxima( 	unsigned int i );
	
			
	
	
	};
	
	
	
	
	
#endif
