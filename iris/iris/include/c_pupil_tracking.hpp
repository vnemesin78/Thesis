/**@file c_pupil_tracking.hpp"
 * 
 **/
#ifndef _C_PUPIL_TRACKING_HPP_
	#define _C_PUPIL_TRACKING_HPP_
	#include "c_tracking.hpp"
	/**@class c_pupil_tracking
	 * @brief
	 * Cette classe permet de gérer le suivi de la pupille dans l'image
	 * 
	 **/
	class c_pupil_tracking : public c_tracking
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 **/
			c_pupil_tracking  ();
			
			/**@fn
			 * @param image_width : largeur de l'image
			 * @param image_height : hauteur de l'image
			 * @param t_0 : espérance de l'état initial ( si vous ne connaissez pas, mettre (width / 2, height / 2, r0 ) )
			 * @param sqrt_q_0 : racine de la matrice de covariance de l'état initial 
			 * @param f : Matrice de transtion 
			 * @param sqrt_q : racine de la matrice de covariance du bruit
			 * @brief size_x : dimension de X (==3)
			 * @brief
			 * Constructeur.
			 * @exception
			 * std::bad_argument s'il y a une erreur sur les arguments...
			 **/
			c_pupil_tracking (	unsigned int image_width,
									unsigned int image_height,
									const gsl_vector * t_0,
									const gsl_matrix * sqrt_q_0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									const double & err_sigma,
									ostream * _err_stream = NULL );
			/**@fn
			 * @param image_width : largeur de l'image
			 * @param image_height : hauteur de l'image
			 * @param t_0 : espérance de l'état initial ( si vous ne connaissez pas, mettre (width / 2, height / 2, r0 ) )
			 * @param sqrt_q_0 : racine de la matrice de covariance de l'état initial 
			 * @param f : Matrice de transtion 
			 * @param sqrt_q : racine de la matrice de covariance du bruit
			 * @brief size_x : dimension de X (==3)
			 * @brief
			 * Setup
			 * @return
			 * - 1 s'il y a un mauvais argument
			 **/
			virtual int setup (	unsigned int image_width,
									unsigned int image_height,
									const gsl_vector * t_0,
									const gsl_matrix * sqrt_q_0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									const double & err_sigma,
									ostream * _err_stream = NULL );
				
			int default_setup ( unsigned int width, 
								unsigned int height,
								ostream * _err_stream = NULL );	
			/**@fn
			 * @param image_width : largeur de l'image
			 * @param image_height : hauteur de l'image
			 * @param params : variables en mémoire
			 * @param n_space : Espace de nom des différentes variables
			 * @param t_0_name : nom de la variable contenant t_0
			 * @param sqrt_q_0_nam : nom de la variable contenant sqrt_Q_0
			 * @param f_name : nom de la variable contenant F
			 * @param sqrt_q_name : nom de la variable contenant sqrt_Q
			 * @param n_name : nom de la variable contenant n_sigma
			 * @return
			 * - 0 si les paramètres sont valides
			 * - 1 si des données sont manquantes
			 * - 2 si les paramètres sont invalides
			 * @brief
			 * Setup à partir d'un fichier
			 */
			virtual int setup (	api_parameters & params,
									unsigned int image_width,
									unsigned int image_height,
									ostream * _err_stream = NULL,								
									const char * n_space = "tracking",
									const char * t_0_name = "t_0",
									const char * sqrt_q_0_name = "sqrt_Q_0",
									const char * f_name = "F",
									const char * sqrt_q_name = "sqrt_Q",
									const char * n_name = "n_sigma" );
								
			/**@fn
			 * @param current_x : abscisse de la pupille
			 * @param current_y : ordonnée de la pupille
			 * @param current_r : rayon de la pupille
			 * @param stride : distance mémoire entre 2 composantes de X
			 * @brief
			 * Prédit la futur position de la pupille.
			 */
			virtual void predict_next_position( 	double current_x,
													double current_y,
													double current_r,
													unsigned int frame_id );
				
			/**@fn 
			 * @brief
			 * Reset le tracking ( Utile quand la pupille est perdue ).
			 * 
			 **/
			virtual void reset_tracking ( );
			virtual void reset( );	
			//Accesseur
			/**@fn
			 * @return
			 * abscisse du premier point du roi.
			 * 
			 **/
			inline unsigned int x_roi() const
			{
				return _x_roi;
			}
				
			/**@fn
			 * @return
			 * ordonnée du premier point du roi.
			 * 
			 **/
			inline unsigned int y_roi() const
			{
				return _y_roi;
			}
				
			/**@fn
			 * @return
			 * Largeur du roi
			 * 
			 **/
			inline unsigned int width_roi() const
			{
				return _width_roi;
			}
				
			/**@fn
			 * @return
			 * Hauteur du roi.
			 * 
			 **/
			inline unsigned int height_roi() const
			{
				return _height_roi;
			}
		
		
			/**@fn
			 * @return
			 * Largeur de l'image
			 * 
			 **/
			inline unsigned int image_width() const
			{
				return _image_width;
			}
			
			/**@fn
			 * @return
			 * Hauteur de l'image
			 * 
			 **/
			inline unsigned int image_height() const
			{
				return _image_height;
			}
			
		
		
		protected:
			/**@fn
			 * @brief
			 * Ini. classe fille.
			 * 
			 **/
			void initialize();
			
			unsigned int 	_image_width,
							_image_height;
			
			unsigned int _x_roi,
						 _y_roi,
						 _width_roi,
						 _height_roi;
			

	};
	
	
#endif
