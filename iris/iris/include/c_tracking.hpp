/**@file c_tracking.hpp
 * 
 **/
#ifndef _C_TRACKING_2_HPP_
	#define _C_TRACKING_2_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_blas.h>
	
	#include "lib_api.hpp"
	#include "lib_tkalman_nc.hpp"
	#include <sstream>
	#include <iostream>
	using namespace std;
	
	/**@class c_tracking
	 * @brief
	 * Classe de gestion du tracking
	 * 
	 **/
	class c_tracking
	{
		public:
			
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 **/
			c_tracking();
			
			/**@fn
			 * @param t_0 : espérance de l'état initial ( si vous ne connaissez pas, mettre (width / 2, height / 2, r0 ) )
			 * @param sqrt_q_0 : racine de la matrice de covariance de l'état initial 
			 * @param f : Matrice de transtion 
			 * @param sqrt_q : racine de la matrice de covariance du bruit
			 * @param size_x : dimension de X
			 * @brief
			 * Constructeur.
			 * @exception
			 * std::bad_argument s'il y a une erreur sur les arguments...
			 **/
			c_tracking ( 	const gsl_vector * t_0,
							 const gsl_matrix * sqrt_q_0,
							 const gsl_matrix * f,
							 const gsl_matrix * sqrt_q,
							 unsigned int size_x,
							 const double & err_sigma,
							 ostream * _err_stream = NULL );
	
			/**@fn
			 * @param t_0 : espérance de l'état initial ( si vous ne connaissez pas, mettre (width / 2, height / 2, r0 ) )
			 * @param sqrt_q_0 : racine de la matrice de covariance de l'état initial 
			 * @param f : Matrice de transtion 
			 * @param sqrt_q : racine de la matrice de covariance du bruit
			 * @param size_x : dimension de X
			 * @brief
			 * Setup
			 * @return
			 * - 1 s'il y a un mauvais argument
			 **/
			virtual int setup (  	const gsl_vector * t_0,
									const gsl_matrix * sqrt_q_0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									unsigned int size_x,
									const double & err_sigma,
									ostream * _err_stream = NULL );
	
			/**@fn
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
									ostream * _err_stream = NULL,
									const char * n_space = "tracking",
									const char * t_0_name = "t_0",
									const char * sqrt_q_0_name = "sqrt_Q_0",
									const char * f_name = "F",
									const char * sqrt_q_name = "sqrt_Q",
									const char * size_x_name = "size_x",
									const char * n_name = "n_sigma" );
		
			/**@fn
			 * @param error_str : flux d'erreur
			 * @brief
			 * Modifie le flux d'erreur.
			 * 
			 */
			inline void set_error_stream( ostream & error_str = cout )
			{
				err_stream = &error_str;
			}
			
			/**@fn
			 * @param x : observation courante du vecteur X
			 * @param stride : distance mémoire entre 2 composantes de X
			 * @brief
			 * Prédit la futur position de la pupille.
			 */
			virtual void predict_next_position( const double * x,
												unsigned int frame_id,
												unsigned int stride = 1);
			
			/**@fn 
			 * @brief
			 * Reset le tracking ( Utile quand la pupille est perdue ).
			 * 
			 **/
			virtual void reset_tracking ( );
			
			/**@fn 
			 * @brief
			 * Reset le tracking ( Utile quand la pupille est perdue ).
			 * 
			 **/
			virtual void reset( );	
			
			
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_tracking();
			
		//Accesseurs
			/**@fn
			 * @return
			 * Valeurs min du vecteur x
			 * 
			 **/
			inline const gsl_vector * x_min ( ) const
			{
				return _x_min;
			}
			
			/**@fn
			 * @return
			 * Valeurs max du vecteur x
			 * 
			 **/
			inline const gsl_vector * x_max ( ) const
			{
				return _x_max;
			}
			
			/**@fn
			 * @return
			 * Valeur prédite pour x
			 **/
			inline const gsl_vector * x ( ) const
			{
				return x_p;
			}
			
			/**@fn
			 * @return Espérance de l'état initial
			 **/
			inline const gsl_vector * t_0() const
			{
				return _t_0;
			}
		
			/**@fn
			 * @return
			 * Racine de la matrice de covariance de l'état initial
			 */
			inline const gsl_matrix * sqrt_q_0() const
			{
				return _sqrt_q_0;
			}
		
			/**@fn
			 * @return
			 * Racine de la matrice de covariance du bruit
			 */
			inline const gsl_matrix * sqrt_q() const
			{
				return _sqrt_q;
			}
		
			/**@fn
			 * @return
			 * Matrice d'évolution
			 */
			inline const gsl_matrix * f() const
			{
				return _f;
			}
		
			/**@fn
			 * @return
			 * Dimension de x.
			 * 
			 **/
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn
			 * @return
			 * % d'objets en dehors de l'intervalle!
			 * 
			 **/
			inline double err_sigma() const
			{
				return _err_sigma;
			}
			
			/**@fn
			 * @return
			 * - 0 si le dernier objet a été détecté
			 * - 1 sinon.
			 * 
			 **/
			inline bool previous_ok() const
			{
				return _previous_ok;
			}
		protected:
		
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			
			/**@fn
			 * @brief
			 * Ini. objet
			 **/
			void initialize();
			
			//Limite du cube pour la recherche des paramètres
			gsl_vector * _x_min, 
					   * _x_max;
			bool _previous_ok;
			
			//Paramètres
			unsigned int _size_x,
						 _size_y,
						 _size_t;
			gsl_vector * _t_0;
			gsl_matrix * _sqrt_q_0,
					   * _f,
					   * _sqrt_q;
			gsl_matrix f_yt,
					   sqrt_q_yy;
			gsl_vector * x_p,
					   * x_f,
					   * innovation,
					   * t_f;
			gsl_vector * _y,
					   * y;
			gsl_matrix * sqrt_q_f,
					   * sqrt_p_f,
					   * sqrt_p_p,
					   * sqrt_s;
					   
			gsl_vector * t_p,
					   * vect_t,
							vect_t_view_x;
					   
			gsl_matrix * sqrt_q_p,
					   * mat_tt,
							mat_tt_view_xx;
					   
					   
			double _err_sigma;
			tkalman_nc_filtering * filtering;
			tkalman_nc_prediction * prediction;
		
			gsl_vector * vect_one;
		
			ostream * err_stream;
			unsigned int _frame_id;
	};
	
	
#endif
