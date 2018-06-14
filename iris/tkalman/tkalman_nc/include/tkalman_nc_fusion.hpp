/**@file tkalman_nc_fusion.hpp
 * @author Valérian Némesin
 */

#ifndef _TKALMAN_NC_FUSION_HPP_
	#define _TKALMAN_NC_FUSION_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <iostream> //A virer
	#include <cmath>
	/**@class tkalman_nc_fusion
	 * @brief
	 * Cette classe permet de faire la fusion des sommes issues de différents jeux de données de type PKF.
	 * 
	 **/
	class tkalman_nc_fusion
	{
		public:
			/**@fn 
			 * @param size_x : taille du vecteur x
			 * @param size_y : taille du vecteur y
			 * @param n : nombre de jeux de données à fusionner
			 * @brief
			 * Constructeur
			 **/
			tkalman_nc_fusion ( unsigned int size_x,
								unsigned int size_y );
			/**@fn 
			 * @param size_x : taille du vecteur x
			 * @param size_y : taille du vecteur y
			 * @param n : nombre de jeux de données à fusionner
			 * @brief
			 * Setup
			 **/
			void setup ( unsigned int size_x,
						 unsigned int size_y );
		
			/**@fn 
			 * @param[out] sqrt_sum_corr_tx : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 * @param[out] sqrt_sum_corr_ty : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, y_{n|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 * @param[in] sqrt_sums_tx : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 * @param[in]  sqrt_sums_corr_ty : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, y_{n|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 * @param[in] n : nombre de jeux de données à fusionner.
			 * @brief
			 * Fonction qui effectue la fusion robuste des différentes sommes.
			 **/
			void fusion ( gsl_matrix * sqrt_sum_tx,
						  gsl_matrix * sqrt_sum_ty,
						  const gsl_matrix * const * sqrt_sums_tx,
						  const gsl_matrix * const * sqrt_sums_ty,
						  unsigned int n,
						  unsigned int step = 1 );
			
			/**@fn 
			 * @param t_0 : espérance de l'état initial
			 * @param sqrt_q_0 : racine de la matrice de covariance de l'état initial
			 * @param t_s_0 : espérance des états initiaux lissés des différents signaux
			 * @param sqrt_q_s_0 : racine des matrices de covariance des états initiaux lissés des différents signaux
			 * @param nb_signals : nombre de signaux
			 * @brief
			 * Cette méthode calcule l'espérance et la matrice de covariance de l'état initial.
			 * 
			 **/
			void fusion_0 ( gsl_vector * t_0,
							gsl_matrix * sqrt_q_0,
							const gsl_vector * const * t_s_0,
							const gsl_matrix * const * sqrt_q_s_0,
							unsigned int nb_signals );
							
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~tkalman_nc_fusion();
		
			/**@fn inline unsigned int tkalman_nc_sums :: size_x() const
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int tkalman_nc_sums :: size_y() const
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int tkalman_nc_sums :: size_t() const
			 * @return 
			 * Dim. de t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
		
		protected:
			/**@fn
			 * @brief
			 * Initialisation
			 * 
			 */
			void initialize();
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 */
			void free();
			
			gsl_matrix * mat_2t_x,
					   mat_2t_x_view_00,
					   mat_2t_x_view_10,
					   * mat_2t_y,
					   mat_2t_y_view_00,
					   mat_2t_y_view_10;
			gsl_vector * vect_x,
					   * vect_y;
			
			gsl_matrix * mat_2t_t,
						 mat_2t_t_view_0,
						 mat_2t_t_view_1,
						 mat_2t_t_view_tp1_t;
			gsl_vector mat_2t_t_view_vector;
			gsl_vector * vect_t;
			
			
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
	};
	
	
#endif
