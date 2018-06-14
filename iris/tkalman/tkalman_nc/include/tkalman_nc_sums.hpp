/**@file tkalman_nc_sums.hpp
 * @author Valérian Némesin
 */
#ifndef _TKALMAN_NC_SUMS_HPP_
	#define _TKALMAN_NC_SUMS_HPP_
	#include "tkalman_nc_sums_base.hpp"
	using namespace std;
	/**@class tkalman_nc_sums
	 * @brief
	 * Cette classe permet de calculer les sommes nécessaires à l'algorithme EM du filtre de Kalman couple non-corrélé.
	 */
	class tkalman_nc_sums : public tkalman_nc_sums_base
	{
		public:
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_sums
			 */
			tkalman_nc_sums( const gsl_matrix * f_xt,
							 const gsl_matrix * sqrt_q_xx) throw (exception &);
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @brief
			 * Setup de la classe @class tkalman_nc_sums
			 */
			void setup ( const gsl_matrix * f_xt,
					     const gsl_matrix * sqrt_q_xx) throw (exception &); 
		
			/**@fn void compute_sqrt_sum_corr_tx_(gsl_matrix * sqrt_sum_corr_tx_,
											 const gsl_vector * const * z_s,
											 const gsl_vector * const * y,
											 const gsl_matrix * const * sqrt_z_f,
											 const gsl_matrix * const * sqrt_p_s,
											 const gsl_matrix * const * c_s,
											 unsigned int n);
			 * @param[out] sqrt_sum_corr_tx_ : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 * @param[in] z_s : $\hat{x}_{n|N}$, espérances des états lissés
			 * @param[in] y : $y_n$, observations
			 * @param[in] sqrt_z_f : $[Z_{n|n}]^{\frac{1}{2}}$, racines des matrices de covariance des états filtrés.
			 * @param[in] sqrt_p_s : $[Z_{n|N}]^{\frac{1}{2}}$, racines des matrices de covariance des états lissés.
			 * @param[in] c_s : $[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T$
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction calcule $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 */
			void compute_sqrt_sum_corr_tx_(gsl_matrix * sqrt_sum_corr_tx_,
										   const gsl_vector * const * z_s,
										   const gsl_vector * const * y,
										   const gsl_matrix * const * sqrt_z_f,
										   const gsl_matrix * const * sqrt_p_s,
										   const gsl_matrix * const * c_s,
										   unsigned int n);
			
			/**@fn void compute_sqrt_sum_corr_ty ( gsl_matrix * sqrt_sum_corr_ty,
												   const gsl_vector * const * z_s,
												   const gsl_vector * const * y,
												   const gsl_matrix * const * sqrt_p_s,
												   unsigned int n);
			 * @param[out] sqrt_sum_corr_tx_ : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, y_{n|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 * @param[in] z_s : $\hat{x}_{n|N}$, espérances des états lissés
			 * @param[in] y : $y_n$, observations
			 * @param[in] sqrt_p_s : $[Z_{n|N}]^{\frac{1}{2}}$, racines des matrices de covariance des états lissés.
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction calcule $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 */
			void compute_sqrt_sum_corr_ty ( gsl_matrix * sqrt_sum_corr_ty,
										    const gsl_vector * const * z_s,
										    const gsl_vector * const * y,
											const gsl_matrix * const * sqrt_p_s,
											unsigned int n);
		
		protected:
		
			/**@fn void tkalman_nc_sums :: create_views();
			 * @brief 
			 * Cette fonction génère les différentes vues sur les matrices.
			 */
			void create_views();
		
			//Paramètres
				gsl_matrix _f_xx;
				gsl_matrix _f_xy;
			
		//Accesseurs
		public:
			/**@fn inline const gsl_matrix * tkalman_nc_sums :: f_xx() const
			 * @return
			 * F^{x,x}
			 **/
			inline const gsl_matrix * f_xx() const
			{
				return &_f_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_sums :: f_xy() const
			 * @return
			 * F^{x,y}
			 **/
			inline const gsl_matrix * f_xy() const
			{
				return &_f_xy;
			}
			
						 
			
	};
	
	
#endif
