/**@file tkalman_nc_sums.hpp
 * @author Valérian Némesin
 */
#ifndef _TKALMAN_NC_SUMS_MARKOV_HPP_
	#define _TKALMAN_NC_SUMS_MARKOV_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include "tkalman_nc_sums_base.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <iostream> //A virer
	using namespace std;
	/**@class tkalman_nc_sums_markov
	 * @brief
	 * Cette classe permet de calculer les sommes nécessaires à l'algorithme EM du filtre de Kalman couple non-corrélé.
	 */
	class tkalman_nc_sums_markov : public tkalman_nc_sums_base
	{
		public:
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @param[in] p : nombre de classes
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_sums
			 */
			tkalman_nc_sums_markov( const gsl_matrix * f_xt,
									const gsl_matrix * sqrt_q_xx,
									unsigned int p) throw (exception &);
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @param[in] p : nombre de classes
			 * @brief
			 * Setup de la classe @class tkalman_nc_sums
			 */
			virtual void setup ( const gsl_matrix * f_xt,
								 const gsl_matrix * sqrt_q_xx,
								 unsigned int p) throw (exception &); 
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~tkalman_nc_sums_markov();
		
			/**@fn
			 * @param[out] sqrt_sum_corr_tx_ : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 * @param[in] z_s : $\hat{x}_{n|N}$, espérances des états lissés
			 * @param[in] y : $y_n$, observations
			 * @param[in] sqrt_z_f : $[Z_{n|n}]^{\frac{1}{2}}$, racines des matrices de covariance des états filtrés.
			 * @param[in] sqrt_p_s : $[Z_{n|N}]^{\frac{1}{2}}$, racines des matrices de covariance des états lissés.
			 * @param[in] c_s : $[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T$
			 * @param[in] r : processus à saut
			 * @param[in] r_id : id de la somme
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction calcule $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 */
			void compute_sqrt_sum_corr_tx_(	gsl_matrix * sqrt_sum_corr_tx_,
											const gsl_vector * const * z_s,
											const gsl_vector * const * y,
											const gsl_matrix * const * sqrt_z_f,
											const gsl_matrix * const * sqrt_p_s,
											const gsl_matrix * const * c_s,
											const unsigned int * r,
											unsigned int r_id,
											unsigned int n);
			
			/**@fn
			 * @param[out] sqrt_sum_corr_tx_ : $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, y_{n|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 * @param[in] z_s : $\hat{x}_{n|N}$, espérances des états lissés
			 * @param[in] y : $y_n$, observations
			 * @param[in] sqrt_p_s : $[Z_{n|N}]^{\frac{1}{2}}$, racines des matrices de covariance des états lissés.
			 * @param[in] r : processus à saut
			 * @param[in] r_id : id de la somme
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction calcule $\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 */
			void compute_sqrt_sum_corr_ty (	gsl_matrix * sqrt_sum_corr_ty,
											const gsl_vector * const * z_s,
											const gsl_vector * const * y,
											const gsl_matrix * const * sqrt_p_s,
											const unsigned int * r,
											unsigned int r_id,
											unsigned int n);
		

		protected:
		
			/**@fn
			 * @brief
			 * Cette fonction met tous les attributs de l'objet à 0.
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Cette fonction alloue les différents élements de la classe.
			 * @throw
			 * bad_alloc en cas de problème de mémoire
			 */
			void alloc() throw(exception &);
			
			/**@fn
			 * @brief
			 * Cette fonction désalloue tous les attributs alloués.
			 **/
			void free();
			
			/**@fn
			 * @brief 
			 * Cette fonction génère les différentes vues sur les matrices.
			 */
			void create_views();
		
			//Paramètres
				gsl_matrix * _f_xx;
				gsl_matrix * _f_xy;
			unsigned int _p;
			
		//Acesseurs
		public:
			
			/**@fn inline const gsl_matrix * tkalman_nc_sums :: f_xx() const
			 * @return
			 * F^{x,x}
			 **/
			inline const gsl_matrix * f_xx() const
			{
				return _f_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_sums :: f_xy() const
			 * @return
			 * F^{x,y}
			 **/
			inline const gsl_matrix * f_xy() const
			{
				return _f_xy;
			}
			
			inline unsigned int p() const
			{
				return _p;
			}
	};
	
	
	
	
#endif
