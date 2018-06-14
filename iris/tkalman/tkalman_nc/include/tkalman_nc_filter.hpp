/**@file tkalman_nc_filter.hpp
 * @author Valérian Némesin
**/
#ifndef _TKALMAN_NC_FILTER_HPP_
	#define _TKALMAN_NC_FILTER_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include "tkalman_nc_prediction.hpp"
	#include "tkalman_nc_filtering.hpp"
	#include "tkalman_nc_smoothing.hpp"
	#include "tkalman_nc_filter_base.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	using namespace std;
	/**@class tkalman_nc_filter
	 */
	 class tkalman_nc_filter : public tkalman_nc_filter_base
	 {
		public:
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] f : F, matrice d'évolution
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_filter
			 */
			tkalman_nc_filter(const gsl_vector * t0,
							  const gsl_matrix * sqrt_q0,
							  const gsl_matrix * f,
							  const gsl_matrix * sqrt_q,
							  unsigned int size_x,
							  unsigned int n_max) throw(exception &);
		
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] f : F, matrice d'évolution
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @brief
			 * Setup de la classe @class tkalman_nc_filter
			 */ 
			virtual void setup(const gsl_vector * t0,
							   const gsl_matrix * sqrt_q0,
							   const gsl_matrix * f,
							   const gsl_matrix * sqrt_q,
							   unsigned int size_x,
							   unsigned int n_max) throw(exception &);
			
			/**@fn
			 * @brief
			 * Destructeur de la classe @class tkalman_nc_filter
			 */
			~tkalman_nc_filter();
		
			/**@fn
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 */
			void filter(const gsl_vector * const * observations,
					    unsigned int n) throw(exception &);
						
			/**@fn
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 */
			void smooth(const gsl_vector * const * observations,
					    unsigned int n) throw(exception &);
			
			/**@fn
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 */			
			void predict(const gsl_vector * const * observations,
						 unsigned int n,
						 unsigned int p) throw(exception &);

			/**@fn 
			 * @param sq_sum_corr_tx : \f$\left[ \sum_{n = 0}^N corr((t_{n|N}^T, x_{n + 1|N}^T)^T) \right]^{\frac{1}{2}}\f$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, x_{n + 1|N}^T)^T.
			 * @param sq_sum_corr_ty : sqrt_sum_corr_tx_ : \f$\left[ \sum_{n = 0}^N corr((t_{n|N}^T, y_{n|N}^T)^T) \right]^{\frac{1}{2}}\f$, racine de la somme des corrélations des vecteurs (t_{n|N}^T, y_{n|N}^T)^T.
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 * 
			 */
			void compute_sqrt_sums ( gsl_matrix * sq_sum_corr_tx,
									 gsl_matrix * sq_sum_corr_ty,
									 const gsl_vector * const * observations,
									 unsigned int n ) throw(exception &);

		//Methodes internes
		protected:
			/**@fn void tkalman_nc_filter :: initialize();
			 * @brief
			 * Cette fonction met tous les attributs de l'objet à 0.
			 */
			void initialize();
			
			/**@fn void tkalman_nc_filter :: free();
			 * @brief
			 * Cette fonction désalloue tous les attributs alloués.
			 **/
			void free();
			
			/**@fn void tkalman_nc_filter:: create_views();
			 * @brief 
			 * Cette fonction génère les différentes vues sur les matrices.
			 */
			void create_views();
		
			/**@fn void tkalman_nc_filter :: create_object();
			 * @brief
			 * Cette fonction crée les attributs
			 */
			void create_object() throw(exception &);
		
		
		
		//Attributs internes
		protected:
			const gsl_matrix * _f;
				gsl_matrix _f_xx,
						   _f_xy,
						   _f_xt,
						   _f_yx,
						   _f_yy,
						   _f_yt;
			const gsl_matrix * _sqrt_q;
				gsl_matrix _sqrt_q_xx,
						   _sqrt_q_yy;
			tkalman_nc_sums * sums;
		//Accesseurs
		public:
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: sqrt_q() const
			 * @return [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 */
			inline const gsl_matrix * sqrt_q() const
			{
				return _sqrt_q;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f() const
			 * @return F, matrice d'évolution
			 */
			inline const gsl_matrix * f() const
			{
				return _f;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_xx()  const
			 * @return F^{x,x}
			 */
			inline const gsl_matrix * f_xx() const
			{
				return &_f_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_xy()  const
			 * @return F^{x,y}
			 */
			inline const gsl_matrix * f_xy() const
			{
				return &_f_xy;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_xt()  const
			 * @return F^{x,t}
			 */
			inline const gsl_matrix * f_xt() const
			{
				return &_f_xt;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_yx()  const
			 * @return F^{y,x}
			 */
			inline const gsl_matrix * f_yx() const
			{
				return &_f_yx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_yy()  const
			 * @return F^{y,y}
			 */
			inline const gsl_matrix * f_yy() const
			{
				return &_f_yy;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: f_yt()  const
			 * @return F^{y,t}
			 */
			inline const gsl_matrix * f_yt() const
			{
				return &_f_yt;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: sqrt_q_xx()  const
			 * @return [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			 */
			inline const gsl_matrix * sqrt_q_xx() const
			{
				return &_sqrt_q_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filter :: sqrt_q_yy()  const
			 * @return [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure.
			 */
			inline const gsl_matrix * sqrt_q_yy() const
			{
				return &_sqrt_q_yy;
			}

	 };
	
	
#endif
