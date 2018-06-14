
#ifndef _TKALMAN_NC_EM_MARKOV_HPP_
	#define _TKALMAN_NC_EM_MARKOV_HPP_
	#include "tkalman_nc_em_base.hpp"
	#include "auxi_function_tools.hpp"
	#include "tkalman_nc_filter_markov.hpp"
	#include <cstring>
	
	class tkalman_nc_em_markov : public tkalman_nc_em_base
	{
		public:
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] f : F, matrice d'évolution
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @param[in] p_max : nombre maximal de signaux
			 * @param[in] x_mask : Masque sur la matrice Fxt
			 * si x_mask(i) == 0, alors Fxt(:, i) == 0
			 * @param[in] y_mask : Masque sur la matrice Fyt
			 * si y_mask(i) == 0, alors Fyt(:, i) == 0
			 * @brief
			 * Constructeur
			**/
			tkalman_nc_em_markov(	const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									unsigned int size_x,
									unsigned int n_max,
									unsigned int p_max,
									unsigned int nb_signal_max,
									const gsl_vector * x_mask = NULL,
									const gsl_vector * y_mask = NULL,
									bool estimate_initial_state = false ) throw (exception &);
		
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] f : F, matrice d'évolution
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @param[in] p_max : nombre maximal de signaux
			 * @param[in] x_mask : Masque sur la matrice Fxt
			 * si x_mask(i) == 0, alors Fxt(:, i) == 0
			 * @param[in] y_mask : Masque sur la matrice Fyt
			 * si y_mask(i) == 0, alors Fyt(:, i) == 0
			 * @brief
			 * Setup
			 **/
			void setup(	const gsl_vector * t0,
						const gsl_matrix * sqrt_q0,
						const gsl_matrix * f,
						const gsl_matrix * sqrt_q,
						unsigned int size_x,
						unsigned int n_max,
						unsigned int p_max,
						unsigned int nb_signal_max,
						const gsl_vector * x_mask = NULL,
						const gsl_vector * y_mask = NULL,
						bool estimate_initial_state = false ) throw (exception &);
			/**@fn
			 * Destructeur de la classe @class tkalman_nc_em.
			 */	
			~tkalman_nc_em_markov();
						
			/**@fn
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 */
			void filter(const gsl_vector * const * observations,
						const unsigned int * r,
					    unsigned int n) throw(exception &);
					    
			/**@fn
			 * @param[in] observations : observations
			 * @param[in] n : nombre d'observations
			 * @brief
			 * Cette fonction effectue le lissage des observations par le filtre de Kalman couple. Le résultat est stocké dans x_f et sqrt_p_f.
			 */
			void smooth(const gsl_vector * const * observations,
						const unsigned int * r,
					    unsigned int n) throw(exception &);
					    
			/**@fn
			 * @param[in] observations : observations des différents signaux
			 * @param[in] n : nombre d'observations par signal
			 * @param[in] p : nombre de signaux à analyser.
			 * @param[in] nb_iter : nombre d'itérations de l'EM.
			 * @brief
			 * Cette méthode permet d'apprendre les paramètres optimaux pour le jeu de données
			 */
			virtual void learn_parameters(const  gsl_vector * const * const * observations,
										  const unsigned int * const * r,
										  const unsigned int * n,
										  unsigned int p,
										  unsigned int nb_iter) throw(exception &);	
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
		
			/**@fn
			 * @brief
			 * Cette fonction crée les attributs
			 */
			void create_object() throw(exception &);
		
		
			    
			//Objets
			tkalman_nc_filter_markov * filter_obj;
			
			//Sommes
			gsl_matrix ** sqrt_sums_tx;
			gsl_matrix ** sqrt_sums_ty;
			
			gsl_matrix ** _f,
						* _f_xt,
						* _f_yt;
			gsl_matrix ** _sqrt_q;
				gsl_matrix * _sqrt_q_xx;
				gsl_matrix * _sqrt_q_yy;
			
			unsigned int * sum_n;
			unsigned int _p_max;
		
		//Accesseurs
		public:
			/**@fn
			 * @return [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 */
			inline const gsl_matrix * const * sqrt_q() const
			{
				return _sqrt_q;
			}
			
			/**@fn
			 * @return F, matrice d'évolution
			 */
			inline const gsl_matrix * const * f() const
			{
				return _f;
			}
			
			/**@fn
			 * @return F^{x,x}
			 */
			inline const gsl_matrix * f_xx() const
			{
				return filter_obj->f_xx();
			}
			
			/**@fn
			 * @return F^{x,y}
			 */
			inline const gsl_matrix * f_xy() const
			{
				return filter_obj->f_xy();
			}
			
			/**@fn
			 * @return F^{x,t}
			 */
			inline const gsl_matrix * f_xt() const
			{
				return filter_obj->f_xt();
			}
			
			/**@fn
			 * @return F^{y,x}
			 */
			inline const gsl_matrix * f_yx() const
			{
				return filter_obj->f_yx();
			}
			
			/**@fn
			 * @return F^{y,y}
			 */
			inline const gsl_matrix * f_yy() const
			{
				return filter_obj->f_yy();
			}
			
			/**@fn
			 * @return F^{y,t}
			 */
			inline const gsl_matrix * f_yt() const
			{
				return filter_obj->f_yt();
			}
			
			/**@fn
			 * @return [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			 */
			inline const gsl_matrix * sqrt_q_xx() const
			{
				return filter_obj->sqrt_q_xx();
			}
			
			/**@fn
			 * @return [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure.
			 */
			inline const gsl_matrix * sqrt_q_yy() const
			{
				return filter_obj->sqrt_q_yy();
			}
			
			/**@fn
			 * @return
			 * \hat{x}_{n|n - 1}, espérance des états prédits.
			 */
			inline const gsl_vector * const * x_p() const
			{
				return filter_obj->x_p();
			}

			/**@fn
			 * @return
			 * \hat{x}_{n|n}, espérance des états filtrés.
			 */
			inline const gsl_vector * const * x_f() const
			{
				return filter_obj->x_f();
			}

			/**@fn
			 * @return
			 * \hat{x}_{n|N}, espérance des états lissés.
			 */
			inline const gsl_vector * const * x_s() const
			{
				return filter_obj->x_s();
			}

			/**@fn
			 * @return
			 * \hat{t}_{n|N}, Prédiction à l'ordre p
			 */
			inline const gsl_vector * const * t_p() const
			{
				return filter_obj->t_p();
			}

			/**@fn
			 * @return
			 * [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 
			 **/
			inline const gsl_matrix * const * sqrt_p_p() const
			{
				return filter_obj->sqrt_p_p();
			}
			
			/**@fn
			 * @return
			 * [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré 
			 **/
			inline const gsl_matrix * const * sqrt_p_f() const
			{
				return filter_obj->sqrt_p_f();
			}
			
			/**@fn
			 * @return
			 * [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé 
			 **/
			inline const gsl_matrix * const * sqrt_p_s() const
			{
				return filter_obj->sqrt_p_s();
			}
			
			/**@fn
			 * @return
			 * [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
			 */
			inline const gsl_matrix * const * c_s() const
			{
				return filter_obj->c_s();
			}
			
			/**@fn
			 * @return
			 * [Q_{n|n - p}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit à l'ordre p
			 **/
			inline const gsl_matrix * const * sqrt_q_p() const
			{
				return filter_obj->sqrt_q_p();
			}
			
			/**@fn
			 * @return 
			 * Nombre d'états supportés - 1.
			 */
			inline unsigned int n_max() const
			{
				return filter_obj->n_max();
			}
					    
	};
	
#endif
