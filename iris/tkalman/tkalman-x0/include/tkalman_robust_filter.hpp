#ifndef TKALMAN_ROBUST_FILTER_HPP_INCLUDED
    #define TKALMAN_ROBUST_FILTER_HPP_INCLUDED
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "tkalman_base.hpp"
	#include "tkalman_constants.hpp"
	#include "tkalman_filtering.hpp"
	#include "tkalman_smoothing.hpp"
	#include "gsl_triangle_matrix.hpp"
    class tkalman_robust_filter : public tkalman_base
    {
        public:
			/**@fn  tkalman_robust_filter :: tkalman_robust_filter(const gsl_vector * x0,
																	 const gsl_matrix * sqrt_p0,
																	 const gsl_matrix * f,
																	 const gsl_matrix * sqrt_q,
																	 unsigned int n = 0);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] sqrt_p0 : racine de la matrice de covariance de l'état initial
			 * @param[in] f : Matrice d'évolution
			 * @param[in] sqrt_q : racine de la matrice de covariance
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @brief
			 constructeur de l'objet.
			 */
			tkalman_robust_filter(const gsl_vector * x0,
                                  const gsl_matrix * sqrt_p0,
                                  const gsl_matrix * f,
                                  const gsl_matrix * sqrt_q,
                                  unsigned int n = 0);


            /**@fn virtual int tkalman_robust_filter :: setup(const gsl_vector * x0,
															   const gsl_matrix * p0,
															   const gsl_matrix * f,
															   const gsl_matrix * q,
															   unsigned int n = 0);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] sqrt_p0 : racine de la matrice de covariance de l'état initial
			 * @param[in] f : Matrice d'évolution
			 * @param[in] sqrt_q : racine de la matrice de covariance
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @return
			 * - 0 si l'objet est valide
			 * - 1 en cas de problème
			 * @brief
			 * Cette méthode initialise l'objet.
			 */
			virtual int setup(const gsl_vector * x0,
							  const gsl_matrix * sqrt_p0,
							  const gsl_matrix * f,
							  const gsl_matrix * sqrt_q,
							  unsigned int n = 0);

			/**@fn virtual void tkalman_robust_filter :: filter(const gsl_vector * const * observations,
													   unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual bool tkalman_robust_filter :: operator!() const;
			 * @return
			 * - 0 si l'objet est bien initialisé.
			 * - 1 sinon
			 * @brief
			 * Cette fonction controle les allocations mémoires des différents attributs de l'objet.
			 */
			virtual bool operator!() const;


			/**@fn virtual void tkalman_robust_filter :: smooth(const gsl_vector * const * observations,
																unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le lissage des données par le filtre de Kalman Triple.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn tkalman_robust_filter :: ~tkalman_robust_filter();
			 * @brief
			 destructeur de l'objet.
			 */
			virtual ~tkalman_robust_filter();

			/**@fn void tkalman_robust_filter :: get_p0(gsl_vector * p0) const;
			 * @param p0 : Matrice de covariance de l'état initial (Préallouée)
			 * @param mat_xx : matrice temporaire de taille (x.x) (Préallouée)
			 * @brief
			 Cette fonction calcule la matrice de covariance de l'état initial. (P^{x,x}\;P_0\;[P^{x,x}]^{-1})
			 */
			virtual void get_p0(gsl_matrix * p0,
								gsl_matrix * mat_xx) const;


			/**@fn void  tkalman_robust_filter :: get_q(gsl_matrix * f,
											  gsl_matrix * mat_tt) const;
			 * @param q : Matrice de covariance du bruit (Préallouée)
			 * @param mat_tt : matrice de taille (t.t) (Préallouée)
			 * @brief
			 Cette fonction calcule la matrice de covariance. ((P\;Q\;P^{T})
			 */
			virtual void get_q(gsl_matrix * q,
					   		   gsl_matrix * mat_tt) const;







            /**@fn virtual double tkalman_robust_filter :: log_likelihood(const gsl_vector * const * observations,
											   gsl_matrix * mat_yy_1,
                                               gsl_vector * vect_y,
                                               gsl_permutation * perm_y,
                                               gsl_matrix * mat_yy_2) const
             * @param mat_yy_1 : matrice temporaire de taille (y.y) (préallouée)
             * @param vect_y : vecteur temporaire de taille (y) (préalloué)
             * @param perm_y : permutation de taille (y) préallouée
             * @param mat_yy_2 : NULL (arg. non utilisé)
             * @return
             Valeur du log-vraisemblance.
             */
            double log_likelihood(gsl_matrix * mat_yy_1,
                                  gsl_vector * vect_y_1,
                                  gsl_permutation * perm_y) const;


		protected :

			/**@fn void tkalman_robust_filter :: filter_without_equivalents(const gsl_vector * const * observations);
             * @param observations : observations
             * @brief
             * Cette méthode effectue le filtrage avec les paramètres réduits.
             */
            void filter_without_equivalents(const gsl_vector * const * observations);

            /**@fn void tkalman_robust_filter :: smooth_without_equivalents(const gsl_vector * const * observations);
             * @param observations : observations
             * @brief
             * Cette méthode effectue le lissage avec les paramètres réduits.
             */
			void smooth_without_equivalents(const gsl_vector * const * observations);

			/**@fn void tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
			 * @brief
			 * Cette méthode calcule les p_p et p_f équivalents.
			 */
			void compute_equivalents_sqrt_p_p_and_sqrt_p_f();

			/**@fn void tkalman_robust_filter :: compute_equivalents_p_s()
			 * @brief
			 * Cette méthode calcule les p_s équivalents.
			 */
			void compute_equivalents_sqrt_p_s();

			/**@fn void tkalman_robust_filter :: compute_constants()
			 * @brief
			 * Cette méthode calcule les constantes associées aux paramètres du filtre de Kalman couple.
			 *
			 */
			void compute_constants();


			/**@fn void tkalman_robust_filter :: initialize();
			 * @brief
			 Cette méthode met tous les attributs à zéro.
			**/
			void initialize();

			/**@fn void tkalman_robust_filter :: initialize_tmp();
			 * @brief
			 Cette méthode met tous les temporaires à zéros
			**/
			void initialize_tmp();

            /**@fn void tkalman_robust_filter :: initialize_params()
             * @brief
             Cette méthode met les paramètres à zéro.
            **/
            void initialize_params();


			/**@fn int tkalman_robust_filter :: alloc();
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les attributs de l'objet.
			 */
			int alloc();

			/**@fn int tkalman_robust_filter :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

            /**@fn int tkalman_robust_filter :: alloc_params()
             * @brief
             Cette méthode alloue les attributs stockant les paramètres du filtre de Kalman triplet.
            **/
            int alloc_params();


			/**@fn void tkalman_robust_filter :: free();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
			**/
			void free();

			/**@fn void tkalman_robust_filter :: free_tmp();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
			**/
			void free_tmp();

            /**@fn void tkalman_robust_filter :: free_params()
             * @brief
             Cette méthode désalloue les paramètres
            **/
            void free_params();




			/**@fn bool tkalman_robust_filter :: check_tmp() const;
			 * @return
			 * - 0 si les tmps sont bien alloués
			 * - 1 sinon
			 * @brief
			 * Cette fonction controle l'allocation mémoire des tmp.
			 */
			bool check_tmp() const;


		/**@fn bool tkalman_robust_filter :: check_params() const;
			 * @return
			 * - 0 si les tmps sont bien alloués
			 * - 1 sinon
			 * @brief
			 * Cette fonction controle l'allocation mémoire des paramètres.
			 */
			bool check_params() const;


			//Params
			gsl_matrix * sqrt_q_yy;


            //tmp
            gsl_matrix * mat_tt_1;
                gsl_matrix mat_tt_1_yy;
                gsl_matrix mat_tt_1_yx;
                gsl_matrix mat_tt_1_xy;
                gsl_matrix mat_tt_1_xx;
            gsl_matrix * mat_xy_1;
            gsl_matrix * mat_3x2x_1;
                gsl_matrix mat_3x2x_1_view_mat_2xx;
                gsl_matrix mat_3x2x_1_view_00;
                gsl_matrix mat_3x2x_1_view_01;
                gsl_matrix mat_3x2x_1_view_10;
                gsl_matrix mat_3x2x_1_view_11;
                gsl_matrix mat_3x2x_1_view_20;
                gsl_matrix mat_3x2x_1_view_21;

            gsl_permutation * perm_y_1;
            gsl_vector * vect_t_1;
            gsl_vector * vect_2x_1;





    };
#endif // TKALMAN_ROBUST_FILTER_HPP_INCLUDED
