#ifndef TKALMAN_ROBUST_EM_HPP_INCLUDED
#define TKALMAN_ROBUST_EM_HPP_INCLUDED
    #include "tkalman_robust_filter.hpp"
	#include "tkalman_em_algorithm.hpp"
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	class tkalman_robust_em : public tkalman_robust_filter
	{
		public:
			/**@fn tkalman_robust_em :: tkalman_robust_em (const gsl_vector * x0,
															  const gsl_matrix * sqrt_p0,
															  const gsl_matrix * f,
															  const gsl_matrix * sqrt_q,
															  unsigned int n = 0,
															  unsigned int p = 0,
															  bool data = false);
			 * @param[in] x0 : Esp�rance de l'�tat initial
			 * @param[in] sqrt_p0 : Matrice de covariance de l'�tat initial (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'�volution
			 * @param[in] sqrt_q : Matrice de covariance (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par d�faut)
			 * @param[in] p : Nombre d'it�rations de l'EM (0 par d�faut et dans ce cas, cela �quivaut � un filtrage simple)
			 * @param[in] data : Bool�en (True = suivi de l'EM et stockage des param�tres et de la vraisemblance � chaque it�ration, False = Pas de suivi de l'EM)
			 * @brief
			 * Constructeur
			*/
			tkalman_robust_em(const gsl_vector * x0,
                              const gsl_matrix * sqrt_p0,
                              const gsl_matrix * f,
                              const gsl_matrix * sqrt_q,
                              unsigned int n = 0,
                              unsigned int p = 0,
                              bool data = false);

			/**@fn virtual int tkalman_robust_em  :: setup(const gsl_vector * x0,
															const gsl_matrix * p0,
															const gsl_matrix * f,
															const gsl_matrix * q,
															unsigned int n = 0,
															unsigned int p = 0,
															bool data = false);
			 * @param[in] x0 : Esp�rance de l'�tat initial
			 * @param[in] p0 : Matrice de covariance de l'�tat initial (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'�volution
			 * @param[in] q : Matrice de covariance (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par d�faut)
			 * @param[in] p : Nombre d'it�rations de l'EM (0 par d�faut et dans ce cas, cela �quivaut � un filtrage simple)
			 * @param[in] data : Bool�en (True = suivi de l'EM et stockage des param�tres et de la vraisemblance � chaque it�ration, False = Pas de suivi de l'EM)
			 * @brief
			 * Cette m�thode r�initialise l'objet.
			 */
			virtual int setup(const gsl_vector * x0,
							  const gsl_matrix * p0,
							  const gsl_matrix * f,
							  const gsl_matrix * q,
							  unsigned int n = 0,
							  unsigned int p = 0,
							  bool data = false);


			/**@fn virtual void tkalman_robust_em  :: filter(const gsl_vector * const * observations,
															  unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette m�thode effectue le filtrage des donn�es par le filtre de Kalman Triple non supervis�.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual void ttkalman_robust_em :: smooth(const gsl_vector * const * observations,
															  unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette m�thode effectue le lissage des donn�es par le filtre de Kalman Triple non supervis�.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn tkalman_robust_em  :: ~tkalman_robust_em()
			 * @brief
			 * Destructeur
			 *
			 */
			~tkalman_robust_em(void);

			/**@fn bool tkalman_robust_em :: operator!() const
			 * @return
			 - 0 si l'objet est valide
			 - 1 sinon
			 * @brief
			 Cette m�thode teste la validit� de chaque attribut.
			**/
			virtual bool operator!() const;

			//Accesseurs
			/**@fn const gsl_vector * const * tkalman_robust_em  :: get_x0_est()
			 * @return
			 Liste des x0 estim�s
			 * @brief
			 Cette fonction renvoie les diff�rents x0 estim�s par l'algorithme EM.
			**/
			inline const gsl_vector * const * get_x0_est() const
			{
                return _x0_est;
			}

			/**@fn const gsl_matrix * const * tkalman_robust_em :: get_sqrt_p0_est()
			 * @return
			 Liste des p0 estim�s
			 * @brief
			 Cette fonction renvoie les diff�rents p0 estim�s par l'algorithme EM.
			**/
			inline  const gsl_matrix * const * get_sqrt_p0_est() const
			{
                return _p0_est;
			}

			/**@fn const gsl_matrix * const * tkalman_robust_em :: get_f_est()
			 * @return
			 Liste des F estim�s
			 * @brief
			 Cette fonction renvoie les diff�rents F estim�s par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_f_est() const
			{
                return _f_est;
			}

			/**@fn const gsl_matrix * const * tkalman_robust_em :: get_sqrt_q_est()
			 * @return
			 Liste des Q estim�s
			 * @brief
			 Cette fonction renvoie les diff�rents Q estim�s par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_sqrt_q_est() const
			{
                return _q_est;
			}
			/**@fn inline const double * tkalman_robust_em :: get_likelihood()
			 * @return
			 Valeurs de la log-vraisemblance au cours des it�rations de l'EM
			 * @brief
			 Cette fonction renvoie les valeurs de la log vraisemblance au cours des it�rations de l'EM.
			**/
			inline const double * get_likelihood() const
			{
				return _log_likelihood;
			}

			/**@fn inline unsigned int tkalman_robust_em :: get_nb_iter()
			 * @return
			 Nombre d'it�rations de l'algorithme EM
			 @brief
			 Cette fonction renvoie le nombre d'it�rations de l'algorithme EM.
			**/
			inline unsigned int get_nb_iter() const
			{
				return nb_iter;
			}
		protected:
			/**@fn void tkalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations);
			 * @param[in] observations : observations
			 * @brief
			 * Cette m�thode estime les param�tres du filtre de Kalman couple � partir d'un jeu d'observations.
			 */
			void do_em_algorithm(const gsl_vector * const * observations);

			/**@fn void tkalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
																   unsigned int i);
			 * @param[in] observations : observations
			 * @param[in] i : it�ration de l'EM
			 * @brief
			 * Cette m�thode r�-estime les param�tres du filtre de Kalman couple � partir d'un jeu d'observations.
			 */
			void estimate_parameters(const gsl_vector * const * observations,
									 unsigned int i);



			/**@fn void tkalman_robust_em :: follow(unsigned int i);
			 * @param i : it�ration de l'EM
			 * @brief
			 * Cette m�thode enregistre les param�tres de l'it�ration i de l'algorithme EM.
			 **/
			void follow(unsigned int i);


			/**@fn void tkalman_robust_em :: initialize();
			 * @brief
			 Cette m�thode met tous les attributs � z�ro.
			**/
			void initialize();

			/**@fn void tkalman_robust_em :: initialize_tmp();
			 * @brief
			 Cette m�thode met tous les temporaires � z�ros
			**/
			void initialize_tmp();

			/**@fn int tkalman_robust_em :: alloc(unsigned int p, bool data);
			 * @param p : nombre d'it�rations de l'EM
			 * @param data : bool�en
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien d�roul�e
			  - 1 sinon
			  *@brief
			  Cette m�thode alloue les attributs de l'objet.
			 */
			int alloc(unsigned int p, bool data);

			/**@fn int tkalman_robust_em :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien d�roul�e
			  - 1 sinon
			  *@brief
			  Cette m�thode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

			/**@fn int tkalman_robust_em :: alloc_data()
			 * @brief
			 * Cette m�thode alloue les donn�es
			 *
			 */
			int alloc_data();

			/**@fn void tkalman_robust_em :: free();
			 * @brief
			 Cette m�thode lib�re la m�moire utilis�e par les attributs de l'objet.
			**/
			void free();

			/**@fn void tkalman_robust_em :: free_tmp();
			 * @brief
			 Cette m�thode lib�re la m�moire utilis�e par les temporaires de l'objet.
			**/
			void free_tmp();


			/**@fn void tkalman_robust_em :: free_data();
			 *  @brief
			 Cette m�thode lib�re la m�moire utilis�e par les donn�es de suivi de l'EM.
			**/
			void free_data();

			/**@fn bool tkalman_robust_em :: check_tmp()
			 * @return
			 - 0 si les tmp ont �t� bien allou�s
			 - 1 sinon
			 * @brief
			 Cette m�thode teste la bonne allocation des variables temporaires utilis�es par l'objet.
			**/
			bool check_tmp() const;

			/**@fn bool tkalman_robust_em :: check_data()
			 * @return
			 - 0 si les donn�es de suivi de l'EM ont �t� bien allou�s
			 - 1 sinon
			 * @brief
			 Cette m�thode teste la bonne allocation des donn�es de suivi de l'EM.
			**/
			bool check_data() const;





			//Nombre d'it�ration de l'EM
				unsigned int nb_iter;


			//Suivi de l'EM
				double * _log_likelihood;
				gsl_vector ** _x0_est;
				gsl_matrix ** _p0_est;
				gsl_matrix ** _f_est;
				gsl_matrix ** _q_est;


			//Tmp de l'EM
                gsl_matrix * mat_4x4x_1;
                      gsl_matrix mat_4x4x_1_view_00;
                      gsl_matrix mat_4x4x_1_view_11;
                      gsl_matrix mat_4x4x_1_view_02;
                      gsl_matrix mat_4x4x_1_view_20;
                      gsl_matrix mat_4x4x_1_view_22;
                      gsl_matrix mat_4x4x_1_view_23;
                      gsl_matrix mat_4x4x_1_view_32;
                      gsl_matrix mat_4x4x_1_view_33;
                gsl_matrix * mat_4tp12t_1;
                    gsl_matrix mat_4tp12t_1_view_mat_4t2t;
                    gsl_matrix mat_4tp12t_1_view_10b;
                    gsl_matrix mat_4tp12t_1_view_10b_view_00;
                    gsl_matrix mat_4tp12t_1_view_10b_view_02;
                    gsl_matrix mat_4tp12t_1_view_10b_view_22;
                    gsl_vector mat_4tp12t_1_view_10b_view_30;
                    gsl_vector mat_4tp12t_1_view_10b_view_31;
                    gsl_vector mat_4tp12t_1_view_10b_view_32;
                    gsl_vector mat_4tp12t_1_view_10b_view_33;
                        gsl_matrix mat_4tp12t_1_view_sqrt_sum_view_00;
                        gsl_matrix mat_4tp12t_1_view_sqrt_sum_view_01;
                        gsl_matrix mat_4tp12t_1_view_sqrt_sum_view_11;

                  gsl_vector * vect_2t_1;
                  gsl_vector * vect_4x_1;
                  gsl_permutation * perm_t_1;

            //tmp du suivi
                    gsl_permutation * perm_y_1;













	};


#endif // TKALMAN_ROBUST_EM_HPP_INCLUDED
