/**@file tkalman_nc_filtering.hpp
 * @author Valérian Némesin
 * 
 */
#ifndef _TKALMAN_NC_FILTERING_HPP_
	#define _TKALMAN_NC_FILTERING_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	using namespace std;
	/**@class tkalman_nc_filtering
	 * @brief
	 * Cette classe permet de réaliser le filtrage du filtre de Kalman triplet sans corrélation
	 * - @fn void tkalman_nc_prediction :: compute_filtering_0
	 * - @fn void tkalman_nc_prediction :: compute_filtering
	 */
	class tkalman_nc_filtering
	{
		public: 
			/**@fn tkalman_nc_filtering :: tkalman_nc_filtering(const gsl_matrix * f_yt,
																 const gsl_matrix * sqrt_q_yy) throw(exception &);
			 * @param[in] f_yt : F^{y,t}
			 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_filtering
			 * @throw 
			 * Exception (std :: bad_alloc si problème de mémoire ou invalid_argument en cas d'arguments invalides)
			 */
			tkalman_nc_filtering(const gsl_matrix * f_yt,
								 const gsl_matrix * sqrt_q_yy) throw(exception &);
			
			/**@fn void tkalman_nc_filtering :: setup(const gsl_matrix * f_yt,
													  const gsl_matrix * sqrt_q_yy) throw(exception &);
			 * @param[in] f_yt : F^{y,t}
			 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @brief
			 * Cette fonction libère les attributs et les réalloue.
			 */
			void setup(const gsl_matrix * f_yt,
					   const gsl_matrix * sqrt_q_yy) throw(exception &);
			
			/**@fn bool tkalman_nc_filtering :: operator !() const
			 * @return 
			 * - 0 si l'objet est normal
			 * - 1 sinon.
			 * @brief
			 * Check de l'objet.
			 **/
			virtual bool operator !() const;
			
			/**@fn void tkalman_nc_filtering :: compute_filtering_0(gsl_vector * t_0_f,
																	gsl_matrix * sqrt_q_0_f,
																	gsl_vector * innovation,
																	gsl_matrix * sqrt_s_0,
																	const gsl_vector * t_0,
																	const gsl_matrix * sqrt_q_0,
																	const gsl_vector * y_0);
			  * @param t_f_0 : espérance de t0 filtré
			  * @param sqrt_q_f_0, : racine de la matrice de covariance de t0 filtré
			  * @param innovation : espérance de l'innovation
			  * @param sqrt_s_0 : racine de la covariance de l'innovation.
			  * @param[in] t_0, : \hat{t}_{0}, espérance de l'état prédit courant
			  * @param[in] sqrt_q_0 :  [Q_{0}]^{\frac{1}{2}}
			  * @param[in] y_0 : y_0, espérance de l'observation courante
			  * @brief
			  * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
			  Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation_0. Puis dans un second temps, nous calculons les  racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_t_0_f.
			 **/
			void compute_filtering_0(gsl_vector * t_0_f,
									 gsl_matrix * sqrt_q_0_f,
									 gsl_vector * innovation,
									 gsl_matrix * sqrt_s_0,
									 const gsl_vector * t_0,
									 const gsl_matrix * sqrt_q_0,
									 const gsl_vector * y_0);
			
			/**@fn void tkalman_nc_filtering :: compute_filtering(gsl_vector * x_f,
																  gsl_matrix * sqrt_p_f,
																  gsl_vector * innovation,
																  gsl_matrix * sqrt_s,
																  const gsl_vector * x_p,
																  const gsl_matrix * sqrt_p_p,
																  const gsl_vector * y,
																  const gsl_vector * _y);
			   * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
			   * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
			   * @param innovation : \tilde{y}_{n}, innovation
			   * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
			   * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit courant
			   * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
			   * @param[in] y : observation cournate
			   * @param[in] _y : observation précédente
			   * @brief
			   * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
			   Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_pf_sqrt_s_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_x_f.
			 **/
			void compute_filtering(gsl_vector * x_f,
								   gsl_matrix * sqrt_p_f,
								   gsl_vector * innovation,
								   gsl_matrix * sqrt_s,
								   const gsl_vector * x_p,
								   const gsl_matrix * sqrt_p_p,
								   const gsl_vector * y,
								   const gsl_vector * _y);
			
			/**@fn tkalman_nc_filtering :: ~tkalman_nc_filtering();
			 * @brief
			 * Destructeur de la classe @class tkalman_nc_filtering
			 */
			~tkalman_nc_filtering();
			
			
			//Accéseurs
			/**@fn inline unsigned int tkalman_nc_filtering :: size_x() const
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int tkalman_nc_filtering :: size_y() const
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int tkalman_nc_filtering :: size_t() const
			 * @return 
			 * Dim. de t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filtering :: f_yt() const
			 * @return 
			 * F^{y,t}
			 **/
			inline const gsl_matrix * f_yt() const
			{
				return _f_yt;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filtering :: f_yx() const
			 * @return 
			 * F^{y,x}
			 **/
			inline const gsl_matrix * const f_yx() const
			{
				return &_f_yx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filtering :: f_yy() const
			 * @return 
			 * F^{y,y}
			 **/
			inline const gsl_matrix * const f_yy() const
			{
				return &_f_yy;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_filtering :: sqrt_q_yy() const
			 * @return 
			 * [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure.
			 **/
			inline const gsl_matrix * sqrt_q_yy() const
			{
				return _sqrt_q_yy;
			}
		
		protected: 
			
			/**@fn void tkalman_nc_filtering :: initialize();
			 * @brief
			 * Cette fonction met tous les attributs de l'objet à 0.
			 */
			void initialize();
			
			/**@fn void tkalman_nc_filtering :: alloc() throw(exception &);
			 * @brief
			 * Cette fonction alloue les différents élements de la classe.
			 * @throw
			 * bad_alloc en cas de problème de mémoire
			 */
			void alloc() throw(exception &);
			
			/**@fn void tkalman_nc_filtering :: free();
			 * @brief
			 * Cette fonction désalloue tous les attributs alloués.
			 **/
			void free();
			
			/**@fn void tkalman_nc_filtering :: create_views();
			 * @brief 
			 * Cette fonction génère les différentes vues sur les matrices.
			 */
			void create_views();
			
		//Tmps
			//Gains 
			gsl_matrix * mat_xy;
			gsl_matrix * mat_ty;
		
			//Matrice de taille (t, t)
			gsl_matrix * mat_tt,
						 mat_tt_yy,
						 mat_tt_yx,
						 mat_tt_xy,
						 mat_tt_xx;
			//Matrice de taille (y + t, y + t)
			gsl_matrix * mat_tpy_tpy,
						 mat_tpy_tpy_view_00,
						 mat_tpy_tpy_view_01,
						 mat_tpy_tpy_view_10,
					     mat_tpy_tpy_view_11;
			//Vecteurs
			gsl_vector * vect_tpy,
					   * vect_t;
					   
			//Perm
			gsl_permutation * perm_y;
			
		//Paramètres
			unsigned int _size_x;
			unsigned int _size_y;	
			unsigned int _size_t;
			const gsl_matrix * _f_yt;
				gsl_matrix _f_yx;
				gsl_matrix _f_yy;
			const gsl_matrix * _sqrt_q_yy;
	};
	
	/**@fn void tkalman_nc_get_innovation(gsl_vector * innovation,
										  const gsl_vector * x_p,
										  const gsl_vector * _y,
										  const gsl_vector * y,
										  const gsl_matrix * f_yx,
										  const gsl_matrix * f_yy)
	 * @param innovation : \tilde{y}_{n}, innovation
	 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit
	 * @param[in] _y : \hat{y}_{n - 1}, espérance de l'observation précédente
	 * @param[in] y : \hat{y}_n, espérance de l'observation courante
	 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
	 * @param[in] f_yy : F^{y,y}, terme de la matrice d'évolution
	 * @brief
	 Cette fonction calcule l'innovation selon la formule :
	 \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
	 */
	void tkalman_nc_get_innovation(gsl_vector * innovation,
								   const gsl_vector * x_p,
								   const gsl_vector * _y,
								   const gsl_vector * y,
								   const gsl_matrix * f_yx,
								   const gsl_matrix * f_yy);

	/**@fn void tkalman_get_innovation_0(gsl_vector * innovation,
										 const gsl_vector * t_0,
										 const gsl_vector * y_0,
										 const gsl_matrix * f_yt)
	 * @param innovation : \tilde{y}_{0}, innovation
	 * @param[in] t_p : \hat{t}_{0}
	 * @param[in] y_0 : y_0, espérance de l'observation courante
	 * @param[in] f_yt : F^{y,t}, terme de la matrice d'évolution
	 * @brief
	 Cette fonction calcule l'innovation selon la formule :
	 \tilde{y}_{0} = y_{0} - F^{y,t} \hat{t}_{0}
	 */
	void tkalman_nc_get_innovation_0(gsl_vector * innovation,
								  const gsl_vector * t_0,
								  const gsl_vector * y_0,
								  const gsl_matrix * f_yt);
	//Espérance de l'état filtré
	/**@fn void tkalman_nc_get_x_f(gsl_vector * x_f,
										 const gsl_vector * x_p,
										 const gsl_vector * innovation,
										 const gsl_matrix * gain)
	 * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
	 * @param[in] x_p : \hat{x}_{n|n - 1}, espérance de l'état prédit
	 * @param[in] innovation : \tilde{y}_{n}, innovation
	 * @param[in] gain : K_{n|n}, gain de filtrage
	 * @brief
	 Cette fonction calcule l'espérance de l'état filtré selon la formule :
	 \hat{x}_{n|n} = \hat{x|n-1} + K_{n|n} \; \tilde{y}_n
	 */
	void tkalman_nc_get_x_f(gsl_vector * x_f,
							const gsl_vector * x_p,
							const gsl_vector * innovation,
							const gsl_matrix * gain);
	/**@fn void tkalman_nc_get_t_0_f(gsl_vector * t_0_f,
									 const gsl_vector * t_0,
									 const gsl_vector * innovation_0,
									 const gsl_matrix * gain)
	 * @param t_0_f: \hat{t}_{0|0}
	 * @param[in] t_0 : \hat{t}_0,
	 * @param[in] innovation_0 : \tilde{y}_{0}, innovation
	 * @param[in] gain : K^t_{0|0}, gain de filtrage
	 * @brief
	 Cette fonction calcule l'espérance de l'état filtré selon la formule :
	 \hat{t}_{0|0} = \hat{t}_0 + K^t_{0|0} \; \tilde{y}_{0}
	 */
	void tkalman_nc_get_t_0_f(gsl_vector * t_0_f,
							  const gsl_vector * t_0,
							  const gsl_vector * innovation_0,
							  const gsl_matrix * gain);
	//Racine de la matrice de cov. de l'état filtré
	/**@fn void tkalman_nc_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
													   gsl_matrix * sqrt_s,
													   gsl_matrix * gain,
													   const gsl_matrix * sqrt_p_p,
													   const gsl_matrix * f_yx,
													   const gsl_matrix * sqrt_q_yy,
													   gsl_matrix * mat_tt,
													   gsl_matrix * mat_tt_yy,
													   gsl_matrix * mat_tt_yx,
													   gsl_matrix * mat_tt_xy,
													   gsl_matrix * mat_tt_xx,
													   gsl_permutation * perm_y,
													   gsl_vector * vect_t);
	 * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
	 * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
	 * @param gain : K_{n|n} = P_{n|n - 1} \; [F^{y,x}]^T \; S_{n}^{-1},  gain de filtrage
	 * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
	 * @param[in] f_yx : F^{y,x}, terme de la matrice de transition
	 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
	 * @param mat_tt : matrice de taille  (n_t, n_t) préallouée.
	 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y - 1, n_y - 1))
	 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y - 1, n_t - 1))
	 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t - 1, n_y - 1))
	 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t - 1, n_t - 1))
	 * @param perm_y : permutation de taille y préallouée
	 * @param vect_t : vecteur de taille n_t préalloué pour le calcul
	 * @brief
	 Cette fonction calcule les racines de la covariance de l'état filtré courant et de la matrice de covariance de l'innovation. \n
	ECe calcul s'effectue en plusieurs étapes : nous construisons la matrice M : \n
	+---------------------------+
	|sqrt_q_yy      0           |\n
	|                           |\n
	|sqrt_pp.Fyx^T  sqrt_p_p    |\n
	+---------------------------+\n
	puis nous effectuons sa décomposition QR. A partir de la matrice R de la décomposition, nous obtenons :
	+---------------------------------------+
	|sqrt_s      sqrt_s^-1 K^T             |\n
	|                           			|\n
	|0           sqrt_p_f       		|\n
	+---------------------------------------+\n
	 */
	void tkalman_nc_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
												gsl_matrix * sqrt_s,
												gsl_matrix * gain,
												const gsl_matrix * sqrt_p_p,
												const gsl_matrix * f_yx,
												const gsl_matrix * sqrt_q_yy,
												gsl_matrix * mat_tt,
												gsl_matrix * mat_tt_yy,
												gsl_matrix * mat_tt_yx,
												gsl_matrix * mat_tt_xy,
												gsl_matrix * mat_tt_xx,
												gsl_permutation * perm_y,
												gsl_vector * vect_t);
	//Racine de la matrice de cov. de l'état filtré
	/**@fn void tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain(gsl_matrix * sqrt_q_0_f,
															gsl_matrix * sqrt_s_0,
															gsl_matrix * gain,
															const gsl_matrix * sqrt_q_0,
															const gsl_matrix * f_yt,
															const gsl_matrix * sqrt_q_yy,
															gsl_matrix * mat_tpy_tpy,
															gsl_matrix * mat_tpy_tpy_view_00,
															gsl_matrix * mat_tpy_tpy_view_01,
															gsl_matrix * mat_tpy_tpy_view_10,
															gsl_matrix * mat_tpy_tpy_view_11,
															gsl_permutation * perm_y,
															gsl_vector * vect_tpy)
	 * @param sqrt_q_0_f :  [Q_{0|0}]^{\frac{1}{2}}
	 * @param sqrt_s_0 : [S_{0}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
	 * @param gain : K^t_{0|0} = Q_{0} \; [F^{y,t}]^T \; S_{0}^{-1},  gain de filtrage
	 * @param[in] sqrt_q_0 :  [Q_{0}]^{\frac{1}{2}}
	 * @param[in] f_yt : F^{y,t}, terme de la matrice de transition
	 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
	 * @param mat_tpy_tpy : matrice de taille (n_t + n_y, n_t + n_y)  préallouée.
	 * @param mat_tpy_tpy_view_00 : vue sur la matrice mat_tpy_tpy (de (0,0) à (n_y - 1, n_y - 1))
	 * @param mat_tpy_tpy_view_01 : vue sur la matrice mat_tpy_tpy (de (0, n_y) à (n_y - 1, n_t + n_y - 1))
	 * @param mat_tpy_tpy_view_10 : vue sur la matrice mat_tpy_tpy (n_y, 0) à (n_t + n_y - 1, n_y - 1))
	 * @param mat_tpy_tpy_view_11 : vue sur la matrice mat_tpy_tpy (n_y, n_y) à (n_t + n_y - 1,  n_t + n_y - 1))
	 * @param perm_y : permutation de taille y préallouée
	 * @param vect_tpy : vecteur de taille (n_t + n_y) préalloué pour le calcul
	 * @brief
	 Cette fonction calcule les racines de la covariance de l'état filtré courant et de la matrice de covariance de l'innovation. \n
	ECe calcul s'effectue en plusieurs étapes : nous construisons la matrice M : \n
	+---------------------------+
	|sqrt_q_yy      0           |\n
	|                           |\n
	|sqrt_q0.Fyt^T  sqrt_q0    |\n
	+---------------------------+\n
	puis nous effectuons sa décomposition QR. A partir de la matrice R de la décomposition, nous obtenons :
	+---------------------------------------+
	|sqrt_s0      sqrt_s0^-1 K^T             |\n
	|                           			|\n
	|0           sqrt_q0_f       		|\n
	+---------------------------------------+\n
	 */
	void tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain(gsl_matrix * sqrt_q_0_f,
													 gsl_matrix * sqrt_s_0,
													 gsl_matrix * gain,
													 const gsl_matrix * sqrt_q_0,
													 const gsl_matrix * f_yt,
													 const gsl_matrix * sqrt_q_yy,
													 gsl_matrix * mat_tpy_tpy,
													 gsl_matrix * mat_tpy_tpy_view_00,
													 gsl_matrix * mat_tpy_tpy_view_01,
													 gsl_matrix * mat_tpy_tpy_view_10,
													 gsl_matrix * mat_tpy_tpy_view_11,
													 gsl_permutation * perm_y,
													 gsl_vector * vect_tpy);
	//Filtrage complet
	/**@fn void tkalman_nc_do_filtering(gsl_vector * x_f,
										gsl_matrix * sqrt_p_f,
										gsl_vector * innovation,
										gsl_matrix * sqrt_s,
										const gsl_vector * x_p,
										const gsl_matrix * sqrt_p_p,
										const gsl_vector * y,
										const gsl_vector * _y,
										const gsl_matrix * f_yx,
										const gsl_matrix * f_yy,
										const gsl_matrix * sqrt_q_yy,
										gsl_matrix * mat_tt,
										gsl_matrix * mat_tt_yy,
										gsl_matrix * mat_tt_yx,
										gsl_matrix * mat_tt_xy,
										gsl_matrix * mat_tt_xx,
										gsl_matrix * mat_xy,
										gsl_permutation * perm_y,
										gsl_vector * vect_t)
	 * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
	 * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
	 * @param innovation : \tilde{y}_{n}, innovation
	 * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
	 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit courant
	 * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
	 * @param[in] y : observation cournate
	 * @param[in] _y : observation précédente
	 * @param[in] f_yx : terme de la matrice de transition (Fyx)
	 * @param[in] f_yy : terme de la matrice de transition (Fyy)
	 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
	 * @param mat_tt : matrice de taille  (n_t, n_t) préallouée.
	 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y - 1, n_y - 1))
	 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y - 1, n_t - 1))
	 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t - 1, n_y - 1))
	 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t - 1, n_t - 1))
	 * @param mat_xy : matrice de taille (n_x.n_y) préallouée
	 * @param perm_y : permutation de taille n_y préallouée
	 * @param vect_t : vecteur de taille n_t préalloué pour le calcul
	 * @brief
	 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
	 Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_pf_sqrt_s_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_x_f.
	 */
	void tkalman_nc_do_filtering(gsl_vector * x_f,
								 gsl_matrix * sqrt_p_f,
								 gsl_vector * innovation,
								 gsl_matrix * sqrt_s,
								 const gsl_vector * x_p,
								 const gsl_matrix * sqrt_p_p,
								 const gsl_vector * y,
								 const gsl_vector * _y,
								 const gsl_matrix * f_yx,
								 const gsl_matrix * f_yy,
								 const gsl_matrix * sqrt_q_yy,
								 gsl_matrix * mat_tt,
								 gsl_matrix * mat_tt_yy,
								 gsl_matrix * mat_tt_yx,
								 gsl_matrix * mat_tt_xy,
								 gsl_matrix * mat_tt_xx,
								 gsl_matrix * mat_xy,
								 gsl_permutation * perm_y,
								 gsl_vector * vect_t);
								 
	/**@fn void tkalman_nc_do_filtering_0(gsl_vector * t_0_f,
										  gsl_matrix * sqrt_q_0_f,
										  gsl_vector * innovation,
										  gsl_matrix * sqrt_s_0,
										  const gsl_vector * t_0,
										  const gsl_matrix * sqrt_q_0,
										  const gsl_vector * y_0,
										  const gsl_matrix * f_yt,
										  const gsl_matrix * sqrt_q_yy,
										  gsl_matrix * mat_tpy_tpy,
										  gsl_matrix * mat_tpy_tpy_view_00,
										  gsl_matrix * mat_tpy_tpy_view_01,
										  gsl_matrix * mat_tpy_tpy_view_10,
										  gsl_matrix * mat_tpy_tpy_view_11,
										  gsl_matrix * mat_ty,
										  gsl_permutation * perm_y,
										  gsl_vector * vect_tpy);
	 * @param t_f_0 : espérance de t0 filtré
	 * @param sqrt_q_f_0, : racine de la matrice de covariance de t0 filtré
	 * @param innovation : espérance de l'innovation
	 * @param sqrt_s_0 : racine de la covariance de l'innovation.
	 * @param[in] t_0, : \hat{t}_{0}, espérance de l'état prédit courant
	 * @param[in] sqrt_q_0 :  [Q_{0}]^{\frac{1}{2}}
	 * @param[in] y_0 : y_0, espérance de l'observation courante
	 * @param[in] f_yt : F^{y,t}, terme de la matrice de transition
	 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
	 * @param mat_tpy_tpy : matrice de taille (n_t + n_y.n_t + n_y) préallouée.
	 * @param mat_tpy_tpy_view_00 : vue sur la matrice mat_tpy_tpy (de (0,0) à (n_y - 1, n_y - 1))
	 * @param mat_tpy_tpy_view_01 : vue sur la matrice mat_tpy_tpy (de (0,n_y) à (n_y - 1, n_y + n_t - 1))
	 * @param mat_tpy_tpy_view_10 : vue sur la matrice mat_tpy_tpy (de (n_y, 0) à (n_y + n_t - 1, n_y - 1))
	 * @param mat_tpy_tpy_view_11 : vue sur la matrice mat_tpy_tpy (de (n_y, n_y) à (n_y + n_t - 1, n_y + n_t - 1))
	 * @param mat_ty : matrice de taille (n_t.n_y) préallouée
	 * @param perm_y : permutation de taille y préallouée
	 * @param vect_tpy : vecteur de taille (n_t + n_y) préalloué pour le calcul
	 * @brief
	 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
	 Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation_0. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_t_0_f.
	 */
	void tkalman_nc_do_filtering_0(gsl_vector * t_0_f,
								   gsl_matrix * sqrt_q_0_f,
								   gsl_vector * innovation,
								   gsl_matrix * sqrt_s_0,
								   const gsl_vector * t_0,
								   const gsl_matrix * sqrt_q_0,
								   const gsl_vector * y_0,
								   const gsl_matrix * f_yt,
								   const gsl_matrix * sqrt_q_yy,
								   gsl_matrix * mat_tpy_tpy,
								   gsl_matrix * mat_tpy_tpy_view_00,
								   gsl_matrix * mat_tpy_tpy_view_01,
								   gsl_matrix * mat_tpy_tpy_view_10,
								   gsl_matrix * mat_tpy_tpy_view_11,
								   gsl_matrix * mat_ty,
								   gsl_permutation * perm_y,
								   gsl_vector * vect_tpy);
	
	
	
#endif
