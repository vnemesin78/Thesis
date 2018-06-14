/**@file tkalman_nc_smoothing.hpp
 * @author Valérian Némesin
 * 
 */
#ifndef _TKALMAN_NC_SMOOTHING_HPP_
	#define _TKALMAN_NC_SMOOTHING_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	using namespace std;
	/**@class tkalman_nc_smoothing
	 * @brief
	 * Cette classe permet de réaliser le lissage dans le filtre de Kalman triplet non corrélé. \newline
	 * Les méthodes à utiliser sont : 
	 * @fn compute_smoothing_0
	 * @fn compute_smoothing
	 */
	 class tkalman_nc_smoothing
	 {
		public:
			/**@fn tkalman_nc_smoothing :: tkalman_nc_smoothing(const gsl_matrix * f_xt,
															 const gsl_matrix * sqrt_q_xx) throw(exception &);
			* @param[in] f_xt : F^{x,t}
			* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			* @brief
			* Constructeur de la classe @class tkalman_nc_smoothing
			* @throw
			* Exception en cas d'erreur.
			*/
			tkalman_nc_smoothing(const gsl_matrix * f_xt,
								 const gsl_matrix * sqrt_q_xx) throw(exception &);
							  
			/**@fn void tkalman_nc_smoothing :: setup(const gsl_matrix * f_xt,
												   const gsl_matrix * sqrt_q_xx) throw(exception &);
			* @param[in] f_xt : F^{x,t}
			* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			* @brief
			* Réinitialisation de la classe @class tkalman_nc_smoothing
			* @throw
			* Exception en cas d'erreur.
			*/
			void setup(const gsl_matrix * f_xt,
					   const gsl_matrix * sqrt_q_xx) throw(exception &);
			
			/**@fn tkalman_nc_smoothing :: ~ tkalman_nc_smoothing();
			 * @brief
			 * Destructeur de la classe @class tkalman_nc_smoothing.
			 */
			~tkalman_nc_smoothing();
			 
			/**@fn inline unsigned int tkalman_nc_smoothing :: size_x() const
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int tkalman_nc_smoothing :: size_y() const
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int tkalman_nc_smoothing :: size_t() const
			 * @return 
			 * Dim. de t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
		 
			/**@fn inline const gsl_matrix * tkalman_nc_smoothing  :: f_xt() const
			 * @return 
			 * F^{x,t}
			 **/
			inline const gsl_matrix * f_xt() const
			{
				return _f_xt;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_smoothing  :: f_xx() const
			 * @return 
			 * F^{x,x}
			 **/
			inline const gsl_matrix * f_xx() const
			{
				return &_f_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_smoothing  :: sqrt_q_xx() const
			 * @return 
			 * [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			 **/
			inline const gsl_matrix * sqrt_q_xx() const
			{
				return _sqrt_q_xx;
			}
			
			/**@fn void tkalman_nc_smoothing :: compute_smoothing_0( gsl_vector * t_0_s,
																	 gsl_matrix * sqrt_q_0_s,
																	 gsl_matrix * c_s,
																	 const gsl_vector * t_0_f,
																	 const gsl_matrix * sqrt_q_0_f,
																	 const gsl_vector * x_1_p,
																	 const gsl_matrix * sqrt_p_1_p,
																	 const gsl_vector * x_1_s,
																	 const gsl_matrix * sqrt_p_1_s)
			 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
			 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
			 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
			 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
			 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
			 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
			 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
			 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
			 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
			 * @brief
			 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain_0. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_q_0_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_t_0_s .
			 */
			void compute_smoothing_0( gsl_vector * t_0_s,
									  gsl_matrix * sqrt_q_0_s,
									  gsl_matrix * c_s,
									  const gsl_vector * t_0_f,
									  const gsl_matrix * sqrt_q_0_f,
									  const gsl_vector * x_1_p,
									  const gsl_matrix * sqrt_p_1_p,
									  const gsl_vector * x_1_s,
									  const gsl_matrix * sqrt_p_1_s);
			
			/**@fn void tkalman_nc_smoothing  :: compute_smoothing( gsl_vector * x_s,
																	gsl_matrix * sqrt_p_s,
																	gsl_matrix * c_s,
																	const gsl_vector * x_f,
																	const gsl_matrix * sqrt_p_f,
																	const gsl_vector * x_p_,
																	const gsl_matrix * sqrt_p_p_,
																	const gsl_vector * x_s_,
																	const gsl_matrix * sqrt_p_s_)	
			 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
			 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
			 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
			 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
			 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
			 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
			 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
			 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
			 * @param[in] sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
			 * @brief
			 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_p_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_x_s .
			 */
			void compute_smoothing( gsl_vector * x_s,
								    gsl_matrix * sqrt_p_s,
								    gsl_matrix * c_s,
								    const gsl_vector * x_f,
								    const gsl_matrix * sqrt_p_f,
								    const gsl_vector * x_p_,
								    const gsl_matrix * sqrt_p_p_,
								    const gsl_vector * x_s_,
								    const gsl_matrix * sqrt_p_s_);
			/**@fn bool tkalman_nc_smoothing :: operator!() const;
			 * @return
			 * - 0 si l'objet est valide
			 * - 1 sinon
			 * 
			 */			    
			virtual bool operator !() const;	
			
		 protected:
			/**@fn void tkalman_nc_smoothing :: free();
			 * @brief
			 * Cette méthode libère la mémoire occupée par les attributs de l'objet.
			 */
			void free();
			
			/**@fn void tkalman_nc_smoothing :: alloc() throw(exception &);
			 * @brief
			 * Cette méthode alloue la mémoire nécessaire pour les attributs de l'objet
			 * @throw
			 * Exception en cas de prob.
			 */
			void alloc() throw(exception &);
			
			/**@fn void tkalman_nc_smoothing :: initialize();
			 * @brief
			 * Cette méthode met tous les attributs de l'objet à 0.
			 * 
			 */
			void initialize();
			
			/**@fn void tkalman_nc_smoothing :: create_views();
			 * @brief
			 * Cette méthode crée les vues...
			 */
			void create_views();
			
			//Params
			unsigned int _size_x,
						 _size_y,
						 _size_t;
			const gsl_matrix * _f_xt;
				gsl_matrix _f_xx;
			const gsl_matrix * _sqrt_q_xx;
			
			//Tmp
			gsl_matrix * mat_xx;
			gsl_matrix * mat_tx;
			gsl_matrix * mat_2xpt_xpt;
				gsl_matrix mat_2xpt_xpt_view_00;
				gsl_matrix mat_2xpt_xpt_view_01;
				gsl_matrix mat_2xpt_xpt_view_10;
				gsl_matrix mat_2xpt_xpt_view_11;
				gsl_matrix mat_2xpt_xpt_view_20;
				gsl_matrix mat_2xpt_xpt_view_21;

			gsl_matrix * mat_3x2x;
				gsl_matrix mat_3x2x_view_00;
				gsl_matrix mat_3x2x_view_01;
				gsl_matrix mat_3x2x_view_10;
				gsl_matrix mat_3x2x_view_11;
				gsl_matrix mat_3x2x_view_20;
				gsl_matrix mat_3x2x_view_21;
				 
			 gsl_permutation * perm_x;
			 gsl_vector * vect_2x;	
			 gsl_vector * vect_xpt;
		 
	 };
	
	/**@fn void tkalman_nc_get_x_s(gsl_vector * x_s,
								   const gsl_vector * x_f,
								   const gsl_vector * x_p_,
								   const gsl_vector * x_s_,
								   const gsl_matrix * gain)
	 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
	 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
	 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
	 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
	 * @param[in] gain : K_{n|N}, gain de lissage
	 * @brief
	 * Cette fonction calcule l'espérance de l'état lissé selon la formule \n
	 * \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n})
	 * .
	 */
	void tkalman_nc_get_x_s(gsl_vector * x_s,
							const gsl_vector * x_f,
							const gsl_vector * x_p_,
							const gsl_vector * x_s_,
							const gsl_matrix * gain);
	/**@fn void tkalman_nc_get_t_0_s(gsl_vector * t_0_s,
									 const gsl_vector * t_0_f,
									 const gsl_vector * x_1_p,
									 const gsl_vector * x_1_s,
									 const gsl_matrix * gain)
	 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
	 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
	 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
	 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
	 * @param[in] gain : K^t_{0|N}, gain de lissage
	 * @brief
	 * Cette fonction calcule l'espérance de l'état lissé selon la formule \n
	 * \hat{t}_{0|N} = \hat{t}_{0|n} + K^t_{0|N} (\hat{x}_{1|N} - \hat{x}_{1|0})
	 */
	void tkalman_nc_get_t_0_s(gsl_vector * x_s,
							  const gsl_vector * x_f,
							  const gsl_vector * x_p_,
							  const gsl_vector * x_s_,
							  const gsl_matrix * gain);
							  
	/**@fn void tkalman_nc_get_smoothing_gain_0(gsl_matrix * s_gain,
											 const gsl_matrix * sqrt_q_0_f,
											 const gsl_matrix * sqrt_p_1_p,
											 const gsl_matrix * f_xt,
											 gsl_matrix * mat_xx,
											 gsl_matrix * mat_tx,
											 gsl_permutation * perm_x)
	 * @param s_gain : K^t_{0|N}, gain de Kalman triplet pour le lissage
	 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
	 * @param[in] f_xt : F^{x,t}
	 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
	 * @param mat_tx : matrice de taille (n_t.n_x) préallouée
	 * @param perm_x : permutation de taille (n_x) préallouée
	 * @brief
	 * Cette fonction calcule le gain de lissage du filtre de Kalman
	 * triplet : \newline
	 * K^t_{0|N} = Q_{0|0} [F^{x,t}]^T P_{1|0}^{-1}  
	 */
	void tkalman_nc_get_smoothing_gain_0(gsl_matrix * s_gain,
										 const gsl_matrix * sqrt_q_0_f,
										 const gsl_matrix * sqrt_p_1_p,
										 const gsl_matrix * f_xt,
										 gsl_matrix * mat_xx,
										 gsl_matrix * mat_tx,
										 gsl_permutation * perm_x);
										 
	/**@fn void tkalman_nc_get_smoothing_gain(gsl_matrix * s_gain,
										   const gsl_matrix * sqrt_p_f,
										   const gsl_matrix * sqrt_p_p_,
										   const gsl_matrix * f_xx,
										   gsl_matrix * mat_xx,
										   gsl_permutation * perm_x)
	 * @param s_gain : K_{n|N}, gain de Kalman triplet pour le lissage
	 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
	 * @param[in] f_xx : F^{x,x}
	 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
	 * @param perm_x : permutation de taille (n_x) préallouée
	 * @brief
	 * Cette fonction calcule le gain de lissage du filtre de Kalman
	 * triplet : \newline
	 * K_{n|N} = P_{n|n} [F^{x,x}]^T P_{n + 1|n} 
	 */
	void tkalman_nc_get_smoothing_gain(gsl_matrix * s_gain,
									   const gsl_matrix * sqrt_p_f,
									   const gsl_matrix * sqrt_p_p_,
									   const gsl_matrix * f_xx,
									   gsl_matrix * mat_xx,
									   gsl_permutation * perm_x);
									   
	/**@fn void tkalman_nc_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
												gsl_matrix * c_s,
												const gsl_matrix * sqrt_p_f,
												const gsl_matrix * sqrt_p_s_,
												const gsl_matrix * f_xx,
												const gsl_matrix * sqrt_q_xx,
												const gsl_matrix * s_gain,
												gsl_matrix * mat_3x2x,
												gsl_matrix * mat_3x2x_view_00,
												gsl_matrix * mat_3x2x_view_01,
												gsl_matrix * mat_3x2x_view_10,
												gsl_matrix * mat_3x2x_view_11,
												gsl_matrix * mat_3x2x_view_20,
												gsl_matrix * mat_3x2x_view_21,
												gsl_vector * vect_2x)
	 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
	 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
	 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in]  sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
	 * @param[in] f_xx : F^{x,x}
	 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	 * @param[in] s_gain :  K_{n|N}, gain de Kalman triplet pour le lissage
	 * @param mat_3x2x : matrice de taille (3n_x.2n_x) préallouée
	 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x - 1 , n_x - 1)
	 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x - 1, 2 n_x - 1)
	 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x - 1 , n_x - 1)
	 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x - 1, 2n_x - 1)
	 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x - 1, n_x - 1)
	 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x - 1, 2n_x - 1)
	 * @param vect_2x : vecteur de taille (2n_x) préalloué
	 * @brief
	 Cette fonction calcule la racine de la matrice de covariance de l'état lissé courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
	 * On calcule la matrice
	 * \begin{pmatrix}
	 [Q^{x,x}]^{\frac{1}{2}}	&	0 \newline
	 [P_{n|n}]^{\frac{1}{2}}	&	[ F^{x,x}]^T \newline
	 0							&	[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
	 * \end{pmatrix}
	 * On effectue la décomposition QR :
	 * \begin{pmatrix}
	[P_{n + 1|N}]^{\frac{1}{2}}	&	0 \newline
	 0							&	[P_{n|N}]^{\frac{1}{2}} K_{n|N}^T
	 0							&	0
	 * \end{pmatrix}
	**/
	void tkalman_nc_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
										 gsl_matrix * c_s,
										 const gsl_matrix * sqrt_p_f,
										 const gsl_matrix * sqrt_p_s_,
										 const gsl_matrix * f_xx,
										 const gsl_matrix * sqrt_q_xx,
										 const gsl_matrix * s_gain,
										 gsl_matrix * mat_3x2x,
										 gsl_matrix * mat_3x2x_view_00,
										 gsl_matrix * mat_3x2x_view_01,
										 gsl_matrix * mat_3x2x_view_10,
										 gsl_matrix * mat_3x2x_view_11,
										 gsl_matrix * mat_3x2x_view_20,
										 gsl_matrix * mat_3x2x_view_21,
										 gsl_vector * vect_2x);
										 
	/**@fn void tkalman_nc_get_sqrt_q_0_s_and_c_s(gsl_matrix * sqrt_q_0_s,
												  gsl_matrix * c_s,
												  const gsl_matrix * sqrt_q_0_f,
												  const gsl_matrix * sqrt_p_1_s,
												  const gsl_matrix * f_xt,
												  const gsl_matrix * sqrt_q_xx,
												  const gsl_matrix * s_gain,
												  gsl_matrix * mat_2xpt_xpt,
												  gsl_matrix * mat_2xpt_xpt_view_00,
												  gsl_matrix * mat_2xpt_xpt_view_01,
												  gsl_matrix * mat_2xpt_xpt_view_10,
												  gsl_matrix * mat_2xpt_xpt_view_11,
												  gsl_matrix * mat_2xpt_xpt_view_20,
												  gsl_matrix * mat_2xpt_xpt_view_21,
												  gsl_vector * vect_xpt)
	 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
	 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
	 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
	 * @param[in] f_xt : F^{x,t}
	 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	 * @param[in] s_gain :  K_{n|N}, gain de Kalman triplet pour le lissage
	 * @param mat_2xpt_xpt : matrice de taille (2n_x + n_t, n_x + n_t) préallouée
	 * @param mat_2xpt_xpt_view_00 : vue sur la matrice mat_2xpt_xpt allant de (0,0) à (n_x - 1 , n_x - 1)
	 * @param mat_2xpt_xpt_view_01 : vue sur la matrice mat_2xpt_xpt allant de (0,n_x) à (n_x - 1, n_x + n_t - 1)
	 * @param mat_2xpt_xpt_view_10 : vue sur la matrice mat_2xpt_xpt allant de (n_x,0) à (n_x + n_t - 1 , n_x - 1)
	 * @param mat_2xpt_xpt_view_11 : vue sur la matrice mat_2xpt_xpt allant de (n_x,n_x) à (n_x + n_t - 1, n_x + n_t - 1)
	 * @param mat_2xpt_xpt_view_20 : vue sur la matrice mat_2xpt_xpt allant de (n_x + n_t,0) à (2n_x + n_t - 1, n_x - 1)
	 * @param mat_2xpt_xpt_view_21 : vue sur la matrice mat_2xpt_xpt allant de (n_x +n_t,n_x) à (2n_x + n_t - 1, n_x + n_t - 1)
	 * @param vect_xpt : vecteur de taille (n_x + n_t) préalloué
	 * @brief
	 Cette fonction calcule la racine de la matrice de covariance de l'état lissé courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
	 * On calcule la matrice
	 * \begin{pmatrix}
	 [Q^{x,x}]^{\frac{1}{2}}	&	0 \newline
	 [P_{n|n}]^{\frac{1}{2}}	&	[ F^{x,x}]^T \newline
	 0							&	[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
	 * \end{pmatrix}
	 * On effectue la décomposition QR
	 * \begin{pmatrix}
	[P_{n + 1|N}]^{\frac{1}{2}}	&	0 \newline
	 0							&	[P_{n|N}]^{\frac{1}{2}} K_{n|N}^T
	 0							&	0
	 * \end{pmatrix}
	**/
	void tkalman_nc_get_sqrt_q_0_s_and_c_s(gsl_matrix * sqrt_q_0_s,
										   gsl_matrix * c_s,
										   const gsl_matrix * sqrt_q_0_f,
										   const gsl_matrix * sqrt_p_1_s,
										   const gsl_matrix * f_xt,
										   const gsl_matrix * sqrt_q_xx,
										   const gsl_matrix * s_gain,
										   gsl_matrix * mat_2xpt_xpt,
										   gsl_matrix * mat_2xpt_xpt_view_00,
										   gsl_matrix * mat_2xpt_xpt_view_01,
										   gsl_matrix * mat_2xpt_xpt_view_10,
										   gsl_matrix * mat_2xpt_xpt_view_11,
										   gsl_matrix * mat_2xpt_xpt_view_20,
										   gsl_matrix * mat_2xpt_xpt_view_21,
										   gsl_vector * vect_xpt);
										   
	/**@fn void tkalman_nc_do_smoothing ( gsl_vector * x_s,
										  gsl_matrix * sqrt_p_s,
										  gsl_matrix * c_s,
										  const gsl_vector * x_f,
										  const gsl_matrix * sqrt_p_f,
										  const gsl_vector * x_p_,
										  const gsl_matrix * sqrt_p_p_,
										  const gsl_vector * x_s_,
										  const gsl_matrix * sqrt_p_s_,
										  const gsl_matrix * f_xx,
										  const gsl_matrix * sqrt_q_xx,
										  gsl_matrix * mat_xx,
										  gsl_matrix * mat_3x2x,
										  gsl_matrix * mat_3x2x_view_00,
										  gsl_matrix * mat_3x2x_view_01,
										  gsl_matrix * mat_3x2x_view_10,
										  gsl_matrix * mat_3x2x_view_11,
										  gsl_matrix * mat_3x2x_view_20,
										  gsl_matrix * mat_3x2x_view_21,
										  gsl_permutation * perm_x,
										  gsl_vector * vect_2x )
	 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
	 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
	 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
	 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
	 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
	 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
	 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
	 * @param[in] sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
	 * @param[in] f_xx : F^{x,x}
	 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
	 * @param mat_3x2x : matrice de taille (3n_x.2n_x) préallouée
	 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x - 1 , n_x - 1)
	 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x - 1, 2 n_x - 1)
	 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x - 1 , n_x - 1)
	 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x - 1, 2n_x - 1)
	 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x - 1, n_x - 1)
	 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x - 1, 2n_x - 1)
	 * @param perm_x : permutation de taille (n_x) préallouée
	 * @param vect_2x : vecteur de taille (2n_x) préalloué
	 * @brief
	 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_p_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_x_s .
	 */
	void tkalman_nc_do_smoothing ( gsl_vector * x_s,
								   gsl_matrix * sqrt_p_s,
								   gsl_matrix * c_s,
								   const gsl_vector * x_f,
								   const gsl_matrix * sqrt_p_f,
								   const gsl_vector * x_p_,
								   const gsl_matrix * sqrt_p_p_,
								   const gsl_vector * x_s_,
								   const gsl_matrix * sqrt_p_s_,
								   const gsl_matrix * f_xx,
								   const gsl_matrix * sqrt_q_xx,
								   gsl_matrix * mat_xx,
								   gsl_matrix * mat_3x2x,
								   gsl_matrix * mat_3x2x_view_00,
								   gsl_matrix * mat_3x2x_view_01,
								   gsl_matrix * mat_3x2x_view_10,
								   gsl_matrix * mat_3x2x_view_11,
								   gsl_matrix * mat_3x2x_view_20,
								   gsl_matrix * mat_3x2x_view_21,
								   gsl_permutation * perm_x,
								   gsl_vector * vect_2x );
								   
	/**@fn void tkalman_nc_do_smoothing_0 ( gsl_vector * t_0_s,
											gsl_matrix * sqrt_q_0_s,
											gsl_matrix * c_s,
											const gsl_vector * t_0_f,
											const gsl_matrix * sqrt_q_0_f,
											const gsl_vector * x_1_p,
											const gsl_matrix * sqrt_p_1_p,
											const gsl_vector * x_1_s,
											const gsl_matrix * sqrt_p_1_s,
											const gsl_matrix * f_xt,
											const gsl_matrix * sqrt_q_xx,
											gsl_matrix * mat_tx,
											gsl_matrix * mat_2xpt_xpt,
											gsl_matrix * mat_2xpt_xpt_view_00,
											gsl_matrix * mat_2xpt_xpt_view_01,
											gsl_matrix * mat_2xpt_xpt_view_10,
											gsl_matrix * mat_2xpt_xpt_view_11,
											gsl_matrix * mat_2xpt_xpt_view_20,
											gsl_matrix * mat_2xpt_xpt_view_21,
											gsl_permutation * perm_x,
											gsl_vector * vect_xpt)
	 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
	 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
	 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
	 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
	 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
	 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
	 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
	 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
	 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
	 * @param[in] f_xt : F^{x,t}
	 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	 * @param mat_tx : matrice de taille (n_t.n_x) préallouée
	 * @param mat_2xpt_xpt : matrice de taille (2n_x + n_t, n_x + n_t) préallouée
	 * @param mat_2xpt_xpt_view_00 : vue sur la matrice mat_2xpt_xpt allant de (0,0) à (n_x - 1 , n_x - 1)
	 * @param mat_2xpt_xpt_view_01 : vue sur la matrice mat_2xpt_xpt allant de (0,n_x) à (n_x - 1, n_x + n_t - 1)
	 * @param mat_2xpt_xpt_view_10 : vue sur la matrice mat_2xpt_xpt allant de (n_x,0) à (n_x + n_t - 1 , n_x - 1)
	 * @param mat_2xpt_xpt_view_11 : vue sur la matrice mat_2xpt_xpt allant de (n_x,n_x) à (n_x + n_t - 1, n_x + n_t - 1)
	 * @param mat_2xpt_xpt_view_20 : vue sur la matrice mat_2xpt_xpt allant de (n_x + n_t,0) à (2n_x + n_t - 1, n_x - 1)
	 * @param mat_2xpt_xpt_view_21 : vue sur la matrice mat_2xpt_xpt allant de (n_x +n_t,n_x) à (2n_x + n_t - 1, n_x + n_t - 1)
	 * @param vect_xpt : vecteur de taille (n_x + n_t) préalloué
	 * @param perm_x : permutation de taille (n_x) préallouée
	 * @brief
	 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain_0. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_q_0_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_t_0_s .
	 */
	void tkalman_nc_do_smoothing_0 ( gsl_vector * t_0_s,
									 gsl_matrix * sqrt_q_0_s,
									 gsl_matrix * c_s,
									 const gsl_vector * t_0_f,
									 const gsl_matrix * sqrt_q_0_f,
									 const gsl_vector * x_1_p,
									 const gsl_matrix * sqrt_p_1_p,
									 const gsl_vector * x_1_s,
									 const gsl_matrix * sqrt_p_1_s,
									 const gsl_matrix * f_xt,
									 const gsl_matrix * sqrt_q_xx,
									 gsl_matrix * mat_tx,
									 gsl_matrix * mat_2xpt_xpt,
									 gsl_matrix * mat_2xpt_xpt_view_00,
									 gsl_matrix * mat_2xpt_xpt_view_01,
									 gsl_matrix * mat_2xpt_xpt_view_10,
									 gsl_matrix * mat_2xpt_xpt_view_11,
									 gsl_matrix * mat_2xpt_xpt_view_20,
									 gsl_matrix * mat_2xpt_xpt_view_21,
									 gsl_permutation * perm_x,
									 gsl_vector * vect_xpt);
		
#endif
