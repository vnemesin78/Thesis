/**@file tkalman_nc_sums.hpp
 * @author Valérian Némesin
 */
#ifndef _TKALMAN_NC_SUMS_BASE_HPP_
	#define _TKALMAN_NC_SUMS_BASE_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <iostream> //A virer
	using namespace std;
	/**@class tkalman_nc_sums_base
	 * @brief
	 * Cette classe permet de calculer les sommes nécessaires à l'algorithme EM du filtre de Kalman couple non-corrélé.
	 */
	class tkalman_nc_sums_base
	{
		public:
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_sums
			 */
			tkalman_nc_sums_base( const gsl_matrix * f_xt,
								  const gsl_matrix * sqrt_q_xx) throw (exception &);
			
			/**@fn
			 * @param[in] f_xt : F^{x,t}.
			 * @param[in] sqrt_q_xx, [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process
			 * @brief
			 * Setup de la classe @class tkalman_nc_sums
			 */
			virtual void setup ( const gsl_matrix * f_xt,
								 const gsl_matrix * sqrt_q_xx) throw (exception &); 
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~tkalman_nc_sums_base();
		

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
			const gsl_matrix * _f_xt;
			const gsl_matrix * _sqrt_q_xx;
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			
			//Tmp
			gsl_matrix * mat_2xpt_2xpt,
						 mat_2xpt_2xpt_view_00,
						 mat_2xpt_2xpt_view_10,
						 mat_2xpt_2xpt_view_11,
						 mat_2xpt_2xpt_view_12,
					     mat_2xpt_2xpt_view_21,
						 mat_2xpt_2xpt_view_22;
			gsl_vector * vect_2xpt ;
			gsl_matrix * mat_3x_3x,
						 mat_3x_3x_view_00,
						 mat_3x_3x_view_10,
						 mat_3x_3x_view_11,
						 mat_3x_3x_view_12,
					     mat_3x_3x_view_21,
						 mat_3x_3x_view_22;
			gsl_vector * vect_3x ;
			
			
			gsl_matrix * mat_2tp2xp1_tpx,
					     mat_2tp2xp1_tpx_view_00,
					     mat_2tp2xp1_tpx_view_tpxp1_tpx,
					     mat_tpxp1_tpx_view_00,
					     mat_tpxp1_tpx_view_00_bis,
					     mat_tpxp1_tpx_view_02, 
					     mat_tpxp1_tpx_view_02_bis, 
					     mat_tpxp1_tpx_view_22;
			gsl_vector mat_tpxp1_tpx_view_30,
					   mat_tpxp1_tpx_view_31,
					   mat_tpxp1_tpx_view_32;
			gsl_vector * vect_tpx;
			
			gsl_matrix * mat_2tp2yp1_tpy,
					     mat_2tp2yp1_tpy_view_00,
					     mat_2tp2yp1_tpy_view_tpyp1_tpy,
					     mat_tpyp1_tpy_view_00,
					     mat_tpyp1_tpy_view_00_bis;
			gsl_vector mat_tpyp1_tpy_view_30,
					   mat_tpyp1_tpy_view_31,
					   mat_tpyp1_tpy_view_32;
			gsl_vector * vect_tpy;
			
			
		//Accesseur
		public:
			/**@fn inline const gsl_matrix * tkalman_nc_sums :: f_xt() const
			 * @return
			 * F^{x,t}
			 **/
			inline const gsl_matrix * f_xt() const
			{
				return _f_xt;
			}
			
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
		
			
			
						 
			
	};
	
	
	
	
	
	
	
	
	
	
	
	/**@fn void tkalman_nc_get_cov_zx_(const gsl_matrix * sqrt_z_f,
									   const gsl_matrix * sqrt_p_s_,
									   const gsl_matrix * c_s,
									   const gsl_matrix * f_xz,
									   const gsl_matrix * sqrt_q_xx,
									   gsl_matrix * mat_2xpz_2xpz,
									   gsl_matrix * mat_2xpz_2xpz_view_00,
									   gsl_matrix * mat_2xpz_2xpz_view_10,
									   gsl_matrix * mat_2xpz_2xpz_view_11,
									   gsl_matrix * mat_2xpz_2xpz_view_21,
									   gsl_matrix * mat_2xpz_2xpz_view_22,
									   gsl_vector * vect_2xpz)
	 * @param[in] sqrt_z_f : [Z_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré courant.
	 * @param[in] sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant.
	 * @param[in] c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
	 * @param[in] f_xz : F^{x,z}
	 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	 * @param mat_2xpz_2xpz : matrice de taille (2n_x + n_z, 2 n_x + n_z) allouée.
	 * @param mat_2xpz_2xpz_view_00 : vue sur la matrice mat_2xpz_2xpz allant de (0,0) à (n_x - 1, n_x - 1).
	 * @param mat_2xpz_2xpz_view_10 : vue sur la matrice mat_2xpz_2xpz allant de (n_x, 0) à (n_z + n_x - 1, n_x - 1).
	 * @param mat_2xpz_2xpz_view_11 : vue sur la matrice mat_2xpz_2xpz allant de (n_x, n_x) à (n_z + n_x - 1, n_z + n_x - 1).
	 * @param mat_2xpz_2xpz_view_21 : vue sur la matrice mat_2xpz_2xpz allant de (n_z + n_x, n_x) à (n_z + 2 n_x - 1, n_z + n_x - 1).
	 * @param mat_2xpz_2xpz_view_22 : vue sur la matrice mat_2xpz_2xpz allant de (n_z + n_x, n_z + n_x) à (n_z + 2 n_x - 1, n_z + 2 n_x - 1).
	 * @param vect_2xpz : vecteur de taille (2n_x + n_z) alloué
	 * @brief
	 * Cette fonction calcule la racine de la matrice de covariance du vecteur aléatoire (z_{n|N}; x_{n + 1 |N}). (z = x ou t). Pour réaliser ce calcul, on réalise la décomposition QR de la matrice suivante : 
	 * M = 
	 * \begin{pmatrix}
	 * \left[Q_2^{x,x}\right]^{\frac{1}{2}}							&	0_{n_x, n_t}	&	0_{n_x, n_x}	\newline
	 * \left[Z_{n|n}\right]^{\frac{1}{2}} \left[F_2^{x} \right]^T	&	\left[Z_{n|n}\right]^{\frac{1}{2}}	&	0_{n_t, n_x}	\newline
	 * 0_{n_x, n_x}													&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}	&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}
	 * \end{pmatrix}
	 * Elle s'écrit : 
	 * M =
	 * Q
	 * \begin{pmatrix}
	 * 	M_{0,0}	&	M_{0,1}	&	M_{0,2} \newline
	 * 	O	    	&		M_{1,1}	&	M_{1,2} \newline 
	 * 	O   	 	&		O	       &	M_{2,2}
	 * \end{pmatrix}
	 *  
	 * La racine recherchée est la matrice
	 * \left[Z_{n|N}\right]^{\frac{1}{2}} = 
	 * \begin{pmatrix}
	 * 	M_{1,1}	&	M_{1,2} \newline 
	 * 	0    	&	M_{2,2}
	 * \end{pmatrix}
	 *
	 * Pour pouvoir utiliser cette fonction, il est nécessaire de construire avant une vue sur la matrice M allant de (n_x + n_z, n_x + n_z) à (2 n_x + n_z - 1, 2 n_x + n_z - 1). (Le résultat est stocké dedans.). 
	 */
	void tkalman_nc_get_cov_zx_ ( const gsl_matrix * sqrt_z_f,
								  const gsl_matrix * sqrt_p_s_,
								  const gsl_matrix * c_s,
								  const gsl_matrix * f_xz,
								  const gsl_matrix * sqrt_q_xx,
								  gsl_matrix * mat_2xpz_2xpz,
								  gsl_matrix * mat_2xpz_2xpz_view_00,
								  gsl_matrix * mat_2xpz_2xpz_view_10,
								  gsl_matrix * mat_2xpz_2xpz_view_11,
								  gsl_matrix * mat_2xpz_2xpz_view_21,
								  gsl_matrix * mat_2xpz_2xpz_view_22,
								  gsl_vector * vect_2xpz             );
	/**@fn void tkalman_nc_get_corr_tx_ ( const gsl_vector * x_s,
										  const gsl_vector * _y,
										  const gsl_vector * x_s_,
										  const gsl_matrix * sqrt_cov_xx_view_00,
										  const gsl_matrix * sqrt_cov_xx_view_01,
										  const gsl_matrix * sqrt_cov_xx_view_11,
										  gsl_matrix * mat_tpxp1_tpx,
										  gsl_matrix * mat_tpxp1_tpx_view_00, // Modif
										  gsl_matrix * mat_tpxp1_tpx_view_02, 
										  gsl_matrix * mat_tpxp1_tpx_view_22, // Modif
										  gsl_vector * mat_tpxp1_tpx_view_30,
										  gsl_vector * mat_tpxp1_tpx_view_31,
										  gsl_vector * mat_tpxp1_tpx_view_32,
										  gsl_vector * vect_tpx )
	 * @param[in] x_s : \hat{x}_{n|N}, espérance de l'état lissé courant
	 * @param[in] _y : y_{n - 1}, observation précédente
	 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
	 * @param[in] sqrt_cov_xx_view_00 : [Z_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant.
	 * @param[in] sqrt_cov_xx_view_01 : [[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{0,1}, vue sur la racine de la matrice de covariance du vecteur aléatoire (x_{n|N}^T, x_{n + 1|N}^T,)^T partant de (0, n_x) à (n_x - 1, 2 n_x - 1).
	 * @param[in] sqrt_cov_xx_view_11 : [[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{1,1}, vue sur la racine de la matrice de covariance du vecteur aléatoire (x_{n|N}^T, x_{n + 1|N}^T,)^T partant de (n_x, n_x) à (2 n_x - 1, 2 n_x - 1).
	 * @param mat_tpxp1_tpx : matrice de taille (n_t + n_x + 1, n_t + n_x) allouée
	 * @param mat_tpxp1_tpx_view_00 : vue sur la matrice mat_tpxp1_tpx partant de (0, 0) à (n_x - 1, n_x - 1).
	 * @param mat_tpxp1_tpx_view_02 : vue sur la matrice mat_tpxp1_tpx partant de (0, n_t) à (n_t + n_x - 1, n_x - 1).
	 * @param mat_tpxp1_tpx_view_22 : vue sur la matrice mat_tpxp1_tpx partant de (n_t, n_t) à (n_t + n_x - 1, n_t + n_x - 1).
	 * @param mat_tpxp1_tpx_view_30 : vue sur la matrice mat_tpxp1_tpx partant de (2 n_t, 0) à (2 n_t, n_x - 1).
	 * @param mat_tpxp1_tpx_view_31 : vue sur la matrice mat_tpxp1_tpx partant de (2 n_t, n_x) à (2 n_t, n_t - 1).
	 * @param mat_tpxp1_tpx_view_32 : vue sur la matrice mat_tpxp1_tpx partant de (2 n_t, n_t) à (2 n_t, n_t + n_x - 1).
	 * @param vect_tpx : vecteur de taille (n_t + n_x) alloué
	 * @brief
	 * Cette fonction calcule la racine de la matrice de corrélation du vecteur (t_{n|N}; x_{n + 1|N}). \n
	 * Le résultat sera stocké dans la matrice mat_tpxp1_tpx. Pour le récupérer, il suffit de créer une vue partant de (0, 0) à (n_t + n_x - 1, n_t + n_x - 1)
	 * Ce calcul s'effectue en trois étapes :
	 * - On reconstruit dans un premier temps la racine de la matrice de covariance du vecteur (t_{n|N}; x_{n + 1|N}) :
	 * $$
	 * M =  
	 * \begin{pmatrix}
	 * [P_{n|N}]^{\frac{1}{2}}	&	0	&	[[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{0,1}	\newline
	 * 0						&	0	&	0														 	\newline
	 * 0						&	0	&	[[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{1,1}	\newline
	 * 0						&	0	&	0														
	 * \end{pmatrix}
	 * $$
	 * - Puis dans un second temps, on complète cette matrice M : \n
	 * $$
	 * N = 
	 * \begin{pmatrix}
	 * [P_{n|N}]^{\frac{1}{2}}	&	0			&	[[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{0,1}		\newline
	 * 0						&	0			&	0																	\newline
	 * 0						&	0			&	[[cov((x_{n|N}^T, x_{n + 1|N}^T,)^T)]^{\frac{1}{2}}]^{1,1}		\newline
	 * \hat{x}_{n|N}^T			&	y_{n - 1}^T	&	\hat{x}_{n + 1|N}^T											
	 * \end{pmatrix}
	 * $$
	 * - Ensuite, nous effectuons la décomposition QR de cette matrice qui nous permet d'obtenir la matrice de corrélation recherchée : \n
	 * $$
	 * N = Q 
	 * \begin{pmatrix}
	 * [corr(t_{n|N}; x_{n + 1|N})]^{\frac{1}{2}} \newline
	 * 0
	 * \end{pmatrix}
	 * $$
	 */
	void tkalman_nc_get_corr_tx_ ( const gsl_vector * x_s,
								   const gsl_vector * _y,
								   const gsl_vector * x_s_,
								   const gsl_matrix * sqrt_cov_xx_view_00,
								   const gsl_matrix * sqrt_cov_xx_view_01,
								   const gsl_matrix * sqrt_cov_xx_view_11,
								   gsl_matrix * mat_tpxp1_tpx,
								   gsl_matrix * mat_tpxp1_tpx_view_00,
								   gsl_matrix * mat_tpxp1_tpx_view_02, 
								   gsl_matrix * mat_tpxp1_tpx_view_22,
								   gsl_vector * mat_tpxp1_tpx_view_30,
								   gsl_vector * mat_tpxp1_tpx_view_31,
								   gsl_vector * mat_tpxp1_tpx_view_32,
								   gsl_vector * vect_tpx );
	/**@fn void tkalman_nc_get_corr_ty ( const gsl_vector * x_s,
										 const gsl_vector * _y,
										 const gsl_vector * y,
										 const gsl_matrix * sqrt_z_s,
										 gsl_matrix * mat_tpyp1_tpy,
										 gsl_matrix * mat_tpyp1_tpy_view_00,
										 gsl_vector * mat_tpyp1_tpy_view_30,
										 gsl_vector * mat_tpyp1_tpy_view_31,
										 gsl_vector * mat_tpyp1_tpy_view_32,
										 gsl_vector * vect_tpy )
	 * @param[in] x_s : \hat{x}_{n|N}, espérance de l'état lissé courant
	 * @param[in] _y : y_{n - 1}, observation précédente
	 * @param[in] y :  y_{n}, observation courante
	 * @param[in] sqrt_z_s : [Z_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant.
	 * @param mat_tpyp1_tpy : matrice de taille (n_t + n_y + 1, n_t + n_y) allouée
	 * @param mat_tpyp1_tpy_view_00 : vue sur la matrice mat_tpyp1_tpy partant de (0, 0) à (n_z - 1, n_z - 1).).
	 * @param mat_tpyp1_tpy_view_30 : vue sur la matrice mat_tpyp1_tpy partant de (2 n_t, 0) à (2 n_t, n_x - 1).
	 * @param mat_tpyp1_tpy_view_31 : vue sur la matrice mat_tpyp1_tpy partant de (2 n_t, n_x) à (2 n_t, n_t - 1).
	 * @param mat_tpyp1_tpy_view_32 : vue sur la matrice mat_tpyp1_tpy partant de (2 n_t, n_t) à (2 n_t, n_t + n_y - 1).
	 * @param vect_tpy : vecteur de taille (n_t + n_y) alloué
	 * @brief
	 * Cette fonction calcule la racine de la matrice de corrélation du vecteur (t_{n|N}; y_{n|N}). \n
	 * Le résultat sera stocké dans la matrice mat_tpxp1_tpx. Pour le récupérer, il suffit de créer une vue partant de (0, 0) à (n_t + n_y - 1, n_t + n_y - 1)
	 * Ce calcul s'effectue en trois étapes :
	 * - On reconstruit dans un premier temps la racine de la matrice de covariance du vecteur (t_{n|N}; y_{n|N}) :
	 * $$
	 * M =  
	 * \begin{pmatrix}
	 * [P_{n|N}]^{\frac{1}{2}}	&	0	&	0	\newline
	 * 0						&	0	&	0	\newline
	 * 0						&	0	&   0   \newline
	 * 0						&	0	&	0	
	 * \end{pmatrix}
	 * $$
	 * - Puis dans un second temps, on complète cette matrice M : \n
	 * $$
	 * N = 
	 * \begin{pmatrix}
	 * [P_{n|N}]^{\frac{1}{2}}	&	0			&	0	\newline
	 * 0						&	0			&	0	\newline
	 * 0						&	0			&   0   \newline
	 * \hat{x}_{n|N}^T			&	y_{n - 1}^T	&	y_{n}^T										
	 * \end{pmatrix}
	 * $$
	 * - Ensuite, nous effectuons la décomposition QR de cette matrice qui nous permet d'obtenir la matrice de corrélation recherchée : \n
	 * $$
	 * N = Q 
	 * \begin{pmatrix}
	 * [corr(t_{n|N}; y_{n|N})]^{\frac{1}{2}} \newline
	 * 0
	 * \end{pmatrix}
	 * $$
	 */
	void tkalman_nc_get_corr_ty ( const gsl_vector * x_s,
								  const gsl_vector * _y,
								  const gsl_vector * y,
								  const gsl_matrix * sqrt_z_s,
								  gsl_matrix * mat_tpyp1_tpy,
								  gsl_matrix * mat_tpyp1_tpy_view_00,
								  gsl_vector * mat_tpyp1_tpy_view_30,
								  gsl_vector * mat_tpyp1_tpy_view_31,
								  gsl_vector * mat_tpyp1_tpy_view_32,
								  gsl_vector * vect_tpy
								 );
		
	
#endif

