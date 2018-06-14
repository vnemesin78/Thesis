/**@file c_iris_segmentation.hpp
 *
 **/
#ifndef _C_IRIS_SEGMENTATION_HPP_
	#define _C_IRIS_SEGMENTATION_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <cmath>
	#include "lib_image.hpp"
	#include "lib_api.hpp"
	#include <iostream>
	#include <sstream>
	/**@class c_iris_segmentation
	 * @brief
	 * Cette classe permet la segmentation de l'iris.
	 */
	class c_iris_segmentation
	{
		public:
			/**@fn
			 * @brief
			 * Constructor
			 **/
			c_iris_segmentation ();

			/**@fn
			 * @param width : largeur
			 * @param height : hauteur
			 * @param n_r_min : Ratio minimum entre pupille et iris
			 * @param n_r_max : Ratio maximum entre pupille et iris
			 * @param n_r_Daugman : Ratio pour l'algo de Mouloud
			 * @param dx, dy, dr : zone de recherche pour l'algo de Daugman
			 * @param sigma_iris, nu_bg, sigma_bg : moments statistique du background et de l'iris
			 * @brief
			 * Constructor
			 **/
			c_iris_segmentation (	unsigned int width, //Largeur
									unsigned int height, //Hauteur
									unsigned int nb_directions,
									unsigned int kernel_size,
									unsigned int nb_samples,
									unsigned int dx,
									unsigned int nb_x,
									unsigned int dy,
									unsigned int nb_y,
									unsigned int r_min,
									unsigned int r_max,
									double r_ratio,
									unsigned int nb_r_ellipse,
									unsigned int dx_ellipse,
									unsigned int nb_x_ellipse,
									unsigned int dy_ellipse,
									unsigned int nb_y_ellipse,
									unsigned int nb_thetas,
									ostream * _err_stream = NULL
									 );
			int default_setup( 	unsigned int width, //Largeur
								unsigned int height,
								ostream * _err_stream = NULL );
			/**@fn
			 * @param width : largeur
			 * @param height : hauteur
			 * @param n_r_min : Ratio minimum entre pupille et iris
			 * @param n_r_max : Ratio maximum entre pupille et iris
			 * @param n_r_Daugman : Ratio pour l'algo de Mouloud
			 * @param dx, dy, dr : zone de recherche pour l'algo de Daugman
			 * @param sigma_iris, nu_bg, sigma_bg : moments statistique du background et de l'iris
			 * @brief
			 * Setup
			 * @return
			 * - 1 en cas de prob.
			 **/
			int setup (	unsigned int width, //Largeur
						unsigned int height, //Hauteur
						unsigned int nb_directions,
						unsigned int kernel_size,
						unsigned int nb_samples,
						unsigned int dx,
						unsigned int nb_x,
						unsigned int dy,
						unsigned int nb_y,
						unsigned int r_min,
						unsigned int r_max,
						double r_ratio,
						unsigned int nb_r_ellipse,
						unsigned int dx_ellipse,
						unsigned int nb_x_ellipse,
						unsigned int dy_ellipse,
						unsigned int nb_y_ellipse,
						unsigned int nb_thetas,
						ostream * _err_stream = NULL
						 );

			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param prefix : prefix devant les noms de variabels
			 * @param *_name : noms des différentes variables
			 * @return
			 * - 0 si tout se déroule bien
			 * - autre sinon
			 * @brief
			 * Setup
			 **/
			int setup( 	api_parameters & params,
						unsigned int width,
						unsigned int height,
						ostream * _err_stream = NULL,
						const char * n_space = "iris",
						const char * nb_directions_name = "nb_directions",
						const char * kernel_size_name = "kernel_size",
						const char * nb_samples_name = "nb_samples",
						const char * dx = "dx",
						const char * nb_x = "nb_x",
						const char * dy = "dy",
						const char * nb_y = "nb_y",
						const char * r_min = "r_min",
						const char * r_max = "r_max",
						const char * r_ratio = "r_ratio",
						const char * nb_r_ellipse = "nb_samples_ellipse",
						const char * dx_ellipse = "dx_ellipse",
						const char * nb_x_ellipse = "nb_x_ellipse",
						const char * dy_ellipse = "dy_ellipse",
						const char * nb_y_ellipse = "nb_y_ellipse",
						const char * nb_thetas	= "nb_thetas"
					);




		/**@fn
		 * @param image : image prétraitée
		 * @param mask : masque des spots ici
		 * @param x_p, y_p : centre de la pupille
		 * @param r_p : rayon de la pupille : sqrt(a² + b²)
		 * @param nu_iris : mode de l'iris
		 * @brief
		 * Segmentation de l'iris.
		 */
		int segment_iris (	const IplImage * image,
							const IplImage * preprocessed_image,
							const IplImage * mask,
							const IplImage * s_mask,
							const double & x_p,
							const double & y_p,
							const double & a_p,
							const double & b_p,
							const double & theta_p );

		int segment_iris (	const IplImage * image,
							const IplImage * mask,
							const double & x_p,
							const double & y_p,
							double r_p = 0 );
		/**@fn
		 * @brief
		 * Destructor
		 *
		**/
		~c_iris_segmentation();
	//Méthodes
	protected:
		/**@fn
		 * @brief
		 * Cette méthode libère la mémoire allouée pour les attributs.
		 **/
		void free();
		/**@fn
		 * @brief
		 * Cette méthode met les attributs à 0.
		 **/
		void initialize();
		/**@fn
		 * @brief
		 * Alloc. mémoire.
		 *
		 */
		void alloc();

	//Paramètres
	protected:
		unsigned int 	_width,
						_height;
		double _spot_threshold;
	public:
		inline unsigned int width() const
		{
			return _width;
		}
		inline unsigned int height() const
		{
			return _height;
		}
	protected:
		unsigned int _nb_points;
	public:
		inline unsigned int nb_points() const
		{
			return _nb_points;
		}

	protected:
		unsigned int _kernel_size;
	public:
		inline unsigned int kernel_size() const
		{
			return _kernel_size;
		}




	//Paramètre cercle
	protected:
		unsigned int _dx;
	public:
		inline unsigned int dx() const
		{
			return _dx;
		}

	protected:
		unsigned int _nb_x;
	public:
		inline unsigned int nb_x() const
		{
			return _nb_x;
		}

	protected:
		unsigned int _dy;
	public:
		inline unsigned int dy() const
		{
			return _dy;
		}

	protected:
		unsigned int _nb_y;
	public:
		inline unsigned int nb_y() const
		{
			return _nb_y;
		}

	protected:
		unsigned int _r_min;
	public:
		inline unsigned int r_min() const
		{
			return _r_min;
		}

	protected:
		unsigned int _r_max;
	public:
		inline unsigned int r_max() const
		{
			return _r_max;
		}

	protected:
		unsigned int _nb_r;
	public:
		inline unsigned int nb_r() const
		{
			return _nb_r;
		}


	//Paramètres ellipse
	protected:
		double _r_ratio;
	public:
		inline double r_ratio() const
		{
			return _r_ratio;
		}

	protected:
		unsigned int _nb_r_ellipse;
	public:
		inline unsigned int nb_r_ellipse() const
		{
			return _nb_r_ellipse;
		}

	protected:
		unsigned int _dx_ellipse;
	public:
		inline unsigned int dx_ellipse() const
		{
			return _dx_ellipse;
		}

	protected:
		unsigned int _nb_x_ellipse;
	public:
		inline unsigned int nb_x_ellipse() const
		{
			return _nb_x_ellipse;
		}

	protected:
		unsigned int _dy_ellipse;
	public:
		inline unsigned int dy_ellipse() const
		{
			return _dy_ellipse;
		}

	protected:
		unsigned int _nb_y_ellipse;
	public:
		inline unsigned int nb_y_ellipse() const
		{
			return _nb_y_ellipse;
		}
	protected:
		unsigned int _nb_thetas;
	public:
		inline unsigned int nb_thetas() const
		{
			return _nb_thetas;
		}

	//Segmenation
	protected:
		c_integro_differential_polar integro_diff_ellipse_2;
	//Données
		//Segmentation par Daugman
		protected:
			double _x_Daugman,
					_y_Daugman,
					_r_Daugman;
		public:
			inline double x_Daugman() const
			{
				return _x_Daugman;
			}

			inline double y_Daugman() const
			{
				return _y_Daugman;
			}

			inline double r_Daugman() const
			{
				return _r_Daugman;
			}

			inline double x() const
			{
				return _x_Daugman;
			}

			inline double y() const
			{
				return _y_Daugman;
			}

			inline double r() const
			{
				return _r_Daugman;
			}






		//Amélioration elliptique
		protected:
			double _x_ellipse,
					_y_ellipse,
					_a_ellipse,
					_b_ellipse,
					_theta_ellipse;
		public:
			inline double x_ellipse() const
			{
				return _x_ellipse;
			}

			inline double y_ellipse() const
			{
				return _y_ellipse;
			}

			inline double a_ellipse() const
			{
				return _a_ellipse;
			}

			inline double b_ellipse() const
			{
				return _b_ellipse;
			}

			inline double theta_ellipse() const
			{
				return _theta_ellipse;
			}
		protected:
			IplImage * _new_image,
					 * _new_mask;

			double _new_x,
					_new_y,
					_new_r;

			double _new_x_pupil,
					_new_y_pupil,
					_new_a_pupil,
					_new_b_pupil,
					_new_theta_pupil;

		public:
			inline IplImage * new_image ()
			{
				return _new_image;
			}
			inline const IplImage * new_mask () const
			{
				return _new_mask;
			}
			inline const double & new_x() const
			{
				return _new_x;
			}
			inline const double & new_y() const
			{
				return _new_y;
			}
			inline const double & new_r() const
			{
				return _new_r;
			}
			inline const double & new_x_pupil() const
			{
				return _new_x_pupil;
			}
			inline const double & new_y_pupil() const
			{
				return _new_y_pupil;
			}


			inline const double & new_a_pupil() const
			{
				return _new_a_pupil;
			}
			inline const double & new_b_pupil() const
			{
				return _new_b_pupil;
			}
			inline const double & new_theta_pupil() const
			{
				return _new_theta_pupil;
			}

		protected:
			double * pupil_points;
			fitting_ellipse ellipse;
		//Flux d'erreur
		protected:
			ostream * err_stream;
		public:
			/**@fn
			 * @param error_str : flux d'erreur
			 * @brief
			 * Modifie le flux d'erreur.
			 *
			 */
			inline void set_error_stream( ostream & error_str = cout )
			{
				err_stream = &error_str;
			}

	};



#endif
