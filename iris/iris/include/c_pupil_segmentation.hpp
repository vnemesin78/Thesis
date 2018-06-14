/**@file c_pupil_segmentation.hpp
 * @author 
 * Valérian Némesin
 * @brief
 * Ce fichier contient la classe qui permet de segmenter la pupille dans une image.
 * 
 **/
#ifndef _C_PUPIL_SEGMENTATION_HPP_
	#define _C_PUPIL_SEGMENTATION_HPP_
	#include "lib_image.hpp"
	#include "lib_api.hpp"
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <sstream>
	#include <iostream>
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <gsl/gsl_rng.h>
	using namespace std;
	/**@class
	 * @brief
	 * Cette classe gère la segmentation de la pupille.
	 * 
	 * 
	 **/
	class c_pupil_segmentation
	{
		public:
			
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 **/
			c_pupil_segmentation ( void );
			
		
			/**@fn
			 * @param width : largeur d'une image
			 * @param height : hauteur d'une image
			 * @param sigma : sigma du noyau gaussien qui filtre l'hist.
			 * @param nb_min_max : nombre maximum de minima sur l'histogramme
			 * @param radius_min : rayon min
			 * @param radius_max : rayon max
			 * @param diff_surf_max : différence maximale tolérée entre l'ellipse fitté et l'objet connexe
			 * @param nb_points_pupil : nombre de points du contour de la pupil pris en compte
			 * @param nb_iterations_ellipse : nombre d'itération de l'algorithme de fitting d'ellipse
			 * @brief
			 * Constructeur.
			 * 
			 **/
			c_pupil_segmentation(	unsigned int width, //Largeur
									unsigned int height, //Hauteur
									unsigned int window_size,
									double sigma_hist_image,
									double sigma_hist_diff_image,
									unsigned int pupil_limit_threshold,
									double radius_min,
									double radius_max,
									unsigned int nb_points, // Nombre de points du contour rééchantillonné de la pupille
									unsigned int nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
									unsigned int nb_points_ransac,
									unsigned int nb_iter_ransac,
									double ransac_threshold,
									ostream * _err_stream = NULL );
			
			int default_setup( 	unsigned int width,
								unsigned int height,
								ostream * _err_stream = NULL );
			
			/**@fn
			 * @param width : largeur d'une image
			 * @param height : hauteur d'une image
			 * @param sigma : sigma du noyau gaussien qui filtre l'hist.
			 * @param nb_min_max : nombre maximum de minima sur l'histogramme
			 * @param radius_min : rayon min
			 * @param radius_max : rayon max
			 * @param diff_surf_max : différence maximale tolérée entre l'ellipse fitté et l'objet connexe
			 * @param nb_points_pupil : nombre de points du contour de la pupil pris en compte
			 * @param nb_iterations_ellipse : nombre d'itération de l'algorithme de fitting d'ellipse
			 * @brief
			 * Setup
			 * @return
			 * - 1 en cas de problème
			 **/
			int setup (	unsigned int width, //Largeur
						unsigned int height, //Hauteur
						unsigned int window_size,
						double sigma_hist_image,
						double sigma_hist_diff_image,
						unsigned int pupil_limit_threshold,
						double radius_min,
						double radius_max,
						unsigned int nb_points, // Nombre de points du contour rééchantillonné de la pupille
						unsigned int nb_iterations_ellipse, //Nombre d'itérations de l'algorithme de fitting d'ellipse 
						unsigned int nb_points_ransac,
						unsigned int nb_iter_ransac,
						double ransac_threshold,
						ostream * _err_stream = NULL );
		
			/**@fn 
			 * @param params : variables en mémoire
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param n_space : nom de l'espace de nom des variables de l'objet
			 * @param *_name : nom des variables associées aux paramètres
			 * @brief
			 * Setup
			 * @return
			 * - 1 si problème sur les valeurs des paramètres
			 * - -1 si problème sur les variables
			 * 
			 * 
			 **/
			int  setup (	api_parameters & params,
							unsigned int width,
							unsigned int height,
							ostream * _err_stream = NULL,
							const char * n_space = "pupil",
							const char * window_size_name = "filter_window_size",
							const char * sigma_image_name = "sigma_image",
							const char * sigma_diff_name = "sigma_diff",
							const char * pupil_limit_threshold_name = "limit_threshold",
							const char * radius_min_name = "radius_min",
							const char * radius_max_name = "radius_max",
							const char * nb_points_name = "nb_points",
							const char * nb_iter_ellipse_name = "nb_iter_ellipse",
							const char * nb_points_ransac_name = "nb_points_ransac",
							const char * nb_iter_ransac_name = "nb_iter_ransac",
							const char * ransac_threshold_name = "ransac_threshold"
							);
		
			/**@fn
			 * @param img_data : pixels de l'image
			 * @param x : abscisse du début du roi sur l'image
			 * @param y : ordonnée du début du roi sur l'image
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param width_step : largeur réelle d'une ligne de l'image
			 * @return
			 * - 0 si la pupille a été segmentée
			 * - 1 sinon
			 * - -1 en cas de problème sur la taille de l'image allouée
			 * @brief
			 * Segmentation de la pupille
			 */
			int segment(	const unsigned char * img_data,
							unsigned int x,
							unsigned int y, 
							unsigned int width,
							unsigned int height,
							unsigned int width_step);
		
			
			/**@fn
			 * @param image : image à segmenter
			 * @return
			 * - 0 si la pupille a été segmentée
			 * - 1 sinon
			 * - -1 en cas de problème sur la taille de l'image allouée
			 * @brief
			 * Segmentation de la pupille
			 */
			int segment(	const IplImage * image );
		
		
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_pupil_segmentation();
		
		protected:
		
			//Images temporaires
			unsigned int _width, 
						 _height;
			IplImage * _abs_diff_image_8b;
		
			//Taille de la fenêtre pour le filtre moyen.
			double _window_size;
		public:
			inline double window_size() const
			{
				return _window_size;
			}
			
		protected:
			/**
			 * 
			 */
			void get_abs_diff_image( const unsigned char * img_data,
									 unsigned int x,
									 unsigned int y, 
									 unsigned int width,
									 unsigned int height,
									 unsigned int width_step );
					
			double _sigma_hist_1;
			double _sigma_hist_2;

		public:
			inline double sigma_image_hist( ) const
			{
				return _sigma_hist_1;
			}
			
			inline double sigma_d_image_hist( ) const
			{
				return _sigma_hist_2;
			}
		protected:
			c_histogram hist_1, hist_2;
			void compute_histograms( const unsigned char * img_data,
									 unsigned int x,
									 unsigned int y, 
									 unsigned int width,
									 unsigned int height,
									 unsigned int width_step );
		
			/**
			 * 
			 */
			unsigned int 	nb_modes_image,
							nb_modes_d_image;
			unsigned int * m_image,
						 * m_d_image;
			unsigned int _pupil_limit_threshold;
						  
			void mode_detection(	unsigned int width,
									unsigned int height );
			unsigned int * labels_1,
						 * labels_2;
			unsigned int _nb_labels;
		
			void get_labels( );
		
		
		
			IplImage * mask_1, * mask_2;
			void mask( 	 const unsigned char * img_data,
						 unsigned int x,
						 unsigned int y, 
						 unsigned int width,
						 unsigned int height,
						 unsigned int width_step,
						 unsigned int s1, 
						 unsigned int s2  );


			c_label * bg_regions,
					* connex_regions;
			int label( );
		
			unsigned int * _valid_region_labels;
			unsigned int _nb_valid_regions;
			
			
			double _radius_min;
			double _radius_max;
			int erase_small_regions();
		
			c_contour * contour_obj;
			double * _contour;
			unsigned int _nb_points_contour;
			unsigned char * _contour_mask;
			
			double _a,
				   _b,
				   _x,
				   _y,
				   _theta,
				   _r;
			public:
				/**@fn
				 * @return r
				 */
				inline double r() const
				{
					return _r;
				}
				
				/**@fn
				 * @return a
				 */
				inline double a() const
				{
					return _a;
				}
				
				/**@fn
				 * @return b
				 */
				inline double b() const
				{
					return _b;
				}
				
				/**@fn
				 * @return x
				 */
				inline double x() const
				{
					return _x;
				}
				
				/**@fn
				 * @return y
				 */
				inline double y() const
				{
					return _y;
				}

				/**@fn
				 * @return theta
				 */
				inline double theta() const
				{
					return _theta;
				}
			
				   
			protected:
				double _score;
			public:
				inline const double  & score() const
				{
					return _score;
				}
			protected:
				
			unsigned int _nb_points_ransac;
			unsigned int _nb_iter_ransac;
		
			void fit_and_select_best_region();
				convex_boundary<double> * _convex_boundary;
				
				double * _convex_contour;
				double * _convex_contour_r;
				
				
				gsl_rng * rng;
				fitting_ellipse ellipse;
				unsigned int _nb_iterations_ellipse;
				double fit_ellipse( double & x, 
									double & y, 
									double & a, 
									double & b, 
									double & theta );
		
				double get_distance( 	const double & x, 
										const double & y, 
										const double & a, 
										const double & b, 
										const double & theta ) const;
		
			double _ransac_threshold;
		
		
			/**@fn
			 * @brief
			 * Lib. mémoire
			 **/
			void free();

			/**@fn
			 * @brief
			 * Ini. Objet
			 **/
			void initialize();

			/**@fn
			 * @brief
			 * Allocation mémoire
			 * 
			 */
			void alloc();


			//Message d'erreurs
			ostream * err_stream;
			
			IplConvKernel 	* structuring_element_1,
							* structuring_element_2;
							
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
			
		protected:
			double _pupil_threshold;
		public:
			inline double pupil_threshold() const
			{
				return _pupil_threshold;
			}
			
			
	};
	
#endif
