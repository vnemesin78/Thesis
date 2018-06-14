/**@file c_preprocessing.hpp
 * @author 
 * Valérian Némesin
 * 
 * @brief
 * Ce fichier contient le prototype de la classe de gestion du prétraitement de l'image. Le code suivant permet d'initialiser cette classe à partir du fichier argv[1] et de traiter l'image argv[2]
 * @code
#include <opencv/highgui.h>
#include ".../API/lib_api.hpp" //Mettre le bon chemin
#include ".../FresnIsis/include/c_preprocessing.hpp" //idem
int main( int argc, char ** argv)
{
	api_parameters params;
	if ( argc < 3 ) // 
	{
		cout << "Erreur: Pas assez d'arguments!" << endl;
		return 1;
	}
	//Chargement du fichier de config.
	if ( params.load(argv[1]) )
	{
		cout << "Erreur: Lecture de " << argv[1] <<"!" << endl;
		return 1;
	}
	//Lecture de l'image
	IplImage * image = cvLoadImage( argv[2] );
	if ( image == NULL )
	{
		cout << "Erreur: Lecture de " << argv[2] <<"!" << endl;
		return 1;
	}
	//Création de l'objet
	c_preprocessing preprocessing;
	if ( preprocessing.setup( params, image->width, image->height )
	{
		return 1;
	}

	if ( preprocessing.process( image ) < 0 )
	{
		cout << "Erreur : Problème impossible" << endl;
		return 1;
	}
	cvShowImage ("image originale", image );
	cvShowImage ("image traitée", preprocessing.image() );
	//Pause
	cvWaitKey(0);
	
	cvReleaseImage( &image);
	return 0;
}
* @brief Le fichier de configuration est le suivant
* @code
preprocessing::NAME = ""

%gamma (0 - 1.0)
pre_processing::gamma = 0

%epsilon ( 0 - 1.0 )
pre_processing::epsilon = 0
 **/
#ifndef _C_PREPROCESSING_HPP_
	#define _C_PREPROCESSING_HPP_
	#include "lib_api.hpp"
	#include "lib_image.hpp"
	#include <opencv/cv.h>
	#include <sstream>
	#include <iostream>
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>

	
	using namespace std;
	/**@class c_preprocessing
	 * @brief
	 * Classe contenant les méthodes de prétraitement des images.
	 **/
	class c_preprocessing
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 **/
			c_preprocessing();
			
			/**@fn 
			 * @param width : largeur maximale de travail
			 * @param height: hauteur maximale de travail
			 * @param opening : Ouverture dans le masque de la pupille
			 * @param erosion : Erosion dans le masque de l'iris
			 * @param median_filter : taille du filtre median pour l'image de la pupille
			 * @param err_stream : flux d'erreur ( NULL pour désactiver )
			 * @brief
			 * Constructeur
			 **/
			c_preprocessing( 	unsigned int width,
								unsigned int height,
								unsigned int closing_1,
								unsigned int opening_1,
								unsigned int closing_2,
								unsigned int opening_2,
								unsigned int median_filter,
								ostream * _err_stream = NULL );
			
			/**@fn 
			 * @param width : largeur maximale de travail
			 * @param height: hauteur maximale de travail
			 * @param opening : Ouverture dans le masque de la pupille
			 * @param erosion : Erosion dans le masque de l'iris
			 * @param median_filter : taille du filtre median pour l'image de la pupille
			 * @param err_stream : flux d'erreur ( NULL pour désactiver )
			 * @brief
			 * Constructeur
			 **/
			int setup ( 	unsigned int width,
							unsigned int height,
							unsigned int closing_1,
							unsigned int opening_1,
							unsigned int closing_2,
							unsigned int opening_2,
							unsigned int median_filter,
							ostream * _err_stream = NULL  );
							
			/**@fn 
			 * @param width : largeur maximale de travail
			 * @param height: hauteur maximale de travail
			 * @param err_stream : flux d'erreur ( NULL pour désactiver )
			 * @brief
			 * Constructeur
			 **/
			int default_setup( 	unsigned int width,
								unsigned int height,
								ostream * _err_stream = NULL  );
			

			/**@fn 
			 * @param params : paramètres
			 * @param n_space : espace de nom des différentes variables
			 * @param gamma_name : nom de la variable contenant le gamme
			 * @param threshold_name : nom de la variable pour le seuillage des spots
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @brief
			 * Setup.
			 * @return
			 * - -1 si problème dans les variables
			 * - 1 si problème dans les valeurs
			 * - 0 si réussite de la construction de l'objet
			 **/
			int setup ( api_parameters & params,
						unsigned int width,
						unsigned int height,
						ostream * _err_stream = NULL,
						const char * n_space = "pre_processing",
						const char * closing_1_name = "closing_1",
						const char * opening_1_name = "opening_1",
						const char * closing_2_name = "closing_2",
						const char * opening_2_name = "opening_2",
						const char * median_name = "median"
						);
			
			
			/**@fn
			 * @param image : image
			 * @brief
			 * Cette méthode effectue :
			 * - 1/ Un réhaussement d'histogramme
			 * - 2/ Une correction du gamma.
			 * 
			 **/
			int process ( const IplImage * image );
								
			/**@fn 
			 * @brief
			 * Destructeur.
			 **/
			~c_preprocessing();
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
			//Paramètres
			
			unsigned int _width_max;
		public:
			inline unsigned int width_max() const
			{
				return _width_max;
			}
			
			
		protected:
			unsigned int _height_max;
		public:
			inline unsigned int height_max() const
			{
				return _height_max;
			}
		
		protected:
			unsigned int _opening_1;
		public:
			inline unsigned int opening_1() const
			{
				return _opening_1;
			}

		protected:
			unsigned int _closing_1;
		public:
			inline unsigned int closing_1() const
			{
				return _closing_1;
			}
		protected:
			unsigned int _closing_2;
		public:
			inline unsigned int closing_2() const
			{
				return _closing_2;
			}

		protected:
			unsigned int _opening_2;
		public:
			inline unsigned int opening_2() const
			{
				return _opening_2;
			}


		protected:
			unsigned int _median_filter;
		public:
			inline unsigned int median_filter() const
			{
				return _median_filter;
			}

		protected:
			//Paramètres
			unsigned int _width;
		public:
			inline unsigned int width() const
			{
				return _width;
			}
			
			
		protected:
			unsigned int _height;
		public:
			inline unsigned int height() const
			{
				return _height;
			}
		
		protected:
			IplImage * _smoothed_image;
			
		public:
			const IplImage * smoothed_image() const
			{
				return _smoothed_image;
			}
			
		protected:
			//Tmp
			//Op. Morphologiques
			IplConvKernel 	* structuring_element_1,
							* structuring_element_2,
							* structuring_element_3,
							* structuring_element_4;
			
			//Methodes classiques		
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			
			/**@fn 
			 * @brief
			 * Ini. de l'objet.
			 **/
			void initialize();
		
			/**@fn
			 * @brief
			 * Alloc. mémoier
			 * 
			 **/
			void alloc();

			//Paramètres
			//Gestion des erreurs
		protected:
			ostream * err_stream;
	};
	
#endif
