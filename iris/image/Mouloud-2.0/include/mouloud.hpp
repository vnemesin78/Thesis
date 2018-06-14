/**@file Mouloud.hpp
 * 
 */
#ifndef _MOULOUD_HPP_
#define _MOULOUD_HPP_
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "typedef_pdf.hpp"
#include <iostream>
using namespace std;

#define MOULOUD_COMPUTE_EDGES(type)\
template \
int Mouloud :: compute_edges( 	double * radii,\
								void * bg_params,\
								void * obj_params,\
								const type * data,\
								const unsigned char * data_mask,\
								unsigned int width,\
								unsigned int height,\
								unsigned int width_step,\
								unsigned int mask_width_step );



/**@class 
 * @brief
 * Classe pour l'application de la méthode de Mouloud
 * 
 */
class Mouloud
{
public:
	/**@fn
	 * @brief
	 * Constructeur
	 * 
	 */
	Mouloud ( void );
	
	/**@fn
	 * @param nb_iter_gem : nombre d'itération de l'algorithme GEM
	 * @param nb_iter_icm : nombre d'itération de l'algorithme ICM
	 * @param bg_function : function du background
	 * @brief
	 * Constructeur
	 * 
	 */
	Mouloud (	unsigned int nb_iter_gem,
				unsigned int nb_iter_icm,
				pdf_type bg_function,
				estimation_type est_bg_function,		
				const void * est_bg_params,
				unsigned int size_est_bg_params,
				pdf_type obj_function,
				estimation_type est_obj_function,
				const void * est_obj_params,
				unsigned int size_est_obj_params,			
				pdf_field_type edge_function,
				const void * edge_params,
				unsigned int size_edge_params,
				unsigned int width = 0,
				unsigned int height = 0 );
	/**@fn
	 * @brief
	 * Setup
	 * 
	 */				
	int setup ( );

	/**@fn
	 * @brief
	 * Setup
	 * 
	 */		
	int setup (	unsigned int nb_iter_gem,
				unsigned int nb_iter_icm,
				pdf_type bg_function,
				estimation_type est_bg_function,		
				const void * est_bg_params,
				unsigned int size_est_bg_params,
				pdf_type obj_function,
				estimation_type est_obj_function,
				const void * est_obj_params,
				unsigned int size_est_obj_params,			
				pdf_field_type edge_function,
				const void * edge_params,
				unsigned int size_edge_params,
				unsigned int width = 0,
				unsigned int height = 0 );
	
	/**@fn
	 * @brief
	 * Destructeur
	 * 
	 */
	~Mouloud();
	
	
	template<class type> 
	int compute_edges( 	double * radii,
						void * bg_params,
						void * obj_params,
						const type * data,
						const unsigned char * data_mask,
						unsigned int width,
						unsigned int height,
						unsigned int width_step,
						unsigned int mask_width_step );
	
	
	
	
	/**@fn
	 * @param radii : rayons
	 * @param bg_
	 * 
	 */
	int compute_edges( 	double * radii,
						void * bg_params,
						void * obj_params,
						const double * data,
						const unsigned char * data_mask,
						unsigned int width,
						unsigned int height,
						unsigned int width_step,
						unsigned int mask_width_step );
	
	/**@fn
	 * @param radii : rayons
	 * @param bg_
	 * 
	 */
	int compute_edges( 	double * radii,
						void * bg_params,
						void * obj_params,
						const IplImage * image,
						const IplImage * mask );
	
	
	
	//Accesseur
	inline const IplImage * mask() const
	{
		return _mask;
	}
	
	inline const IplImage * log_bg_pdf() const
	{
		return _log_bg_pdf;
	}
	
	inline const IplImage * log_obj_pdf() const
	{
		return _log_obj_pdf;
	}
	
	inline const IplImage * log_radii_pdf() const
	{
		return _log_r_pdf;
	}
	
	inline unsigned int nb_iter_gem() const
	{
		return _nb_iter_gem;
	}
	
	inline unsigned int nb_iter_icm() const
	{
		return _nb_iter_icm;
	}
	
	
	
	
	
	
	
	
	
	
	
protected:
	ostream * err_stream;

	unsigned int  _width,
				    _height;

	IplImage * _log_bg_pdf,
			 * _log_obj_pdf,
			 * _log_r_pdf,
			 * _mask,
			 * _tmp_image;

	//Params
	unsigned int 	_nb_iter_icm,
					_nb_iter_gem;

	pdf_type _bg_function;
	estimation_type _est_bg_function;		
		void * _est_bg_params;
		unsigned int _size_est_bg_params;
	pdf_type _obj_function;
	estimation_type _est_obj_function;
		void * _est_obj_params;
		unsigned int _size_est_obj_params;
		
	//Edge
	pdf_field_type _edge_function;
		void * _edge_params;
		unsigned int _size_edge_params;

	/**@fn
	 * @brief
	 * Calcul de la pdf des rayons
	 * 
	 */
	int compute_radii_pdf( 	double * radii,
							void * bg_params,
							void * obj_params,
							const double * data,
							const unsigned char * data_mask,
							unsigned int width,
							unsigned int height,
							unsigned int width_step,
							unsigned int mask_width_step );
	
	/**@fn
	 * @brief
	 * Calcul des pdf du background
	 * 
	 * 
	 */
	int compute_bg_pdf( 	double * radii,
							void * bg_params,
							void * obj_params,
							const double * data,
							const unsigned char * data_mask,
							unsigned int width,
							unsigned int height,
							unsigned int width_step,
							unsigned int mask_width_step );
	/**@fn
	 * @brief
	 * Calcul de la pdf des obj.
	 * 
	 */
	int compute_obj_pdf( 	double * radii,
							void * bg_params,
							void * obj_params,
							const double * data,
							const unsigned char * data_mask,
							unsigned int width,
							unsigned int height,
							unsigned int width_step,
							unsigned int mask_width_step );
	/**@fn
	 * @brief
	 * Estime les pdf.
	 * 
	 */
	int estimate_pdf( 	double * radii,
						void * bg_params,
						void * obj_params,
						const double * data,
						const unsigned char * data_mask,
						unsigned int width,
						unsigned int height,
						unsigned int width_step,
						unsigned int mask_width_step );
	/**@fn
	 * @brief
	 * Estime les rayons.
	 * 
	 */
	int estimate_radii ( 	double * radii,
							void * bg_params,
							void * obj_params,
							const double * data,
							const unsigned char * data_mask,
							unsigned int width,
							unsigned int height,
							unsigned int width_step,
							unsigned int mask_width_step );
	
	
	
	
	
	
	/**@fn
	 * @brief
	 * Free memory
	 */
	void _free();
	
	/**@fn
	 * @brief
	 * Set to 0
	 */
	void initialize();
	
	
	/**@fn
	 * @brief
	 * Memory allocation.
	 */
	void alloc();













	
};

	
	
	
	
	
	
	
	
	
	
	
	
#endif
