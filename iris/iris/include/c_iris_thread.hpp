#ifndef _C_IRIS_THREAD_HPP_
#define _C_IRIS_THREAD_HPP_
#include <opencv/cv.h>
#include "c_pupil_thread.hpp"
#include "c_focus_score.hpp"
#include "c_iris_segmentation.hpp"
#include "c_eyelids_segmentation.hpp"
#include "c_eyelid_segmentation.hpp"
#include "c_pupil_thread.hpp"
#include "c_polar_iris.hpp"
#include "c_iris_code.hpp"
#include <sys/stat.h>
#include "iris_data.hpp"


/**@class
 * @brief
 * Cette classe gère les différentes opérations pour la segmentation de l'iris.
 */
class c_iris_thread
{
public:
	/**@fn
	 * @brief
	 * Constructeur
	 */
	c_iris_thread();

	/**@fn
	 * @brief
	 * Reset.
	 * 
	 */
	int setup ( void );


	int default_setup(  ostream * _err_stream = NULL );

	/**@fn
	 * @brief
	 * Setup.
	 */
	int setup( char ** argv,
				unsigned int argc,
				ostream * _err_stream = NULL );

	/**@fn
	 * @brief
	 * Routine
	 */
	int segment_iris( const pupil_data & p_data );


	/**@fn
	 * @brief
	 * Destructeur
	 */
	~c_iris_thread();
	
	/**@fn
	 * @brief
	 * Renvoie les données
	 * 
	 */
	int get_data ( iris_data & data );
	
	/**@fn
	 * @brief
	 * Reset des objets
	 * 
	 */
	void reset();
	
	inline const iris_data & iris_seg_data() const
	{
		return *i_data;
	}
	/**@fn
	 * @param error_str : flux d'erreur
	 * @brief
	 * Modifie le flux d'erreur.
	 * 
	 */
	inline void set_error_stream( ostream & error_str = cout )
	{
		err_stream = &error_str;
		iris_segmentation->set_error_stream( error_str );
		eyelid_segmentation->set_error_stream( error_str );
		polar_iris->set_error_stream( error_str );
		//~ focus_score->set_error_stream( error_str );
	}
protected:

	/**@fn
	 * @brief
	 * Lib. mémoire
	 * 
	 */
	void free();
	
	/**@fn
	 * @brief
	 * Ini. mémoire
	 * 
	 */
	void initialize();
	
	c_iris_segmentation * iris_segmentation;
	c_eyelids_segmentation * eyelid_segmentation;
	c_eyelid_segmentation * eyelid_obj;
	c_polar_iris * polar_iris;
		IplImage * tmp_p_img,
				 * tmp_p_mask;
			
	c_focus_score * focus_score;
	int _score_type;
	double _validness_factor;
	
	
	c_iris_code * iris_code;

	iris_data * i_data;

	//Flux d'erreur
	ostream * err_stream;
	
	double _spot_threshold;
	IplImage * mask_1,
			 * mask_2,
			 * mask_3;
	
	IplConvKernel 	* structuring_element_1,
					* structuring_element_2;
					
	double _n_sigma_up,
		   _n_sigma_low;
		   
	c_label * label_obj;
};







#endif
