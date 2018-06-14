#ifndef GSL_TRIANGLE_MATRIX_HPP_INCLUDED
#define GSL_TRIANGLE_MATRIX_HPP_INCLUDED
	#include <gsl/gsl_matrix.h>
    /**\fn void gsl_triangle_matrix(gsl_matrix * matrix);
     * @param matrix : matrice dont on doit supprimer le tiangle inférieur.
     * @brief
     * Cette fonction supprime le triangle inférieur de la matrice, la rendant triangulaire supérieur.
     */
    void gsl_triangle_matrix(gsl_matrix * matrix);

#endif // GSL_TRIANGLE_MATRIX_HPP_INCLUDED
