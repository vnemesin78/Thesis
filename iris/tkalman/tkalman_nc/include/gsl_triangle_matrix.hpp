#ifndef GSL_TRIANGLE_MATRIX_HPP_INCLUDED
#define GSL_TRIANGLE_MATRIX_HPP_INCLUDED
	#include <gsl/gsl_matrix.h>
    /**\fn void gsl_triangle_matrix(gsl_matrix * matrix);
     * @param matrix : matrice dont on doit supprimer le tiangle inf�rieur.
     * @brief
     * Cette fonction supprime le triangle inf�rieur de la matrice, la rendant triangulaire sup�rieur.
     */
    void gsl_triangle_matrix(gsl_matrix * matrix);

#endif // GSL_TRIANGLE_MATRIX_HPP_INCLUDED
