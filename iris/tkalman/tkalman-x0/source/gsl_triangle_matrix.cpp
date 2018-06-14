#include "gsl_triangle_matrix.hpp"
/**\fn void gsl_triangle_matrix(gsl_matrix * matrix);
 * @param matrix : matrice dont on doit supprimer le tiangle inf�rieur.
 * @brief
 * Cette fonction supprime le triangle inf�rieur de la matrice, la rendant triangulaire sup�rieur.
 */
void gsl_triangle_matrix(gsl_matrix * matrix)
{
	unsigned int i,
		j;
	unsigned int step;
	for (i = 0; i < matrix->size1; i ++)
	{
		if (i > matrix->size2)
			step = matrix->size2;
		else
			step = i;
		for (j = 0; j < step; j ++)
		{
			matrix->data[i * matrix->tda + j] = 0;

		}

	}
}
