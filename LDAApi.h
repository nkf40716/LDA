#ifndef __LDAApi_h__
#define __LDAApi_h__

#include "base_types.h"

#ifdef __cplusplus
extern "C"{
#endif

HANDLE LDA_Create(INT d, INT q);
INT LDA_Release(HANDLE hLDA);
INT LDA_Add(HANDLE hLDA, double *v, INT k);
INT LDA_Solve(HANDLE hLDA, double *eigenvector, double *eigenvalue);

#ifdef __cplusplus
}
#endif


#endif	//__LDAApi_h__
