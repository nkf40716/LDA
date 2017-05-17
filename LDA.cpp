#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "base_types.h"
#include "LDAApi.h"

#ifdef _WIN32
	typedef __int64 INT64;
#else
	typedef long long INT64;
#endif

int GeneralizedEigenvalueDecomposition(int n, double *a, double *b, double *eigenvector, double *eigenvalRe, double *eigenvalIm);

static int mat_mul(const double x[], const int xx, const int xy, const double y[], const int yx, const int yy, double a[])
{
	int i, j, k;
	double *wx, *wy, *wa;
	double *t = NULL;

	if (xy != yx || xx == 0 || yy == 0) {
		fprintf(stderr, "Invalid matrix size x= %d*%d,y= %d*%d\n", xx, xy, yx, yy);
		return -1;
	}

	wx = (double *)x;
	wa = a;
	if (x == a || y == a) {
		t = new double [xx * yy];
		wa = t;
	}

	for (i = 0; i < xx; i++) {
		for (j = 0; j < yy; j++) {
			wy = (double *)&y[j];
			*wa = 0;
			for (k = 0; k < xy; k++) {
				*wa += (double)(*wx) * (double)(*wy);
				wx++;
				wy += yy;
			}
			wx -= xy;
			wa++;
		}
		wx += xy;
	}

	if (t) {
		memcpy(a, t, sizeof(double) * (xx * yy));
		delete [] t;
	}

	return 0;
}

static int mat_add(const double x[], const double y[], const int xx, const int xy, double a[])
{
	int i, n;
	n = xx * xy;
	for (i = 0; i < n; i++)
		a[i] = x[i] + y[i];
	return 0;
}

static int mat_scal(const double x[], const int xx, const int xy, const double scaler, double a[])
{
	int i, n;
	n = xx * xy;
	for (i = 0; i < n; i++)
		a[i] = x[i] * scaler;
	return 0;
}

//=============================================================================

#ifndef MAX
#define MAX(a,b)	((a) >= (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)	((a) <= (b) ? (a) : (b))
#endif

typedef struct LDA {
	INT count;
	INT d;
	INT q;
	INT *N;
	double **S;
	double **C;
	double *mean;
	BOOL bTrained;
} LDA;

HANDLE LDA_Create(INT d, INT q)
{
	LDA *pLDA = NULL;
	INT i;

	if (d <= 0 || q <= 0 || q > d)
		return NULL;

	pLDA = (LDA *)malloc(sizeof(*pLDA));
	if (pLDA == NULL)
		return NULL;

	memset(pLDA, 0, sizeof(*pLDA));

	pLDA->d = d;
	pLDA->q = q;
	pLDA->count = 0;
	pLDA->N = (INT *)calloc(q, sizeof(INT));
	pLDA->S = (double **)calloc(q, sizeof(double *));
	pLDA->C = (double **)calloc(q, sizeof(double *));
	if (pLDA->N == NULL || pLDA->S == NULL || pLDA->C == NULL){
		LDA_Release((HANDLE)pLDA);
		return NULL;
	}
	pLDA->S[0] = (double *)calloc(q * d * d, sizeof(double));
	pLDA->C[0] = (double *)calloc(q * d, sizeof(double));
	if (pLDA->S[0] == NULL || pLDA->C[0] == NULL){
		LDA_Release((HANDLE)pLDA);
		return NULL;
	}
	for (i = 1; i < q; i++){
		pLDA->S[i] = pLDA->S[i-1] + (d * d);
		pLDA->C[i] = pLDA->C[i-1] + d;
	}

    pLDA->mean = (double *)calloc(d, sizeof(double));
	if (pLDA->mean == NULL){
		LDA_Release((HANDLE)pLDA);
		return NULL;
	}

	pLDA->bTrained = FALSE;

    return pLDA;
}

INT LDA_Release(HANDLE hLDA)
{
	LDA *pLDA = (LDA *)hLDA;
	if (pLDA == NULL)
		return -1;
	if (pLDA->N)
		free(pLDA->N);
	if (pLDA->S[0])
		free(pLDA->S[0]);
	if (pLDA->S)
		free(pLDA->S);
	if (pLDA->C[0])
		free(pLDA->C[0]);
	if (pLDA->C)
		free(pLDA->C);
	if (pLDA->mean)
		free(pLDA->mean);
    free(pLDA);
	return 0;
}

INT LDA_Add(HANDLE hLDA, double *v, INT k)
{
	LDA *pLDA = (LDA *)hLDA;
    INT d, q;
	INT i, j;

	if (pLDA == NULL)
		return -1;
	if (pLDA->bTrained)
		return -1;

	d = pLDA->d;
	q = pLDA->q;
	if (k >= q)
		return -1;

	for (i = 0; i < d; i++){
        pLDA->mean[i] += v[i];
		pLDA->C[k][i] += v[i];
        for(j = i; j < d; j++)
            pLDA->S[k][j + i*d] += v[i]*v[j];
    }
	pLDA->N[k]++;
    pLDA->count++;
	return pLDA->count;
}

INT LDA_Solve(HANDLE hLDA, double *eigenvector, double *eigenvalue)
{
	LDA *pLDA = (LDA *)hLDA;
	INT d, q;
    INT i, j, k;
	double *Sw = NULL;
	double *Sb = NULL;
	double *t_n = NULL;
	double *t_n_n = NULL;

	if (pLDA == NULL)
		return -1;
	if (pLDA->bTrained)
		return -1;

	d = pLDA->d;
	q = pLDA->q;
	pLDA->bTrained = TRUE;

	for (j = 0; j < d; j++)
		pLDA->mean[j] /= pLDA->count;

	Sw = (double *)calloc(d * d, sizeof(double));
	Sb = (double *)calloc(d * d, sizeof(double));
	t_n = (double *)calloc(d, sizeof(double));
	t_n_n = (double *)calloc(d * d, sizeof(double));
	if (Sw == NULL || Sb == NULL || t_n == NULL || t_n_n == NULL)
		goto L_ERROR;

	for (k = 0; k < q; k++){
		for (j = 0; j < d; j++){
			pLDA->C[k][j] /= pLDA->N[k];
			t_n[j] = pLDA->C[k][j] - pLDA->mean[j];
			for (i = 0; i <= j; i++){
				pLDA->S[k][j + i*d] -= (pLDA->C[k][i] * pLDA->C[k][j]) * pLDA->N[k];
				pLDA->S[k][i + j*d] = pLDA->S[k][j + i*d];
			}
		}
		mat_add(Sw, pLDA->S[k], d, d, Sw);
		mat_mul(t_n, d, 1, t_n, 1, d, t_n_n);
		mat_scal(t_n_n, d, d, pLDA->N[k], t_n_n);
		mat_add(Sb, t_n_n, d, d, Sb);
	}

	// eigenvector & eigenvalue
	if (GeneralizedEigenvalueDecomposition(d, Sb, Sw, eigenvector, eigenvalue, NULL) != 0)
		goto L_ERROR;

	free(Sw);
	free(Sb);
	free(t_n);
	free(t_n_n);

    return 0;

L_ERROR:

	if (Sw)
		free(Sw);
	if (Sb)
		free(Sb);
	if (t_n)
		free(t_n);
	if (t_n_n)
		free(t_n_n);

	return -1;
}
