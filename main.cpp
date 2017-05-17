#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LDAApi.h"

int DimReduction(const double *eigenvector, const double *v, const int d, double *u, const int k)
{
	double *W = (double *)eigenvector;
	double *y = NULL;
	int i, j;

	if (W == NULL)
		return -1;
	if (v == NULL || u == NULL)
		return -1;
	if (d <= 0 || k <= 0 || d < k)
		return -1;

	y = u;
	if (v == u){
		y = new double [k];
	}

	// Y=X*W;
	for (i = 0; i < k; i++){
		y[i] = 0;
		for (j = 0; j < d; j++)
			y[i] += v[j] * W[i+j*d];
	}

	if (v == u){
		memcpy(u, y, sizeof(double) * k);
		delete [] y;
	}

	return 0;
}

int main(int argc, char *argv[])
{
    HANDLE hLDA;
    int i, j, k;
#define LEN 4
    double eigenvector[LEN*LEN];
    double eigenvalue[LEN];

    hLDA= LDA_Create(LEN, 3);

	FILE *fp = fopen(argv[1], "r");
	if (fp == NULL){
		printf("ERROR: failed to open file [%s].\n", argv[1]);
		return -1;
	}

	float f;
	double v[LEN+1];
	while (!feof(fp)){
		for (j = 0; j < LEN+1; j++){
			if (fscanf(fp, "%f", &f) != 1)
				break;
			v[j] = f;
		}
		if (j == LEN+1)
			LDA_Add(hLDA, v, (int)v[LEN] - 1);
	}

    LDA_Solve(hLDA, eigenvector, eigenvalue);
	LDA_Release(hLDA);

    for(i=0; i<LEN; i++){
        for(j=0; j<LEN; j++){
            printf("%9.4f ", eigenvector[i + j*LEN]);
        }
        printf("  %9.4f %f\n", eigenvalue[i], eigenvalue[i]/eigenvalue[0]);
    }
	printf("\n");

	fseek(fp, 0L, SEEK_SET);
	i = 0;
	while (!feof(fp)){
		for (j = 0; j < LEN+1; j++){
			if (fscanf(fp, "%f", &f) != 1)
				break;
			v[j] = f;
		}
		if (j == LEN+1){
			DimReduction(eigenvector, v, LEN, v, LEN);
			printf("%d.", i+1);
			for (j = 0; j < LEN; j++)
				printf(" %.2f,", v[j]);
			printf("\n");
			i++;
		}
	}

    return 0;
}
