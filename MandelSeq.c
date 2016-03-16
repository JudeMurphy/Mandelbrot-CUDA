#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <Time.h>

int main(int argc, char **argv)
{
	int height, width;
	double xcenter, ycenter;
	double resolution;
	double gamma;
	int max_iter;
	FILE *outfp;
	double **matrixR;
	double **matrixG;
	double **matrixB;
	int i, j, r, c;

	height = atoi(argv[1]);
	width = atoi(argv[2]);
	xcenter = atof(argv[3]);
	ycenter = atof(argv[4]);
	resolution = atof(argv[5]);
	gamma = atof(argv[6]);
	max_iter = atoi(argv[7]);
	outfp = fopen(argv[8], "w");

	clock_t start = clock(), diff;

	matrixR = (double **) malloc(height*sizeof(int *));
	for (i=0; i<height; i++) 
	{
		matrixR[i] = (double*) malloc(width*sizeof(double));
	}

	matrixG = (double **) malloc(height*sizeof(double *));
	for (i=0; i<height; i++) 
	{
		matrixG[i] = (double*) malloc(width*sizeof(double));
	}

	matrixB = (double **) malloc(height*sizeof(int *));
	for (i=0; i<height; i++) 
	{
		matrixB[i] = (double*) malloc(width*sizeof(double));
	}

	double xoffset = -(width-1)/2.0;
	double yoffset = (height-1)/2.0;

	for (r=0; r<height; r++)
	{
		double y = ycenter + (yoffset - r)/resolution;
		for (c=0; c<width; c++)
		{
			double x = xcenter + (xoffset + c)/resolution;
			int iter=0;
			double a=0.0, b=0.0, a_old=0.0, b_old=0.0;
			double dist_sqr = 0.0;
			while (iter<max_iter && dist_sqr<=4.0)
			{
				iter++;
				a = a_old*a_old - b_old*b_old + x;
				b = 2.0*a_old*b_old + y;
				dist_sqr = a*a + b*b;
				a_old = a;
				b_old = b;
			}
			if (iter == max_iter)
			{
				matrixR[r][c] = 0.0f;
				matrixG[r][c] = 0;
				matrixB[r][c] = 0;
			}
			else
			{
				matrixR[r][c] = (double)pow(((double) iter)/((double)max_iter), gamma);
				matrixG[r][c] = 1.0;
				matrixB[r][c] = 1.0;
			}
		}
	}

	diff = clock() - start;
	int msec = (int)(diff * 1000 / CLOCKS_PER_SEC);
	printf("EXECUTION TIME: %d seconds and %d milliseconds\n", msec/1000, msec%1000);


	for (i=0; i<height; i++)
	{
		for (j=0; j<width; j++)
			fprintf(outfp, " %lf %lf %lf", matrixR[i][j], matrixG[i][j], matrixB[i][j]);
		fprintf(outfp,"\n");
	}



	fclose(outfp);
}