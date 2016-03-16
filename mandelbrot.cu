//JUDE MURPHY
//PARALLEL AND SCIENTIFIC COMPUTING
//ASSIGNMENT 2

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 25000000

__global__ void plotMandelbrotSet(int width, int height, double xcenter, double ycenter, double resolution, double gamma, int max_iter, double *matrixR, double *matrixG, double *matrixB);

int main(int argc, char** argv)
{	
	//ALLOCATE ALL VARIABLES
	int height, width;
	double xcenter, ycenter;
	double resolution;
	double gamma;
	int max_iter;
	FILE *outfp;

	//GET ARGUMENTS FROM RUNTIME ARGS
	height = atoi(argv[1]);
	width = atoi(argv[2]);
	xcenter = atof(argv[3]);
	ycenter = atof(argv[4]);
	resolution = atof(argv[5]);
	gamma = atof(argv[6]);
	max_iter = atoi(argv[7]);
	outfp = fopen(argv[8], "w");

	//MAKE TOTAL GRID SPACE BE FROM THE HEIGHT TIMES THE WIDTH
	int totalGridSpace = height * width;

	//USED TO DETERMINE THE TIME TAKEN TO COMPLETE THE PROGRAM
	cudaEvent_t cudaStart, cudaEnd;
	float elapsedTime;

	//NEED FOR TIME TO START
	srand(time(NULL));

	//_________________________________________________________________________
	//ALLOCATES MEMORY FOR THREE ARRAYS ON THE CPU	
	double * matrixR = (double *) malloc(totalGridSpace * sizeof(double));
	double * matrixG = (double *) malloc(totalGridSpace * sizeof(double));
	double * matrixB = (double *) malloc(totalGridSpace * sizeof(double));

	//_________________________________________________________________________
	//CREATE CUDA EVENTS AND START RECORDING TIME
	cudaEventCreate(&cudaStart);
	cudaEventCreate(&cudaEnd);
	cudaEventRecord(cudaStart, 0);

	//_________________________________________________________________________
	//ALLOCATES MEMORY FOR THREE ARRAYS ON THE GPU
	double * dev_R, *dev_G, *dev_B;
	cudaMalloc((void**)&dev_R, totalGridSpace * sizeof(double));
	cudaMalloc((void**)&dev_G, totalGridSpace * sizeof(double));
	cudaMalloc((void**)&dev_B, totalGridSpace * sizeof(double));

	//COMPUTE MANDELBROT SET IN PARALLEL
	plotMandelbrotSet<<<256, 256>>>(width, height, xcenter, ycenter, resolution, gamma, max_iter, dev_R, dev_G, dev_B);

	//COPY DATA BACK TO CPU AFTER FINISHING COMPUTATION
	cudaMemcpy(matrixR, dev_R, (totalGridSpace * sizeof(double)), cudaMemcpyDeviceToHost);
	cudaMemcpy(matrixG, dev_G, (totalGridSpace * sizeof(double)), cudaMemcpyDeviceToHost);
	cudaMemcpy(matrixB, dev_B, (totalGridSpace * sizeof(double)), cudaMemcpyDeviceToHost);

	//GET THE END OF THE TIME IT TOOK TO CALCULATE
	cudaEventRecord(cudaEnd, 0);
	cudaEventSynchronize(cudaEnd);	
	cudaEventElapsedTime(&elapsedTime, cudaStart, cudaEnd);

	//CLOSE OUT EVENTS
	cudaEventDestroy(cudaStart);
	cudaEventDestroy(cudaEnd);

	printf("CUDA Elapsed Time: %3.3f sec\n", elapsedTime/1000);

	//WRITING OUT TO THE FILE
	for (int i = 0; i < totalGridSpace; i++)
	{
		fprintf(outfp, " %lf %lf %lf\n", matrixR[i], matrixG[i], matrixB[i]);
	}

	//FREE MEMORY
	cudaFree(dev_R);
	cudaFree(dev_G);
	cudaFree(dev_B);
	free(matrixR);
	free(matrixG);
	free(matrixB);	

	//RESET CUDA DEVICE
	cudaDeviceReset();

	return 0;
}

__global__ void plotMandelbrotSet(int width, int height, double xcenter, double ycenter, double resolution, double gamma, int max_iter, double *matrixR, double *matrixG, double *matrixB)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	int totalNumberOfBoxes = height * width;

	while (id < totalNumberOfBoxes)
	{
		int c = id / width;
		int r = id % width;
		
		int currentIndex = c + (r * width);

		double xoffset = -(width-1)/2.0;
		double yoffset = (height-1)/2.0;
		
		double x = xcenter + (xoffset + c)/resolution;
		double y = ycenter + (yoffset - r)/resolution;

		int iter = 0;
		double a = 0.0, b = 0.0, a_old = 0.0, b_old = 0.0;
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
			matrixR[currentIndex] = 0.0f;
			matrixG[currentIndex] = 0;
			matrixB[currentIndex] = 0;
		}
		else
		{
			matrixR[currentIndex] = (double)pow(((double) iter)/((double)max_iter), gamma);
			matrixG[currentIndex] = 1.0;
			matrixB[currentIndex] = 1.0;
		}

		id += blockDim.x * gridDim.x;
	}	
}
