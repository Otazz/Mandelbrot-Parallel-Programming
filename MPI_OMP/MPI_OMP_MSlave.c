#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <float.h>
typedef struct {
    int x;
    int y;
} int2;
typedef struct{
    double x;
    double y;
} double2;

int2 coordinatesConversion( double x,double y, int nColumns,int nLines){
    
    int2 ret;

    //error return code=================================================
    int2 retError;
    retError.x=-1;
    retError.y=-1;
    //end===============================================================
    
    

    ret.x=round(((2.0+x)/3.5) *((double)(nColumns-1)));
    ret.y=round(((1.5+y)/3.5) *((double)(nLines-1)));

    //invalid parameters for  x or y arguments==========================
    if(ret.x<0 || ret.x>=nColumns) return retError;
    if(ret.y<0 || ret.y>=nLines) return retError;
    //end===============================================================
    return ret;
}
int printMatrixToFilePGM(float *mat,int nCol, int nLines, char *srcFile){
    FILE *arq=fopen(srcFile,"w");

    int cont, cont2;
    float min,max; 
    min=mat[nCol*0+0];
    max=mat[nCol*0+0];
    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<nCol;cont2++){
            if(min>mat[cont*nCol + cont2]) min=mat[cont*nCol + cont2];
            if(max<mat[cont*nCol + cont2]) max=mat[cont*nCol + cont2];
        }
    }
    max=max*0.35;
    float delta=max-min;
    fprintf(arq,"P2 \n");
    fprintf(arq,"#comentario qualquer \n");
    fprintf(arq,"%d\n%d \n",nCol,nLines);
    fprintf(arq,"255\n");
    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<nCol;cont2++){ 
            int valpixel=((mat[cont*nCol + cont2]-min)/delta)*255.0f;
            if(valpixel>255) valpixel=255;
            fprintf(arq,"%d \n", valpixel);
        } 
    } 
    fclose(arq);
}
float* mallocFloatMatrix(int nCol, int nLines, int nProc, float defaultValueOfTheElementsAtMatrix){
    float *mat;
    int i;
    mat=malloc(sizeof(float)*nLines*nCol/nProc);
	#pragma omp parallel for
    for(i = 0; i < nLines * nCol/nProc; i++){
    	mat[i] = defaultValueOfTheElementsAtMatrix;
	}
	return mat;
}

int iteration(double x,double y, int nColumns,int nLines, int ite,int2 *iterationPath){

    int cont;    
    int condInvalidPointer=1;
    double2 z;
    z.x=0.0;
    z.y=0.0;
    double2 c;
    c.x=x;
    c.y=y;
    double2 zt;

    for(cont=0;cont<ite;cont++){
        //This  does z=z^2+c==================================================
        zt.x=((z.x*z.x)-(z.y*z.y))+c.x;
        zt.y=(2.0*(z.x*z.y))+c.y;
        z=zt;
        //end=================================================================
        if(((z.x*z.x)+(z.y*z.y))>4.0){
            if(cont>100)
                condInvalidPointer=0;
            break;
        }
        iterationPath[cont]=coordinatesConversion(z.x,z.y,nColumns,nLines);
    }

    if(condInvalidPointer)
        return 0;
    
    return cont;

}

int main(int argc, char *argv[]){
	omp_set_num_threads(2);
	int my_rank, p, provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	//MPI_Init(&argc, &argv);
  	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  	MPI_Comm_size (MPI_COMM_WORLD, &p);
  	MPI_Status status;
  	int nProc = p;
    //image size=========================================================
    int nColumns=2048;
    int nLines=2048;
    //end================================================================
    //size points========================================================
    double dt=0.001;//quantity of points going to increase with the  decrease  of the dt value
    int size=round(4.0/dt);//sizeOfPoints=size*size
    //end================================================================
    int ite=600;
    int len = nColumns*nLines;

    float *mat = mallocFloatMatrix(nColumns,nLines,1,0.0f);
    float *matS = mallocFloatMatrix(nColumns,nLines,1,0.0f);
    
    
    
    int i,j,k,target;   
    double x,y;
    int progress=0;
    if(my_rank == 0){
    	for(i=0;i<size;i++){//real component of C at $z_{n+1}=z_n+C$
	
      	if (i < size)
        	x=-2.0+((double)i*dt);
        
      	MPI_Recv (&target, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      	MPI_Send (&x, 1, MPI_DOUBLE, target, 0, MPI_COMM_WORLD);
	}
	for(i = 1; i < nProc; i++){
	x = -DBL_MAX;
		MPI_Recv (&target, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      	MPI_Send (&x, 1, MPI_DOUBLE, target, 0, MPI_COMM_WORLD);
}
}
	else{
		while (1) {
		    MPI_Send (&my_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      		MPI_Recv (&x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

      		if (x == -DBL_MAX) break;

			double y;
			for(y=-2.0;y<2.0;y=y+dt){//imaginary component of C at $z_{n+1}=z_n+C$

        	int2* iterationPath=(int2 *)malloc(sizeof(int2)*ite);
        	if(iterationPath==0x0) return 0x0;
        	int completedIterations = iteration(x,y,nColumns,nLines,ite, iterationPath);//completedIterations= quantity of elements at vector iterationPath

			for(k=0;k<completedIterations;k++){
            	if(iterationPath[k].x!=-1 && iterationPath[k].y!=-1)//test if a point z in the iteration k may be normalized to coordinates at matrix mat. 
                	matS[iterationPath[k].x*nColumns + iterationPath[k].y]=matS[iterationPath[k].x*nColumns+ iterationPath[k].y]+1.0f;//increments a point in matrix, this point is pointed by z with  z points normalized.
        	}

        	free(iterationPath);
        }

        progress++;
        if(progress%100 ==0)//print at screen information about progrees of the operation
           printf("[%d] %lf \n",my_rank,x);
        }
        
	
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Reduce(matS, mat, len, MPI_FLOAT,MPI_SUM, 0, MPI_COMM_WORLD);
	

	if(my_rank == 0)
		printMatrixToFilePGM(mat,nColumns,nLines,"saida3.pgm");

	
	
 MPI_Finalize(); 
    return 0;
}
