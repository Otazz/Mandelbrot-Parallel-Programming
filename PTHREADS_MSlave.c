#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <float.h>
#define N_THREADS 4
#define N_ESCRAVOS 4

typedef struct data{
    int tid;
} data;

typedef struct {
    int x;
    int y;
} int2;

typedef struct{
    double x;
    double y;
} double2;



float **mat;
double x = -DBL_MAX, dt=0.001;
int nColumns=2048, nLines=2048, i = 0, size;

pthread_mutex_t mutex;
pthread_cond_t condM;
pthread_cond_t condE;

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

int printMatrixToFilePGM(float **mat,int tamx, int nLines, char *srcFile){

    FILE *arq=fopen(srcFile,"w");
    int cont, cont2;
    float min,max; 
    min=mat[0][0];
    max=mat[0][0];

    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<tamx;cont2++){
            if(min>mat[cont][cont2]) min=mat[cont][cont2];
            if(max<mat[cont][cont2]) max=mat[cont][cont2];

        }

    }

    max=max*0.35;
    float delta=max-min;
    fprintf(arq,"P2 \n");
    fprintf(arq,"#comentario qualquer \n");
    fprintf(arq,"%d\n%d \n",tamx,nLines);
    fprintf(arq,"255\n");

    for(cont=0;cont<nLines;cont++){
        for(cont2=0;cont2<tamx;cont2++){ 
            int valpixel=((mat[cont][cont2]-min)/delta)*255.0f;
            if(valpixel>255) valpixel=255;
            fprintf(arq,"%d \n", valpixel);
        } 
    } 

    fclose(arq);

}

float** mallocFloatMatrix(int tamx, int nLines, float defaultValueOfTheElementsAtMatrix){

    float **errorCodeReturn=0x0;
    float **mat;
    int i,j;
    int condErrorMalloc=0;//error indicator at malloc
    mat=malloc(sizeof(float *)*nLines);

    if(mat==0x0) return errorCodeReturn;//error at malloc return null pointer

    for(i=0;i<tamx;i++)
        mat[i]=malloc(sizeof(float )*tamx);

    for(i=0;i<tamx;i++){//detect error at malloc
        if(mat[i]==0x0){
            condErrorMalloc=1;
            break;
        }
    }



    if(condErrorMalloc==0){  
        return mat;
    }

    for(i=0;i<nLines;i++){
        for(j=0;j<tamx;j++)
            mat[i][j]=defaultValueOfTheElementsAtMatrix;
    }

    for(i=0;i<tamx;i++)
        if(mat[i]!=0x0) free(mat[i]);
        
    free(mat);
    
    return errorCodeReturn;

}

void freeFloatMatrix(float **mat,int tamx, int nLines){

    int i;
    
    for(i=0;i<nLines;i++){
        if(mat[i]!=0x0) free(mat[i]);
    }
    free(mat);

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



void *mestre(void *param){

	for(i=0;i<size;i++){//real component of C at $z_{n+1}=z_n+C$
		pthread_mutex_lock(&mutex);
        x=-2.0+((double)i*dt);
        pthread_cond_signal(&condE);
    	pthread_cond_wait(&condM, &mutex);
    	pthread_mutex_unlock(&mutex);

    }
}

void *escravo(void *param){

	data *dados = (data *)param;
	int k, ite=600,progress=0;
	
	while (1) {
      if (i >= size) {
		printf("%d saiu. \n", dados->tid);
        pthread_mutex_unlock(&mutex);
        pthread_exit(0);
      }

      while (x == -DBL_MAX) {
         pthread_cond_signal(&condM);
         pthread_cond_wait(&condE, &mutex);  
      }

      double x1 = x;
	  x = -DBL_MAX;
	  pthread_cond_signal(&condM);
      pthread_mutex_unlock(&mutex);

		double y;
		for(y=-2.0;y<2.0;y=y+dt){//imaginary component of C at $z_{n+1}=z_n+C$

        	int2* iterationPath=(int2 *)malloc(sizeof(int2)*ite);
        	if(iterationPath==0x0) return 0x0;
        	int completedIterations = iteration(x1,y,nColumns,nLines,ite, iterationPath);//completedIterations= quantity of elements at vector iterationPath

			for(k=0;k<completedIterations;k++){
            	if(iterationPath[k].x!=-1 && iterationPath[k].y!=-1)//test if a point z in the iteration k may be normalized to coordinates at matrix mat. 
                	mat[iterationPath[k].x][iterationPath[k].y]=mat[iterationPath[k].x][iterationPath[k].y]+1.0f;//increments a point in matrix, this point is pointed by z with  z points normalized.
        }
        	free(iterationPath);
        }

        progress++;
        if(progress%100 ==0)//print at screen information about progrees of the operation
            printf("[%d] %lf \n",dados->tid,x1);
        }
}

int main(void){
	
	pthread_t mestre_t;
    pthread_t threads[N_ESCRAVOS];
    pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&condM, NULL);
	pthread_cond_init(&condE, NULL);

	size = round(4.0/dt);
    mat=mallocFloatMatrix(nColumns,nLines,0.0f);
    if(mat==0x0) return 0;   

    int t;
    pthread_create(&mestre_t, NULL, mestre, NULL);
    data *dados = (data *)malloc(N_ESCRAVOS * sizeof(data));

    for (t=0; t<N_ESCRAVOS; t++){
    	dados[t].tid = t;
      	pthread_create(&threads[t], NULL, escravo, &dados[t]);
    }

    for (t=0; t<N_ESCRAVOS; t++)
    	pthread_join(threads[t], NULL);
    	
    printMatrixToFilePGM(mat,nColumns,nLines,"saida3.pgm");
    freeFloatMatrix(mat,nColumns,nLines);

    return 0;
    
    pthread_exit(NULL);
}
