/********************************************************
 ***IMPORTANT NOTE***
 If you have problems with the provided sample code,
 part of the reason might be due to the function "fopen".
 Please try changing "r/w" to "rb/wb" or the other way
 when you use this function.
 *********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <sys/time.h>
#include "mpi.h"

#define max(x, y) ((x>y) ? (x):(y))
#define min(x, y) ((x<y) ? (x):(y))
#define winSize 20
#define seg1 5
#define seg2 80

int xdim;
int ydim;
int maxraw;
int rows;
unsigned char *image;

void ReadPGM(FILE*);
void WritePGM(FILE*);


int main(int argc, char **argv)
{
    int numproces, myid;
    FILE *fp;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);
    
    if(myid == 0){
        if (argc != 3){
            printf("Usage: MyProgram <input_ppm> <output_ppm> \n");
            printf("       <input_ppm>: PGM file \n");
            printf("       <output_ppm>: PGM file \n");
            exit(0);
        }
        
        /* begin reading PGM.... */
        printf("begin reading PGM.... \n");
        if ((fp=fopen(argv[1], "r"))==NULL){
            // printf("read error...\n");
            exit(0);
        }
        
        ReadPGM(fp);
        printf("Width=%d, Height=%d \nMaximum=%d\n", xdim, ydim, maxraw);
        
        if(ydim%(numproces-1) == 0)
            rows = ydim/(numproces-1);
        else
            rows = ydim/(numproces-1) + 1;
        
        int dest, index, i;
        unsigned char newImg[rows * xdim];
        struct timeval t1, t2, total;
        
        gettimeofday(&t1, NULL);
        
        for (dest = 1; dest < numproces; dest++) {
            MPI_Recv(&newImg, rows * xdim, MPI_UNSIGNED_CHAR, dest, 2, MPI_COMM_WORLD, NULL);
            for (index = (dest - 1) * rows * xdim, i = 0; i < rows * xdim; i++, index++)
                image[index] = newImg[i];
        }
        
        gettimeofday(&t2, NULL);
    	timersub(&t2, &t1, &total);
        
        printf("The total time cost is %f seconds\n", total.tv_sec + total.tv_usec/1000000.0);
        
        /* Begin writing PGM.... */
        printf("Begin writing PGM.... \n");
        if ((fp=fopen(argv[2], "wb")) == NULL){
            printf("write pgm error....\n");
            exit(0);
        }
        
        WritePGM(fp);
    }
    
    if(myid > 0){
        if ((fp=fopen(argv[1], "r"))==NULL){
            exit(0);
        }
        
        ReadPGM(fp);
        
        if(ydim%(numproces-1) == 0)
            rows = ydim/(numproces-1);
        else
            rows = ydim/(numproces-1) + 1;
        
        double sum, weightSum, weight;
        unsigned char tempImg[rows * xdim];
        double dictionary[winSize * winSize];
        int k = winSize / 2;
        
        for (int i = 0; i < winSize; i++) {
            for (int j = 0; j < winSize; j++) {
                dictionary[i * winSize + j] = (double)((i-k)*(i-k)+(j-k)*(j-k)) / (seg1*seg1);
            }
        }
        
        int i, col, row, indexX, indexY, m, n;
        for (row=(myid - 1) * rows, i = 0; i < rows; i++, row++) {
            for (col = 0; col < xdim; col++) {
                sum = 0.0;
                weightSum = 0.0;
                for (m = 0, indexY = row - k; m < winSize; m++, indexY++) {
                    if (indexY >= 0 && indexY < ydim) {
                        for (n = 0, indexX = col - k; n < winSize; n++, indexX++) {
                            if(indexX >= 0 && indexX < xdim) {
                                weight = exp(-dictionary[m * winSize + n] - pow((image[indexY * xdim +indexX] - image[row * xdim + col]) / seg2, 2.0));
                                sum += weight * image[indexY * xdim + indexX];
                                weightSum += weight;
                            }
                        }
                    }
                }
                tempImg[i * xdim + col] = (int)(sum / weightSum);
            }
        }
        
        MPI_Send(&tempImg, rows * xdim, MPI_UNSIGNED_CHAR, 0, 2, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    
    return (1);
}

void ReadPGM(FILE* fp)
{
    int c;
    int i,j;
    int val;
    unsigned char *line;
    char buf[1024];
    
    
    while ((c=fgetc(fp)) == '#')
        fgets(buf, 1024, fp);
    ungetc(c, fp);
    if (fscanf(fp, "P%d\n", &c) != 1) {
        // printf ("read error ....");
        exit(0);
    }
    if (c != 5 && c != 2) {
        // printf ("read error ....");
        exit(0);
    }
    
    if (c==5) {
        while ((c=fgetc(fp)) == '#')
            fgets(buf, 1024, fp);
        ungetc(c, fp);
        if (fscanf(fp, "%d%d%d",&xdim, &ydim, &maxraw) != 3) {
            // printf("failed to read width/height/max\n");
            exit(0);
        }
        // printf("Width=%d, Height=%d \nMaximum=%d\n",xdim,ydim,maxraw);
        
        image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
        getc(fp);
        
        line = (unsigned char *)malloc(sizeof(unsigned char)*xdim);
        for (j=0; j<ydim; j++) {
            fread(line, 1, xdim, fp);
            for (i=0; i<xdim; i++) {
                image[j*xdim+i] = line[i];
            }
        }
        free(line);
        
    }
    
    else if (c==2) {
        while ((c=fgetc(fp)) == '#')
            fgets(buf, 1024, fp);
        ungetc(c, fp);
        if (fscanf(fp, "%d%d%d", &xdim, &ydim, &maxraw) != 3) {
            // printf("failed to read width/height/max\n");
            exit(0);
        }
        // printf("Width=%d, Height=%d \nMaximum=%d,\n",xdim,ydim,maxraw);
        
        image = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim);
        getc(fp);
        
        for (j=0; j<ydim; j++)
            for (i=0; i<xdim; i++) {
                fscanf(fp, "%d", &val);
                image[j*xdim+i] = val;
            }
        
    }
    
    fclose(fp);
}

void WritePGM(FILE* fp)
{
    int i,j;
    
    fprintf(fp, "P5\n%d %d\n%d\n", xdim, ydim, 255);
    for (j=0; j<ydim; j++)
        for (i=0; i<xdim; i++) {
            fputc(image[j*xdim+i], fp);
        }
    
    fclose(fp);
    
}
