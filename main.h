#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <stdarg.h>
#include<math.h>
#include<time.h>
#define PI 3.1415926


//#define UNUM 2
#include <assert.h>

#define NUMBER_MAX 20000   //Greater than datasetNum
#define ASIZE 9733
//#define ASIZE 6073
//#define ASIZE 6997
//#define ASIZE 4937
//#define ASIZE 7933
//#define ASIZE 5737
//#define ASIZE 9973
#define TEMP_MAX 30

#define DATA_DIMENSION 128 
#define dimension 128	// dimension
#define BLOOM_L 4	// the number of bloom(parameter L in lsh)
#define k_hfun 8	// the number of h_function(parameter K in lsh)
#define W 1 		// The default value for lsh algorithm parameter W


#define SETBIT(a, n) (a[n/CHAR_BIT] |= (1<<(n%CHAR_BIT)))  // set a bit to 1
#define GETBIT(a, n) (a[n/CHAR_BIT] & (1<<(n%CHAR_BIT)))   // get a bit ,1 or 0

#define SQR(a) ((a) * (a))

#define CHAR_BIT 8
// 4294967291 = 2^32-5
#define UH_PRIME_DEFAULT 4294967291U
// 2^29
#define MAX_HASH_RND 536870912U     // lsh_r(1,MAX_HASH_RND)

//unsigned **lsh_r; 	//used to calculate the gindex of lsh, gindex=((lsh_r*a)mod prime)mod tableSize
float **dataset; 	//data set
int dataNum; 		//the size of dataset
float **queryset;
float **project_vector;
float *feature_value;
//float **weight;



#define SQRT sqrt
#define ABS fabs
#define LOG log
#define COS cos
#define CEIL(x) ((int)(ceil(x)))
#define POW(x, y) (pow(x,y))
#define ERFC erfc
#define EXP exp
#define ERF erf

// typedef struct
// {
        //int nfuncs;
        // float **para_a;
        // float *para_b;
// } BLOOM_M;

typedef struct
{
	int asize; 		// the size of a
	unsigned char *a;
	int nfuncs; 	// there are nfuncs function in a bloom(= K)
	float **para_a;// parameter a for lsh
	float *para_b; // parameter b for lsh
        //BLOOM_M **bloom_n;
	unsigned *para_r;//used to calculate the gindex of lsh, gindex=((lsh_r*a)mod prime)mod tableSize
} BLOOM;

typedef struct
{
       unsigned int *index;      //the hashed value of points at each hash table
}DATA_INDEX;

typedef struct
{
       int index;
       int *num;
       int m;
} SIMILAR;

float genGaussianRandom();
float genUniformRandom(float rangeStart, float rangeEnd);
int genRandomInt(int rangeStart, int rangeEnd);

// read datasets
int getds(char *s, int dsNum);
//read querysets
int getqs(char *s, int dsNum);

int getvs(char *s, int dsNum);

int getfs(char *s, int dsNum);

int getindex(BLOOM *bloom, unsigned *temp);					// get g(v) = ((lsh_r*h(v))mod prime)mod tableSize
unsigned int *getvector(BLOOM* bloom, float *f,float R, int n); 	// get temp[0...k-1] = [h1(v)...hk(v)]

// init L bloom filter
BLOOM **bloom_create(int nfuncs, int size);

// set L bloom filter using dataset						
DATA_INDEX **bloom_set(BLOOM **bloom, float **dataset, int dataNum,float R);


int * distance(DATA_INDEX *sourcepoint, unsigned int index,int datasetNum,float *s,float R,int flag[],float **dataset,int *real_num,int *true_num);

float distanceSqr(float *p1,float *p2,float threshold);

// check whether point s is similar to dataset
float *bloom_check_similar(BLOOM **bloom, float *s,float R,DATA_INDEX **dataindex,int dataNum,float **dataset);	

// destory L bloom filter
int bloom_destroy(BLOOM **bloom);	

int data_destroy(int datasetNum,int queryNum,int vectorNum);	

//int bloom_m_destroy(BLOOM_M **bloom_m,int m,int nfuncs);

int sourcepoint_destroy(DATA_INDEX **sourcepoint,int datasetNum);	

float *exact(FILE *fp,int datasetNum,float **dataset,int queryNum,float **queryset,float R);				
