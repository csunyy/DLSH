#include "main.h"

//read dataset
int getds(char *s, int dsNum)
{
	FILE *fp;
	int i, d;
	fp = fopen(s, "r");
	dataset = (float **)malloc(dsNum * sizeof(float*));
	for(i = 0; i < dsNum; i++)
	{
		dataset[i] = (float*)malloc(DATA_DIMENSION * sizeof(float));
		for(d = 0; d < DATA_DIMENSION; d++)
		{
			fscanf(fp, "%f", &dataset[i][d]);
		}
	}
	fclose(fp);
	return 1;
}

//read queryset
int getqs(char *s, int dsNum)
{
	FILE *fp;
	int i, d;
	fp = fopen(s, "r");
	queryset = (float **)malloc(dsNum * sizeof(float*));
	for(i = 0; i < dsNum; i++)
	{
		queryset[i] = (float*)malloc(DATA_DIMENSION * sizeof(float));
		for(d = 0; d < DATA_DIMENSION; d++)
		{
			fscanf(fp, "%f", &queryset[i][d]);
		}
	}
	fclose(fp);
	return 1;
}

//get the project_vector set
int getvs(char *s, int dsNum)
{
	FILE *fp;
	int i, d;
	fp = fopen(s, "r");
	project_vector = (float **)malloc(dsNum * sizeof(float*));
	for(i = 0; i < dsNum; i++)
	{
		project_vector[i] = (float*)malloc(DATA_DIMENSION * sizeof(float));
		for(d = 0; d < DATA_DIMENSION; d++)
		{
			fscanf(fp, "%f", &project_vector[i][d]);
		}
	}
	fclose(fp);
	return 1;
}

// get the feature_value set
/* int getfs(char *s, int dsNum)
{
	FILE *fp;
	float temp;
	int i, j, k;
	fp = fopen(s, "r");
	feature_value = (float*)malloc(dsNum * sizeof(float));
	weight = (float*)malloc(BLOOM_L * sizeof(float));
	for(i = 0; i < dsNum; i ++)
	{
		fscanf(fp, "%f", &feature_value[i]);
	}
	for(j = 0; j < BLOOM_L; j ++)
	{
		for(k = 0; k < k_hfun; k ++)
		{
			temp = 0.0;
			temp += feature_value[k + j*k_hfun];
		}
		weight[j] = temp;
	}
	fclose(fp);
	return 1;
} */

float genUniformRandom(float rangeStart, float rangeEnd)
{
	float r;
	do
	{
		r = rangeStart + ((rangeEnd - rangeStart) * rand() / (float) RAND_MAX);
	} while (r < rangeStart || r > rangeEnd);
	return r;
}

float genGaussianRandom()
{
	// Use Box-Muller transform to generate a point from normal
	// distribution.
	float x1, x2, z;
	do
	{
		x1 = genUniformRandom(0.0, 1.0);
	} while (x1 == 0); // cannot take log of 0.
	x2 = genUniformRandom(0.0, 1.0);
	z = sqrt(-2.0 * log(x1)) * cos(2.0 * PI * x2);
	return z;
}

int genRandomInt(int rangeStart, int rangeEnd)
{
	int r;
	do
	{
		r = rangeStart
				+ (int) ((rangeEnd - rangeStart + 1.0) * rand() / (RAND_MAX + 1.0));
	} while (r < rangeStart || r > rangeEnd);
	return r;
}

//init L bloom filters.
BLOOM **bloom_create(int nfuncs, int size)
{
	BLOOM **bloom;
	int i,j,k,l;
	// init g(v) functions
	// if (!(lsh_r = (unsigned **) malloc(BLOOM_L * sizeof(unsigned))))
		// return NULL;
	// for (k = 0; k < nfuncs; k++)
	// {
		// lsh_r[k] = genRandomInt(1, MAX_HASH_RND);
	// }
	
	// init L bloom filters
	if (!(bloom = (BLOOM **) malloc(BLOOM_L * sizeof(BLOOM *))))
		return NULL;
	for (i = 0; i < BLOOM_L; i++)
	{

		// first malloc space
		if (!(bloom[i] = (BLOOM *) malloc(sizeof(BLOOM))))
			return NULL;

		if (!(bloom[i]->a = calloc((size + CHAR_BIT - 1) / CHAR_BIT,
				sizeof(char))))
		{
			free(bloom[i]);
			free(bloom);
			return NULL;
		}
		bloom[i]->asize = size;
		bloom[i]->nfuncs = nfuncs;
		bloom[i]->para_r = (unsigned *)malloc(nfuncs * sizeof(unsigned));
		bloom[i]->para_a = (float **)malloc(nfuncs * sizeof(float *));
		for(j = 0; j < nfuncs; j ++)
		{
			bloom[i]->para_a[j] = (float *)malloc(dimension *sizeof(float));
		}
		bloom[i]->para_b = (float*)malloc(nfuncs*sizeof(float));
		
		for(k = 0; k < nfuncs; k ++)
		{
			for(l = 0; l < dimension; l ++)
			{
				bloom[i]->para_a[k][l] = project_vector[i*k_hfun + k][l];
			}
			bloom[i]->para_b[k] = genUniformRandom(0, W);
			bloom[i]->para_r[k] = genRandomInt(1, MAX_HASH_RND);
			//bloom[i]->para_r[k] = feature_value[k + i*k_hfun] * MAX_HASH_RND/weight[i];
		}
	}
	return bloom;
}

int getindex(BLOOM* bloom, unsigned *temp)
{
	int i;
	unsigned int index = 0;
	float item = 0.0;
	for (i = 0; i < bloom->nfuncs; i++)
	{
		index += temp[i] * bloom->para_r[i];
	}
	index %= UH_PRIME_DEFAULT;
	index %= bloom->asize;

	return index; 	//gIndex = g(v) =((lsh_r*h(v))mod prime)mod tableSize
	free(temp);
}

// in this function, we get h(v) and store to temp
unsigned int *getvector(BLOOM* bloom, float *f,float R, int n) 
{
	int i,j,k;
    float result = 0;
	unsigned int *vec = (unsigned int *)malloc(bloom->nfuncs*sizeof(unsigned int));
	memset(vec,0,bloom->nfuncs * sizeof(unsigned int));
    for (i = 0; i < bloom->nfuncs; i++)
    {
		result = bloom->para_b[i];
	    for (k = 0; k < dimension; k++)
	    {
		    result += f[k] * ((*(bloom->para_a[i] + k))/R);
	    }
		result/=W/pow(2,n);
        vec[i] = (unsigned) (floor(result)); // h(v) = (a.v+b)/w
	}
       return vec;
}

// set L bloom filter using dataset
DATA_INDEX **bloom_set(BLOOM **bloom, float **dataset, int datasetNum,float R) 
{
	int i, j,k,m;
	unsigned int index, *temp;
        DATA_INDEX **sourcepoint;
        
        if (!(sourcepoint = (DATA_INDEX **) malloc(BLOOM_L * sizeof(DATA_INDEX *))))
		return NULL;
	for (i = 0; i < BLOOM_L; i++)	
	{
        if (!(sourcepoint[i] = (DATA_INDEX *) malloc(sizeof(DATA_INDEX))))
		{
			return NULL;
		}
        sourcepoint[i]->index = (unsigned int *)malloc(datasetNum * sizeof(unsigned int));
		for (j = 0; j < datasetNum; j++)
		{
            temp = getvector(bloom[i], dataset[j],R,i);
			index = getindex(bloom[i], temp);
            sourcepoint[i]->index[j] = index;
			SETBIT(bloom[i]->a, index);
		}
	}
	return sourcepoint;
}

float distanceSqr(float *p1,float *p2,float threshold)
{
      float result = 0;
      int d;
	for (d = 0; d < dimension; d++)
	{
		result += (p1[d] - p2[d]) * (p1[d] - p2[d]);
	}
         result = SQRT(result);
        if(result <= threshold)
        {
        	 return result;
        }
        else
	{
		return -1;
	}
}

int * distance(DATA_INDEX *sourcepoint, unsigned int index,int datasetNum,float *s,float R,int flag[],float **dataset,int *real_num,int *true_num)
{
      int j = 0;
      float result;
      for(j = 0;j < datasetNum; j++)
      {  
	      if(sourcepoint->index[j] == index)
            {
				if(flag[j] == 0)
				{
					*real_num += 1;//result_num
					result = distanceSqr(dataset[j],s,R);
					if(result >= 0)
					{
	//			printf("similar num.%d,  distance = %f\n",j,result);
						*true_num += 1;//result_true_num
					}
					flag[j] += 1;//mark the item compared
				}                    
            }
      }
      return NULL;
}

int sourcepoint_destroy(DATA_INDEX **sourcepoint,int datasetNum)
{
        int i,j;
        for(i = 0;i < BLOOM_L;i ++)
        {
                free(sourcepoint[i]->index);
                free(sourcepoint[i]);
        }
        free(sourcepoint);
        return 0;
}


int bloom_destroy(BLOOM **bloom)
{
	int i, j, k;
	for (i = 0; i < BLOOM_L; i++)
	{
		free(bloom[i]->a);
		for(j = 0;j < k_hfun; j++)
		{
			free(bloom[i]->para_a[j]);
		}
		free(bloom[i]->para_a);
		free(bloom[i]->para_b);
		free(bloom[i]->para_r);
		free(bloom[i]);
	}
	free(bloom);
	printf("*** destroy bloom done! ***\n");
	return 0;
}

int data_destroy(int datasetNum,int queryNum, int vectorNum)
{

    int i = 0;
    for (i = 0; i < datasetNum; i++)
	{
		free(dataset[i]);
	}

	free(dataset);
    for (i = 0; i < queryNum; i++)
	{
		free(queryset[i]);
	}
		free(queryset);
	for(i = 0; i < vectorNum; i ++)
	{
		free(project_vector[i]);
	}
		free(project_vector);
    return 0;
}


float *bloom_check_similar(BLOOM **bloom, float *s,float R,DATA_INDEX **sourcepoint,int datasetNum,float **dataset)
{
	int i,j,k;
  	int *flag1,*flag2;
	flag2 = (int *)malloc(datasetNum * sizeof(int));//allocate a flag for each data
	memset(flag2,0,datasetNum * sizeof(int));//initial 0
	flag1 = (int *)malloc(datasetNum * sizeof(int));
	memset(flag1,0,datasetNum * sizeof(int));
    float result = 0;
	unsigned int *temp, index = 0;
	float *final_result;
	final_result = (float *) malloc(3*sizeof(int));
	int total_num = 0,real_num = 0,true_num = 0;
	float similar_weight = 0.0;
    for (i = 0; i < BLOOM_L; i++)
    {
		temp = getvector(bloom[i], s,R,i);
		index = getindex(bloom[i], temp);
		if(GETBIT(bloom[i]->a, index))
		{
			if(i == 0)
			{
				for(j = 0; j < datasetNum; j ++)
				{
					if(flag2 == 0)
					{
						real_num += 1;//result_num
							result = distanceSqr(dataset[j],s,R);
							if(result >= 0)
							{
								similar_weight += 1.0/(pow(result,2)+10000);
								true_num += 1;
							}
							flag2[j] = 1;//mark the item compared
					}
				}
			}
			//distance(sourcepoint[i],index,datasetNum,s,R,flag2,dataset,&real_num,&true_num);
			for(j = 0; j < datasetNum; j ++)
			{
				if(sourcepoint[i]->index[j] == index)
				{
					flag1[j] += 1;
					if(flag1[j] >= ((BLOOM_L-1) *2)/3)
					{
						if(flag2[j] == 0)
						{
							real_num += 1;//result_num
							result = distanceSqr(dataset[j],s,R);
							if(result >= 0)
							{
								similar_weight += 1.0/(pow(result,2)+10000);
								true_num += 1;
							}
							flag2[j] = 1;//mark the item compared
						}
					}
					
				}
			}
		}
    }
	final_result[0] = real_num;
	final_result[1] = true_num;
	final_result[2] = similar_weight;
	free(flag2);
	free(flag1);
    return final_result;
}


float *exact(FILE *fp,int datasetNum,float **dataset,int queryNum,float **queryset,float R)
{
        int i,j,k;
        float result = 0;
		float *total;
		total = (float*)malloc(2*sizeof(float));
	int total_num = 0;
	float total_weight = 0.0;
	//weight = (float **)malloc(queryNum * sizeof(float*));
       	for(i = 0;i < queryNum;i ++)
        {
              k = 0;
			  //weight[i] = (float *)malloc(datasetNum * sizeof(float));
              for(j = 0;j < datasetNum; j++)
              {
                    result = distanceSqr(queryset[i],dataset[j],R);
                    if(result >= 0)
                    {
                          //  printf("\tsimilar point num.%d\tdistance:%f\n",j,result);
                            k++;
							//weight[i][j] = 1/pow(result,2);
							//total_weight += weight[i][j];
							total_weight += 1.0/(pow(result,2)+10000);
                    }
					// if(result == 0)
					// {
						// k++;
						// total_weight +=1.0;
					// }
              }
        //      printf("Query point %d: found %d right points.\n",i,k);
			  fprintf(fp,"Query point %d: found %d right points.\n",i,k);
	      total_num += k;
        }
		total[0] = total_num;
		total[1] = total_weight;
	return total;
}

int main(int nargs, char **args) {
	if(nargs < 8) {
		printf("Usage: %s #nDatasetPoints #nQueryPoints R datasetFilename querysetFilename vectorFilename successProbability\n", args[0]);
		return 0;
	}
	clock_t start,end;
	BLOOM **bloom;
	int i, j, d,m;
	float *final_result;
	float *total;
	int datasetNum, queryNum;
	float R, P;
	float total_weight = 0.0;
	float similar_weight = 0.0;
	char fname[100];
	//float *query;
	int exact_num = 0,similar_num = 0,result_num = 0;
	float recall, precision, F1_score;
    DATA_INDEX **sourcepoint;
	datasetNum = atoi(args[1]);
	queryNum = atoi(args[2]);
	R = atof(args[3]);
	getds(args[4], datasetNum);
    getqs(args[5], queryNum);
	getvs(args[6], BLOOM_L * k_hfun);
	//getfs(args[7], DATA_DIMENSION);
    P = atof(args[7]); 
	//final_result = (int *)malloc(2*sizeof(int));
	sprintf(fname, "%s_K_%d_L_%d_R_%s_Wpca_result.txt", args[4], k_hfun, BLOOM_L, args[3]);
	// fname = args[4];
	// fname += "_R_";
	// fname += args[3];
	// fname += "_pca_result.txt";
	FILE *fp;
	if(!(fp=fopen(fname,"w+"))){
		printf("Failure to open result.txt!\n");
	}
      
    if(!(bloom = bloom_create(k_hfun,ASIZE)))
	{
		printf("ERROR: Could not create bloom filter\n");
		return -1;
    }
  //  printf("Create bloom succeed!\n");
    sourcepoint=bloom_set(bloom, dataset, datasetNum, R);
  //  printf("Set bloom succeed!\n");
//query
//	printf("\t**********exact query**********\n");
		fprintf(fp,"\t**********exact query**********\n");
		//exact_num = exact(fp,datasetNum,dataset,queryNum,queryset,R,total_weight);
		total = exact(fp,datasetNum,dataset,queryNum,queryset,R);
		exact_num = total[0];
		total_weight = total[1];
		printf("total_weight = %f\n",total_weight);
		
	//	printf("\t**********similar query**********\n");
		fprintf(fp,"\t**********similar query**********\n");
		start = clock();
     	for(i = 0; i < queryNum; i++)
		{
			final_result = bloom_check_similar(bloom, queryset[i],R,sourcepoint,datasetNum,dataset);
			//similar_weight = bloom_check_similar(bloom, queryset[i],R,sourcepoint,datasetNum,dataset,similar_weight);
			//fprintf(fp,"Query Point %d: found %d right points among %d similar points\n", i, final_result[1], final_result[0]);
			//printf("Query Point %d: found %d right points among %d similar points\n", i, final_result[1], final_result[0]);
			similar_num += final_result[1];
			result_num += final_result[0];
			similar_weight += final_result[2];
			free(final_result);
		}
		end = clock();
		printf("The ANN query time is: %f seconds\n",((double)(end - start))/CLOCKS_PER_SEC);
		printf("similar_weight = %f\n",similar_weight);
		fprintf(fp,"\n\n");
		printf("\n\n");
	
		fprintf(fp,"similar points: %d, right points(find): %d, exact points: %d\n",result_num,similar_num,exact_num);
		printf("similar points: %d, right points(find): %d, exact points: %d\n",result_num,similar_num,exact_num);
		//recall = similar_num*1.0/exact_num;
		recall = similar_weight/total_weight;
		precision = similar_num*1.0/result_num;
		F1_score = 2*recall*precision*1.0/(recall + precision);
		fprintf(fp,"when k = %d, L = %d, recall = %f, precision = %f, F1_core = %f\n",k_hfun,BLOOM_L,recall,precision,F1_score);
		printf("when k = %d, L = %d, recall = %f, precision = %f, F1_core = %f\n",k_hfun,BLOOM_L,recall,precision,F1_score);
		//printf("when k = %d, L = %d, R = %f, recall = %f\n",k_hfun,R,BLOOM_L,recall);
		fclose(fp);

		bloom_destroy(bloom);
		sourcepoint_destroy(sourcepoint,datasetNum);
		data_destroy(datasetNum,queryNum,BLOOM_L*k_hfun);
		printf("\t******over******\n");

		return 0;
}
