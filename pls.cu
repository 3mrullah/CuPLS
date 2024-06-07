#include <stdio.h>
#include <curand_kernel.h>
#include <curand.h>

#define MAX_LEN 256
#define MM 1000
#define NN 1000

#define T_NUM 1024
#define B_NUM 1024


__device__ int objFUNC(bool *sol, int m, int n, int Capacity, int max_length, int *H, int *items, int *weight, int *profit){ 
	int i,j,k, cap_val, cap_val_, sum_val;
	bool tempp[NN];
	
    for(j=0; j<n; j++)
        tempp[j] = 0;
	
	for(i=0; i<m; i++)
    {
        k = 0;
        if(sol[i])        
            while(items[k + i * max_length]!= -1)
            {
                tempp[items[k + i * max_length]] = 1;
                k++;
			}        
    }
		
    cap_val = sum_val =0;
    for(j=0; j<n; j++)    
        if(tempp[j])
            cap_val+= weight[j];
    	
    for(i=0; i<m; i++)	
		if(sol[i])		
			sum_val+= profit[i];	
	
    int a,b[MM];
	bool trial[NN], r_temp_sol[MM],r_temp[NN],r_trial[NN];	
	if(cap_val>Capacity)
    {		
		// Repair Start
        
        for(i=0; i<m; i++)
            r_temp_sol[i] = 0;

        for(j=0; j<n; j++)        
            r_temp[j] = r_trial[j] =0;            
        		
		for(k=0; k<max_length; k++)
            b[k] = 0;
		
		for(i=0; i<m; i++)
        {
            a = H[i];
            if(sol[a])
            {
                for(k=0; k<max_length; k++)
                    b[k] = items[k + a * max_length];
					
                for(k=0; k<max_length; k++)
                {
                    if(b[k]!=-1)                    
                        r_trial[b[k]] = 1;                    
                    else
                        break;
                }
                cap_val_ = 0;
                for(j=0; j<n; j++)
                {
                    if(r_trial[j])
                        cap_val_+= weight[j];
                }

                if(cap_val_<=Capacity)
                {
                    for(k=0; k<max_length; k++)
                    {
                        if(b[k]!=-1)                        
                            r_temp[b[k]] = 1;                        
                        else
                            break;
                    }
                    r_temp_sol[a]=1;
                }
                else
                {
                    for(j=0; j<n; j++)
                        r_trial[j] = r_temp[j];
                }
            }
        }
		// Repair End
		
		cap_val_ = 0;
        for(j=0; j<n; j++)        
            if(r_temp[j])
                cap_val_+= weight[j];
			
        sum_val = 0;
        for(i=0; i<m; i++)        
            if(r_temp_sol[i])            
                sum_val+= profit[i];           
		
		// Optimize
		for(i=0; i<m; i++)
            sol[i] = r_temp_sol[i];

        for(j=0; j<n; j++)
            tempp[j] = r_temp[j];

        for(j=0; j<n; j++)
            trial[j] = tempp[j];

        for(k=0; k<max_length; k++)
            b[k] = 0;
		
		for(i=0; i<m; i++)
        {
            a = H[i];
            if(sol[a] ==0)
            {
                for(k=0; k<max_length; k++)
                    b[k] = items[k + a * max_length];

                for(k=0; k<max_length; k++)
                {
                    if(b[k]!=-1)                    
                        trial[b[k]] = 1;                    
                    else
                        break;
                }
                cap_val_ = 0;
                for(j=0; j<n; j++)                
                    if(trial[j])
                        cap_val_+= weight[j];                

                if(cap_val_<=Capacity)
                {
                    for(k=0; k<max_length; k++)                    
						if(b[k]!=-1)
						{
							tempp[b[k]] = 1;
						}
						else
							break;                    
                    sol[a]=1;
                }
                else
                {
                    for(j=0; j<n; j++)
                        trial[j] = tempp[j];
                }
               

            }
        }
        // Optimization End !
	
		
	}
	else
	{
		// Optimize		
        for(j=0; j<n; j++)
            trial[j] = tempp[j];

        for(k=0; k<max_length; k++)
            b[k] = 0;

        for(i=0; i<m; i++)
        {
            a = H[i];
            if(sol[a] ==0)
            {
                for(k=0; k<max_length; k++)
                    b[k] = items[k + a * max_length];

                for(k=0; k<max_length; k++)                
                    if(b[k]!=-1)
                    {
                        trial[b[k]] = 1;
                    }
                    else
                        break;
                
                cap_val_ = 0;
                for(j=0; j<n; j++)                
                    if(trial[j])
                        cap_val_+= weight[j];                

                if(cap_val_<=Capacity)
                {
                    for(k=0; k<max_length; k++)                    
						if(b[k]!=-1)
						{
							tempp[b[k]] = 1;
						}
						else
							break;                    
                    sol[a]=1;
                }
                else
                {
                    for(j=0; j<n; j++)
                        trial[j] = tempp[j];
                }
            }
        }
        // Optimization end
		
		cap_val_ = 0;
		for(j=0; j<n; j++)		
			if(tempp[j])
				cap_val_+= weight[j];
        
		sum_val = 0;
		for(i=0; i<m; i++)		
			if(sol[i])			
				sum_val+= profit[i];		
		
	}
		
	return sum_val;
}


__device__ int penobjFUNC(bool *sol, int m, int n, int Capacity, int max_length, int *H, int *items, int *weight, int *profit, int L)
{ 
	int i,j,k;
	int R = 2;
	int n_x = 0;	
	
	int tempp[NN];
    for(j=0; j<n; j++)
        tempp[j] = 0;
	
	for(i=0; i<m; i++)
    {
        k = 0;
        if(sol[i])
        {
            while(items[k + i * max_length]!= -1)
            {
                tempp[items[k + i * max_length]] = 1;
                k++;
            }
        }
    }
		
	int cap_val = 0;
    for(j=0; j<n; j++)
    {
        if(tempp[j])
            cap_val+= weight[j];
    }
		
	if (cap_val>Capacity)
		n_x = cap_val - Capacity;
	else
		n_x = 0;
	
	int sum_val = 0;
	for(i=0; i<m; i++)
	{
		if(sol[i])
		{
			sum_val+= profit[i];
		}
	}
	
	if(n_x != 0)
	{
		if(sum_val>= L)
		{
			sum_val = L - R * n_x;
		}
		else
		{
			sum_val = sum_val - R * n_x;
		}
		
	}
    return sum_val;	
}

__global__ void sukpSOLVER(int *PARAM, int *p, int *w, int *HHH, int *itms, int *res, bool *sols)
{
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
	int b_tid = threadIdx.x;  // compute tid in a block
	int block_id = blockIdx.x; // Block id
	__shared__ int profits[T_NUM]; // For cost value
	__shared__ int thd_id; // Save the thd id that has a best value in a block
	__shared__ int best_profit_in_a_block ; // Save the best value in a block

	int temp, i, iter=0;
	bool solution[MM], solution_new[MM]; 
	
	curandState rndstate;
    curand_init(tid, 0, 0, &rndstate);
		
	for(i = 0; i < PARAM[0]; i++)	
		solution[i] = __float2uint_rn(curand_uniform(&rndstate)) & __float2uint_rn(curand_uniform(&rndstate));		
	
	profits[b_tid] = objFUNC(solution, PARAM[0], PARAM[1], PARAM[2], PARAM[3], HHH, itms, w, p);		
	__syncthreads();				

	// Find the best solution in the block ___Begin___
	if(b_tid==0)
	{		
		thd_id = 0;	
		best_profit_in_a_block = profits[0];
		for(i = 1; i < T_NUM; i++)
		{
			if(profits[i] > best_profit_in_a_block)
			{
				best_profit_in_a_block = profits[i];
				thd_id = i;				
			}
		}	
	} 	
	__syncthreads();
	
	if(b_tid==thd_id)
	{				
		res[block_id] = (profits[0] > res[block_id]) ? profits[0] : res[block_id];	
		
		for(i = 0; i < PARAM[0]; i++)
				sols[PARAM[0] * block_id + i] = solution[i];
					
	
	}
	// Find the best solution in the block ___End___
	
	
	__syncthreads();
	
	const int total_iter = max(PARAM[0], PARAM[1]) / 5;
	
	while(iter<100)
	{				
		
		
		/***** ##### Crossover ##### *****/
		for(i = 0; i < PARAM[0]; i++)		
			solution_new[i] = (solution[i] ^ ( ( __float2uint_rn(curand_uniform(&rndstate)) & 1) & ( (sols[PARAM[0] * block_id + i] ^ solution[i]) ) ));	
		/***** ##### Crossover ##### *****/
		
		temp = objFUNC(solution_new, PARAM[0], PARAM[1], PARAM[2], PARAM[3], HHH, itms, w, p);	

		if(temp > profits[b_tid])
		{
			profits[b_tid] = temp;
			for(i = 0; i < PARAM[0]; i++)
				solution[i] = solution_new[i];	
		}

		__syncthreads();

		// Find the best solution in the block ___Begin___
		if(b_tid==0)
		{		
			thd_id = 0;	
			best_profit_in_a_block = profits[0];
			for(i = 1; i < T_NUM; i++)
			{
				if(profits[i] > best_profit_in_a_block)
				{
					best_profit_in_a_block = profits[i];
					thd_id = i;				
				}
			}		
		}
		
		__syncthreads();		
		
		
		if(b_tid==thd_id)
		{
			if((best_profit_in_a_block > res[block_id]))
			{
				res[block_id] = best_profit_in_a_block;
				for(i = 0; i < PARAM[0]; i++)
						sols[PARAM[0] * block_id + i] = solution[i];				
			}	
			
		}
		// Find the best solution in the block ___End___
		
		iter++;
		__syncthreads();
	}
}

struct str
{
    float value;
    int index;
};

int cmp(const void *a, const void *b)
{
    struct str *a1 = (struct str *)a;
    struct str *a2 = (struct str *)b;
    if ((*a1).value > (*a2).value)
        return -1;
    else if ((*a1).value < (*a2).value)
        return 1;
    else
        return 0;
}

int main(int argc, char *argv[])
{
	
 	clock_t begin, end;
    double time_spent;
    begin = clock(); 

    if (argc != 2){
	  printf("Usage:\n\t%s input file\n", argv[0]);
	  exit(1);
    }
    char *fname= argv[1];
    int m, n, Capacity;
    int i,j;
    char infile[80],indirectory[]=".\\\\";
	
    strcpy(infile, fname);
    strcat(indirectory,infile);
    FILE *fp = fopen(fname, "r");
	
    char temp[MAX_LEN];
	
    fscanf(fp, "%s", &temp);
	
    char *token;
    const char s[2] = "=";
    token = strtok(temp, s);
    token = strtok(NULL, s);
    m = atoi(token);
    fscanf(fp, "%s", &temp);
    token = strtok(temp, s);
    token = strtok(NULL, s);
    n = atoi(token);
	
    fscanf(fp, "%s", &temp);
    fscanf(fp, "%s", &temp);
    token = strtok(temp, s);
    token = strtok(NULL, s);
    Capacity = atoi(token);

    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    int profit[m];

    for(i=0; i<m; i++)
    {
        fscanf(fp,"%s", &temp);
        profit[i] = atoi(temp);
    }

    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);
	
    int weight[n];

    for(i=0; i<n; i++)
    {
        fscanf(fp,"%s", &temp);
        weight[i] = atoi(temp);
    }

    fscanf(fp,"%s", &temp);
    fscanf(fp,"%s", &temp);

    bool rmatrix[m][n];

    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
        {
            fscanf(fp,"%s", &temp);
            rmatrix[i][j] = atoi(temp);
        }

    int freqs[n];
    for(j=0; j<n; j++)
         freqs[j] = 0;

    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
        {
            if(rmatrix[i][j])
                freqs[j]++;
        }


    double R[m], sum = 0.0;
    for(i=0; i<m; i++)
        R[m] = 0.0;

    for(i=0; i<m; i++)
    {
        sum = 0.0;
        for(j=0; j<n; j++)
        {
            if(rmatrix[i][j])
            {
                sum += weight[j]/ (double)freqs[j];
            }


        }

        R[i] = profit[i] / (double)sum;
    }

    struct str H[m]; 
    for (int i = 0; i < m; i++)
    {
        H[i].value = R[i];
        H[i].index = i;
    }

    //sort objects array according to value maybe using qsort
    qsort(H, m, sizeof(H[0]), cmp);
    int HH[m];
    for(int i = 0; i < m; i++)
        HH[i] = H[i].index;
	
    int max_length = 0, mal_length = 0, k = 0;

    for(i=0; i<m; i++)
    {
        mal_length = 0.0;
        for(j=0; j<n; j++)
        {
            if(rmatrix[i][j])
                mal_length++;
        }

        if(mal_length>max_length)
            max_length = mal_length;
    }

    int items[m][max_length];

    for(i=0; i<m; i++)
        for(j=0; j<max_length; j++)
            items[i][j] = -1;

    for(i=0; i<m; i++)
    {
        k = 0;
        for(j=0; j<n; j++)
        {
            if(rmatrix[i][j])
            {
                items[i][k] = j;
                k++;
            }
        }
    }
	
	// File read completed. Problem is ready to handle it.
	
	int *param; //m, n, Capacity, max_length
	size_t bytes = 4 * sizeof(int);
	
	cudaMallocManaged(&param, bytes);
	
	int *p, *w, *HHH; 
	size_t s_profit = m * sizeof(int);
	size_t s_weight = n * sizeof(int);
	size_t s_HHH = m * sizeof(int);
	cudaMallocManaged(&p, s_profit);
	cudaMallocManaged(&w, s_weight);	
	cudaMallocManaged(&HHH, s_HHH);		
	
	int *itms;
	size_t s_itms = m * max_length * sizeof(int);
	cudaMallocManaged(&itms, s_itms);	
		
	int threadsPerBlock = T_NUM;
    int blocksPerGrid = B_NUM;
	
	
	// Best objective values for each block
	int *res;
	size_t res_size = (blocksPerGrid) * sizeof(int);
	cudaMallocManaged(&res, res_size);
	
	bool *sols;
	size_t sols_size = (blocksPerGrid) * m * sizeof(bool);
	cudaMallocManaged(&sols, sols_size);

	// Get the device ID for prefetching calls
	int id = cudaGetDevice(&id);
	
	cudaMemAdvise(param, bytes, cudaMemAdviseSetPreferredLocation, cudaGetDevice(&id));
	cudaMemAdvise(p, s_profit, cudaMemAdviseSetPreferredLocation, cudaGetDevice(&id));
	cudaMemAdvise(w, s_weight, cudaMemAdviseSetPreferredLocation, cudaGetDevice(&id));
	cudaMemAdvise(HHH, s_HHH, cudaMemAdviseSetPreferredLocation, cudaGetDevice(&id));
	cudaMemAdvise(itms, s_itms, cudaMemAdviseSetPreferredLocation, cudaGetDevice(&id));
	cudaMemPrefetchAsync(res, res_size, id);
	cudaMemPrefetchAsync(sols, sols_size, id);
	
	param[0] = m;
	param[1] = n;
	param[2] = Capacity;
	param[3] = max_length;
	
	for(i=0; i<m; i++)
		p[i] = profit[i];
	
	for(i=0; i<n; i++)
		w[i] = weight[i];
	
	for(i=0; i<m; i++)
		HHH[i] = HH[i];
	
	
	for(i=0; i<m; i++)   
        for(j=0; j<n; j++)
            itms[j +  i*max_length] = items[i][j];
        
    for(i=0; i<blocksPerGrid; i++) 
		res[i] = 0;
	
	cudaMemAdvise(param, bytes, cudaMemAdviseSetReadMostly, cudaGetDevice(&id));
	cudaMemPrefetchAsync(param, bytes, cudaGetDevice(&id));
	
	cudaMemAdvise(p, s_profit, cudaMemAdviseSetReadMostly, cudaGetDevice(&id));
	cudaMemPrefetchAsync(p, s_profit, cudaGetDevice(&id));
	
	cudaMemAdvise(w, s_weight, cudaMemAdviseSetReadMostly, cudaGetDevice(&id));
	cudaMemPrefetchAsync(w, s_weight, cudaGetDevice(&id));
	
	cudaMemAdvise(HHH, s_HHH, cudaMemAdviseSetReadMostly, cudaGetDevice(&id));
	cudaMemPrefetchAsync(HHH, s_HHH, cudaGetDevice(&id));
	
	cudaMemAdvise(itms, s_itms, cudaMemAdviseSetReadMostly, cudaGetDevice(&id));
	cudaMemPrefetchAsync(itms, s_itms, cudaGetDevice(&id));
	

	// Call CUDA kernel
	printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);
    sukpSOLVER<<<blocksPerGrid, threadsPerBlock>>>(param, p, w, HHH, itms, res, sols);
	

	cudaDeviceSynchronize();

	
	cudaMemPrefetchAsync(res, res_size, cudaCpuDeviceId);
	cudaMemPrefetchAsync(sols, sols_size, id);
	
	j = res[0];
	for (int i = 1; i < blocksPerGrid; i++)
    {	
		if(res[i] > j)
		{
			j = res[i];
			k = i;
		}			
        
    }
	printf("\n # %s \t: %d in block %d #\n",argv[1], j, k);

	
	cudaFree(param);
	cudaFree(p);
	cudaFree(w);
	cudaFree(HHH);
	cudaFree(itms);
	cudaFree(res);
	cudaFree(sols);
	
	end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\nDone \t %.3f \t\n", time_spent);

    return 0;
}