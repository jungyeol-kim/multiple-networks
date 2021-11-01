
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <fstream>
#include <iostream> 
#include <vector>
#define SWAP(aa,bb,cc) ((cc)=(aa),(aa)=(bb),(bb)=(cc))

using namespace std;

void QuickSort(int list[], int index[], int sort[], int left, int right);
void rev_QuickSort(int list[], int index[], int sort[], int left, int right);
int partition(int list[], int index[], int sort[], int left, int right);
int rev_partition(int list[], int index[], int sort[], int left, int right);

int main(int argc, char*argv[]){
	
	int total=10000;
	int m=1;
	double sum1=0, sum2=0, sum_temp1, sum_temp2;
	int k1,k2;
	double p1,p2;
	double a1, b1, a2=-0.05, b2;
	
	int *nl11, *nl22, *nl1, *nl2;
	int **ll11, **ll22, **ll;	
	
	double degree;
	double *pk;
	
	int *rank1, *rank2, *index1, *index2, *temp;
	
	double pa, pb, pc, pd, pe, pearson_temp, p_sum, pearson;
	double sa, sb, sc, sd, se, spearman_temp, s_sum, spearman;
	int x,y;
	
	int ensenble=0;
	int ensenble_num=100;
	
	double *error, *error2;
	double errorbar, errorbar2;
	
	double alpha=100, beta=0;
	
	nl22=(int *)malloc((size_t)(total*sizeof(int)));
	ll22=(int **)malloc((size_t)(total*sizeof(int *)));
	nl11=(int *)malloc((size_t)(total*sizeof(int)));
	ll11=(int **)malloc((size_t)(total*sizeof(int *)));
	
	nl1=(int *)malloc((size_t)(total*sizeof(int)));
	nl2=(int *)malloc((size_t)(total*sizeof(int)));

	
	rank1=(int *)malloc((size_t)(total*sizeof(int)));
	index1=(int *)malloc((size_t)(total*sizeof(int)));
	rank2=(int *)malloc((size_t)(total*sizeof(int)));
	index2=(int *)malloc((size_t)(total*sizeof(int)));

	temp=(int *)malloc((size_t)(total*sizeof(int)));

	
	pk=(double *)malloc((size_t)(total*sizeof(double)));
	
	error=(double *)malloc((size_t)(ensenble_num*sizeof(double)));
	error2=(double *)malloc((size_t)(ensenble_num*sizeof(double)));

	
	
	srand48(time(NULL));
	
	
	ofstream output1("k+100,kk_s_pearson");
	ofstream output2("k+100,kk_s_spearman");
	
	
	

	for(int v=0; v<21; v++)
	{
		a1=1;
		b1=1.-a1;
		a2+=0.05;
		b2=1.-a2;
		

		p_sum=0; s_sum=0; ensenble=0; errorbar=0; errorbar2=0;
		
		double a_sum=0, a_summ=0;
		
		for(int ens=0; ens<ensenble_num; ens++)
		{
			for(int i=0; i<total; i++){
				nl11[i]=0;
				nl22[i]=0;
				
				index1[i]=i;
				rank1[i]=i;
				index2[i]=i;
				rank2[i]=i;
				
			}		
			for(int i=0; i<=m; i++){
				nl11[i]=m;
				nl22[i]=m;
				ll11[i]=(int*)realloc(ll11[i],(size_t)(m*sizeof(int)));
				ll22[i]=(int*)realloc(ll22[i],(size_t)(m*sizeof(int)));

				int k=0;
				for(int j=0; j<m; j++){
					if(i==k){	
						k++;
						ll11[i][j]=k;
						ll22[i][j]=k;
					}
					if(i!=k){
						ll11[i][j]=k;
						ll22[i][j]=k;
					}
					k++;
				}
			}
			sum1=0; sum2=0;
			for(int i=0; i<=m; i++){
				sum1+=a1*(nl11[i]+alpha)+b1*(nl22[i]+alpha);
				sum2+=b2*(nl11[i]+beta)+a2*(nl22[i]+beta);
			}
			
			for(int i=m+1; i<total; i++){
				for(int j=0; j<m; j++){
					p1=drand48()*sum1; k1=0; sum_temp1=0;
					do{
						sum_temp1+=a1*(nl11[k1]+alpha)+b1*(nl22[k1]+alpha);
						k1++;
					} while(sum_temp1<p1);
					
					
					p2=drand48()*sum2; k2=0; sum_temp2=0;
					do{
						sum_temp2+=b2*(nl11[k2]+beta)+a2*(nl22[k2]+beta);
						k2++;
					} while(sum_temp2<p2);
					
					k1=k1-1;
					nl11[k1]++;
					nl11[i]++;
					ll11[k1]=(int*)realloc(ll11[k1],(size_t)((nl11[k1])*sizeof(int)));
					ll11[i]=(int*)realloc(ll11[i],(size_t)((nl11[i])*sizeof(int)));		
					ll11[k1][nl11[k1]-1]=i;
					ll11[i][nl11[i]-1]=k1;
					
										
					k2=k2-1;
					nl22[k2]++;
					nl22[i]++;
					ll22[k2]=(int*)realloc(ll22[k2],(size_t)((nl22[k2])*sizeof(int)));
					ll22[i]=(int*)realloc(ll22[i],(size_t)((nl22[i])*sizeof(int)));		
					ll22[k2][nl22[k2]-1]=i;
					ll22[i][nl22[i]-1]=k2;
			
					sum1=0., sum2=0.;
					for(int p=0; p<=i; p++){
						sum1+=a1*(nl11[p]+alpha)+b1*(nl22[p]+alpha);
						sum2+=b2*(nl11[p]+beta)+a2*(nl22[p]+beta);
					}
					//printf("%f %f %f %f\n",p1,sum1,p2,sum2);
				/*	
					sum1=sum1-(a*((nl11[i]-1)+alpha)+b*((nl22[i]-1)+alpha))-(a*((nl11[k1]-1)+alpha)+b*((nl22[k1]-1)+alpha))
					+(a*(nl11[i]+alpha)+b*(nl22[i]+alpha))+(a*(nl11[k1]+alpha)+b*(nl22[k1]+alpha));
					
					sum2=sum2-(b*((nl11[i]-1)+alpha)+a*((nl22[i]-1)+alpha))-(b*((nl11[k2]-1)+alpha)+a*((nl22[k2]-1)+alpha))
					+(b*(nl11[i]+alpha)+a*(nl22[i]+alpha))+(b*(nl11[k2]+alpha)+a*(nl22[k2]+alpha));
				*/
				}
			}
			
			
			pa=0; pb=0; pc=0; pd=0; pe=0; 
			for(int i=0;i<total;i++){ 
				x=nl11[i]; y=nl22[i];
				pa+=x; pb+=y; pc+=x*x; pd+=y*y; pe+=x*y; 
			}
			pa/=total; pb/=total; pc/=total; pd/=total; pe/=total;
			pearson_temp=(pe-pa*pb)/sqrt((pc-pa*pa)*(pd-pb*pb));
			//pearson_temp=(pe-pa*pb)/(bb-4.);
			error[ens]=pearson_temp;
			p_sum+=pearson_temp;
			
		//	for(int i=0; i<total; i++){
		//	printf("%d %d\n",nl11[i],nl22[i]);
		//	}
		
			for(int i=0;i<total;i++){ 
				nl1[i]=nl11[i];
				nl2[i]=nl22[i];
			}
			rev_QuickSort(nl1, index1, rank1, 0, total-1);
			rev_QuickSort(nl2, index2, rank2, 0, total-1);

			sa=0; sb=0; sc=0; sd=0; se=0; 
			for(int i=0;i<total;i++){ 
				x=rank1[i]+1; y=rank2[i]+1;
				sa+=x; sb+=y; sc+=x*x; sd+=y*y; se+=x*y; 
			}
			sa/=total; sb/=total; sc/=total; sd/=total; se/=total;
			spearman_temp=(se-sa*sb)/sqrt((sc-sa*sa)*(sd-sb*sb));
			error2[ens]=spearman_temp;
			s_sum+=spearman_temp;
			
			ensenble++;
			
		}
		pearson=p_sum/ensenble;
		spearman=s_sum/ensenble;
		
		for(int i=0; i<ensenble_num; i++){
			errorbar+=pow((error[i]-pearson),2)/ensenble;
			errorbar2+=pow((error2[i]-spearman),2)/ensenble;
		}
	
		errorbar=sqrt(errorbar);
		errorbar2=sqrt(errorbar2);

		
		
		output1 << b2 << " " << pearson << " " << errorbar << endl;
		output2 << b2 << " " << spearman << " " << errorbar2 << endl;

	}
			
		
		
	ofstream output3("c_pdf");
	double tempdouble1;
	
	for(int i=0; i<total; i++){
		pk[i]=0;
	}  
	
	for(int i=0; i<total; i++){
		pk[nl11[i]]++;
	}
	
	tempdouble1=0;
	for(int i=total-1; i>=0; i--){
		if(pk[i]){
			tempdouble1+=pk[i]/(double)total;
			output3 << i << " " << tempdouble1 << " " << pk[i]/(double)total << endl;
		}
	}
	
	free(nl1); free(nl2); free(nl11); free(nl22);
	free(ll); free(ll11); free(ll22); free(pk);
	
}

void QuickSort(int list[], int index[], int sort[], int left, int right){
	if(left<right){
		int q = partition(list, index, sort, left, right);
		QuickSort(list, index, sort, left, q-1);
		QuickSort(list, index, sort, q+1, right);
	}
}

int partition(int list[], int index[], int sort[], int left, int right){
	int pivot, temp;
	int low, high, median;
	low=left;
	high=right+1;
	median=(int)((left+right)*0.5);
	SWAP(list[left], list[median], temp);
	sort[index[left]]=median; sort[index[median]]=left;
	SWAP(index[left], index[median], temp);
	pivot=list[left];
	do{
		do low++;
		while(list[low]<pivot);
		do high--;
		while(list[high]>pivot);
		if(low<high){
			SWAP(list[low], list[high], temp);
			sort[index[low]]=high; sort[index[high]]=low;
			SWAP(index[low], index[high], temp);
		}
	}
	while(low<high);
	SWAP(list[left],list[high],temp);
	sort[index[left]]=high; sort[index[high]]=left;
	SWAP(index[left],index[high],temp);
	return high;
}

void rev_QuickSort(int list[], int index[], int sort[], int left, int right){
	if(left<right){
		int q = rev_partition(list, index, sort, left, right);
		rev_QuickSort(list, index, sort, left, q-1);
		rev_QuickSort(list, index, sort, q+1, right);
	}
}

int rev_partition(int list[], int index[], int sort[], int left, int right){
	int pivot, temp;
	int low, high, median;
	low=left;
	high=right+1;
	median=(int)((left+right)*0.5);
	SWAP(list[left], list[median], temp);
	sort[index[left]]=median; sort[index[median]]=left;
	SWAP(index[left], index[median], temp);
	pivot=list[left];
	do{
		do low++;
		while(list[low]>pivot);
		do high--;
		while(list[high]<pivot);
		if(low<high){
			SWAP(list[low], list[high], temp);
			sort[index[low]]=high; sort[index[high]]=low;
			SWAP(index[low], index[high], temp);
		}
	}
	while(low<high);
	SWAP(list[left],list[high],temp);
	sort[index[left]]=high; sort[index[high]]=left;
	SWAP(index[left],index[high],temp);
	return high;
}
