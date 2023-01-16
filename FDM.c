#include<stdio.h>
#include<stdlib.h>
void printmatrix(int row, int col, double *a){
    int i,j;
    printf("\nHere is your matrix:\n");
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
		    printf("%lf\t",*(a+i*row+j));
		}
	printf("\n");
	}
}

void reset_matrix(int row, int col, double *a){
	int i,j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
		    *(a+i*row+j)=0;
		}
	}
}

int substitution(int k,int i,int j,int n, double *K, double *R, double *boundary){
	int new_n;
	new_n=(n-2)*(n-2);
	*(K+k+(i-1)+3*(j-1))=-4;
	
	// 1. coner
	// coner1: (1,3)
	if( (i==1) && (j==(n-2)) )
	{		
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(R+k)-=(*boundary)*2;
		return 0;
	}
	// coner2: (1,1)
	if ((i==1) && (j==1)){
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1
		*(R+k)-=(*boundary)*2;
		return 0;
	}
	// coner3: (3,1)
	if ((i==(n-2)) && (j==1)){
		
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1
		
		*(R+k)-=(*boundary)*2;
		return 0;
	}
	// coner4: (3,3)
	if ((i==(n-2)) && (j==(n-2))){
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(R+k)-=(*boundary)*2;
		return 0;
	}
	
	// 2. line 
	if(i==1){
		//*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(R+k)-=(*boundary);
		return 0;
	}
	
	if(i==(n-2)){
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1
		//*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(R+k)-=(*boundary);
		return 0;
	}
	
	if(j==1){
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		//*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(R+k)-=(*boundary);
		return 0;
	}
	
	if(j==(n-2)){
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		//*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1	
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		*(R+k)-=(*boundary);
		return 0;
	}
	
	//3. Middle Point
	else{
		*(K+k+(i-1-1)+3*(j-1))=1; // (i-1,j) 1
		*(K+k+(i-1)+3*(j-1-1))=1; // (i,j-1) 1
		*(K+k+(i-1)+3*(j+1-1))=1; // (i,j+1) 1	
		*(K+k+(i+1-1)+3*(j-1))=1; // (i+1,j) 1
		//*(R+k)-=(*boundary);
		return 0;
	}
	
}

void SOR_method(int n, double *K,double *R,double *U){
	printf("\n");
	printf("================================================ SOR Method ================================================\n");
	printf("\n");
	int i,j,k;
	int flag=1.0;
	int new_n;
	int iter=100;
	double tolerance=0.0001;
	double diff=0.0;
	double sum=0.0;
	double omega=1.25;
	double *y;
	
	y=(double *)calloc(n,sizeof(double));
	new_n=(n-2)*(n-2);
	k=0;
   
   	do{ 
		for(i=0; i<new_n; i++){
		    flag=1;
		    sum=0.0;
		    for(j=0; j<i; j++) sum-=(*(K+i*new_n+j)*(*(y+j))); 
		    for(j=i+1; j<new_n; j++) sum-=(*(K+i*new_n+j)*(*(U+j)));
		    *(y+i)= (1.0-omega) * (*(U+i))+omega*(sum + *(R+i))/ *(K+i*new_n+i);
		    diff= *(y+i)-*(U+i);
		    if( diff<=0) diff*=-1.0;
		    if( diff>=tolerance) flag*=0;
		}

		if(flag ==1) {
		    printf("\n The Answer : x=");
		    for (j=0;j<new_n;j++) printf("%10.6f",*(y+j));
		    printf("\n");
		    printf("The error : e=:");
		    for (j=0; j<new_n;j++) printf("%10.6f", *(y+j)-*(U+j));
		    printf("\n");
		    break;
		}
	
		else{
		    printf("%3d: flag=%2d", k, flag);
		     for (j=0;j<new_n;j++) printf(" %10.6f", *(y+j));
		     printf("\n");
		     for (j=0;j<new_n;j++) *(U+j) = *(y+j);
		     k++;
		  }
		
	}while(flag<3 && k<iter);
}


int main(){
	int i,j,k=0;
	int n, new_n;
	double *K,*R,*U;
	double *a,*c;
	double *boundary;
	
	printf("Enter number of rows(columns): ");
	scanf("%d", &n);
	new_n=(n-2)*(n-2);
	
	a=(double *)calloc(n*n,sizeof(double));
	K=(double *)calloc(new_n*new_n, sizeof(double));
	R=(double *)calloc(new_n,sizeof(double));
	
	U=(double *)calloc(new_n,sizeof(double));
	
	c=(double *)calloc(new_n,sizeof(double));
	
	
	boundary=(double *)calloc(1,sizeof(double));
	*boundary=100;
	
// 1. Reset Matrix
	reset_matrix(new_n, new_n, K);
	reset_matrix(new_n, 1, R);
	reset_matrix(new_n, 1, U);
	reset_matrix(new_n, 1, c);
	//printmatrix(new_n,new_n,K);
	//printf("\n");
	
// 2. K,R Matrix setting
	for(i=1;i<n-1;i++){
		for(j=1;j<n-1;j++){
			substitution(k,j,i,n,K,R,boundary);
			k+=new_n;
		}
	}
	//printmatrix(new_n,new_n,K);
	//printmatrix(new_n, 1, R);
	printf("\n");
	printf("============= K Matrix ============= \n");
	printmatrix(new_n,new_n,K);
	printf("\n");
	printf("============= R Matrix ============= \n");
	printmatrix(new_n, 1, R);
// 3. Adding SOR Method
	SOR_method(n,K,R,U);
	
// 4. Print Temperature distribution
	printf("\n");
	printf("============= Here is your temperature distribution matrix: ============= \n");	
	printf("\n");
	for (i=0;i<n-2;i++){
		printf("\t\t");
		for(j=0;j<n-2;j++){
			printf("%lf\t",*(U+i*(n-2)+j));
		}
		printf("\n");
	}
	
	return 0;
}