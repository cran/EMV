#include <R.h>
#include <Rmath.h>

double code_miss=-9999999;

double stdd(double *vector,int *length, int *is_finite);
double **dmatrix(int nb_row,int nb_col);
void mat_vec(double *array_vec,int* nb_row,int *nb_col,double **array);
void vec_mat(double *array_vec,int* nb_row,int *nb_col,double **array);
double *dvector(int length, int init);
int *ivector(int length, int init);
double  mean_vec(double *vector,int *length);
void free_dmatrix(double **array, int nb_row);
void init_ivector(int *vector, int *length, int  value);
void quicksort2(double *a, double *b, int *p, int *r);
int partition2(double *a, double *b, int p, int r);
int rand_part2(double *a, double*b,  int p, int r);
int uni_rand(int min,int max);
void init_dvector(double *vector, int *length, int value);
double  sum_vec(double *vector,int *length);


/****************************************************************************************************/
/*                                                                   stdd                                                                                                  */
/*  Purpose:   Return  the standard deviation of a vector                                                                                       */
/* Argument description:                                                                                                                                       */
/* vector: The sample vector                                                                                                                                  */
/* length: the length of the vector                                                                                                                          */
/* is_finite : (output) The number of finite values in vector                                                                                   */
/****************************************************************************************************/ 


double stdd(double *vector,int *length, int *is_finite)
{
  /* Compute the standard deviation of a vector */

  int i,count=0;
  double sum=0;
  double x_bar;
  double result;
  
  
  x_bar=mean_vec(vector,length);
  if(x_bar==code_miss)
    return(code_miss);
  else
    {
      for(i=0;i<*length;i++)
	{
	  if(vector[i]!=code_miss)
	    {
	      count=count+1;
	      sum=sum+(vector[i]-x_bar)*(vector[i]-x_bar);
	    }
	}
      *is_finite=count;
      if(count>1)
	{
	  result=(sqrt(sum/((double)count-1)));
	  return(result);
	}
      else
	{
	  return(code_miss);
	}
    }
}

/****************************************************************************************************/
/*                                                                   vec_mat                                                                                            */
/*  Purpose:   Coerce a vector into a matrix                                                                                                           */
/* Argument description:                                                                                                                                       */
/* array_vec: The vector to coerce                                                                                                                          */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/* array :(output) The two dimmensional array                                                                                                     */
/****************************************************************************************************/ 
void vec_mat(double *array_vec,int* nb_row,int *nb_col,double **array)
{
  int i,j;

  for(i=0;i<*nb_row;i++)
    for(j=0;j<*nb_col;j++)
      array[i][j]=array_vec[i**nb_col+j];
}

/****************************************************************************************************/
/*                                                                   mat_vec                                                                                            */
/*  Purpose:   Coerce a matrix into a vector                                                                                                           */
/* Argument description:                                                                                                                                       */
/* array_vec: (outpout)The vector                                                                                                                         */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/* array : The two dimmensional array to coerce                                                                                                    */
/****************************************************************************************************/ 

void mat_vec(double *array_vec,int* nb_row,int *nb_col,double **array)
{
  int i,j;

  for(i=0;i<*nb_row;i++)
    for(j=0;j<*nb_col;j++)
      array_vec[i**nb_col+j]=array[i][j];
}

/****************************************************************************************************/
/*                                                                   init_ivector                                                                                       */
/*  Purpose:   Initialize a vector of type int to zero                                                                                                */
/* Argument description:                                                                                                                                       */
/* vector: The vector to initialize                                                                                                                            */
/* length: the length of the vector                                                                                                                          */
/****************************************************************************************************/ 

void init_ivector(int *vector,int *length, int  value)
{
  int i;

  for(i=0;i<*length;i++)
    vector[i]=value;
}

/****************************************************************************************************/
/*                                                                   dmatrix                                                                                            */
/*  Purpose:  Allocate the memory for a matrix of type double                                                                             */
/* Argument description:                                                                                                                                       */
/* nb_row: The number of row for the matrix                                                                                                       */
/* nb_col: The number of column for the matrix                                                                                                    */
/****************************************************************************************************/ 
double **dmatrix(int nb_row,int nb_col)
{
  double **array;
  int i,j;

  /* Allocate the memory */
  array=Calloc(nb_row,double*);

  for(i=0;i<nb_row;i++)
    array[i]=Calloc(nb_col, double);

  /* Initialize to zero*/
  for(i=0;i<nb_row;i++)
    for(j=0;j<nb_col;j++)
      array[i][j]=0;

  return(array);
}
/****************************************************************************************************/
/*                                                                   dvector                                                                                             */
/*  Purpose:  Allocate the memory for a vector of type double                                                                              */
/* Argument description:                                                                                                                                       */
/* length: The length of the vector                                                                                                                         */
/* init: The value to initialize the vector                                                                                                                 */
/****************************************************************************************************/ 

double *dvector(int length, int init)
{
  int i;
  double *vector;

  /* Allocate the memory */
    vector=Calloc(length, double);
  
  /* Initialize the memory */
  for(i=0;i<length;i++)
    vector[i]=init;

  return(vector);
}

/****************************************************************************************************/
/*                                                                    ivector                                                                                             */
/*  Purpose:  Allocate the memory for a vector of type integer                                                                              */
/* Argument description:                                                                                                                                       */
/* length: The length of the vector                                                                                                                         */
/* init: The value to initialize the vector                                                                                                                 */
/****************************************************************************************************/ 


int *ivector(int length, int init)
{
  int i;
  int *vector;

  /* Allocate the memory */
  vector=Calloc(length, int);

  /* Initialize the memory */
  for(i=0;i<length;i++)
    vector[i]=init;

  return(vector);
}

/****************************************************************************************************/
/*                                                                mean_vec                                                                                             */
/*  Purpose:  Return the mean of vector (remove the missing values)                                                                   */
/* Argument description:                                                                                                                                       */
/* vector : The sample vector                                                                                                                                 */
/* length: The length of the vector                                                                                                                         */
/****************************************************************************************************/ 


double  mean_vec(double *vector,int *length)
{
  int i,count=0;
  double sum=0;

  for(i=0;i<*length;i++)
    {
      if(vector[i]!=code_miss)
	{
	  count=count+1;
	  sum=sum+vector[i];
	}
    }
  if (count>0)
    {
      return(sum/(double)count);
    }
  else
    {
      return(code_miss);
    }
  
}

/****************************************************************************************************/
/*                                                                free_dmatrix                                                                                        */
/*  Purpose:  Free the memory of a matrix of type double                                                                                     */
/* Argument description:                                                                                                                                       */
/* array: the two dimmensional array to free                                                                                                         */
/* nb_row : its number of row                                                                                                                               */
/****************************************************************************************************/ 


void free_dmatrix(double **array, int nb_row)
{

  int i;
  for(i=0;i<nb_row;i++)
    Free(array[i]);

  Free(array);
}
/****************************************************************************************************/
/*                                                                   quicksort2                                                                                        */
/*  Purpose:   Sort a vector using the quicksort algorithm  and move another vector at the same time                */
/* Argument description:                                                                                                                                       */
/* a: The vector to sort  from p to r                                                                                                                       */
/* b: The second vector to sort  from p to r                                                                                                            */
/* p: The first index                                                                                                                                                */
/* r: The last index                                                                                                                                                  */
/*****************************************************************************************************/ 


void quicksort2(double *a, double *b, int *p, int *r)
{
  int q;
  int q_p;

  if (*p<*r)
    {
    q=rand_part2(a,b, *p, *r);
    quicksort2(a,b,p,&q);
    q_p=q+1;
    quicksort2(a,b,&q_p,r);
    }
}
int partition2(double *a, double *b, int p, int r)
{
  double x=a[p];
  int i=p-1;
  int j=r+1;
  double temp;

  for(;;)
    {
      do
	{
	  j--;
	}while(a[j]>x);
      do
	{
	  i++;
	}while(a[i]<x);
      if(i<j)
	{
	  temp=a[i];
	  a[i]=a[j];
	  a[j]=temp;
	  temp=b[i];
	  b[i]=b[j];
	  b[j]=temp;
	}
      else
	{
	  return(j);
	}
    }

}

int rand_part2(double *a, double*b,  int p, int r)
{
  int i;
  double temp;

  i=uni_rand(p,r);
  temp=a[p];
  a[p]=a[i];
  a[i]=temp;
  temp=b[p];
  b[p]=b[i];
  b[i]=temp;

  return(partition2(a,b,p,r));

}

/****************************************************************************************************/
/*                                                                  uni_rand                                                                                            */
/*  Purpose:  Return a random number between min and max                                                                               */
/****************************************************************************************************/ 

int uni_rand(int min,int max)
{
  int rand_nb;


  GetRNGstate();
  rand_nb=(int)(unif_rand()*(max-min)+min);
  PutRNGstate();
  
  return(rand_nb);

}


/****************************************************************************************************/
/*                                                              init_dvector                                                                                           */
/*  Purpose:  Initialize a vector of type double                                                                                                       */
/* Argument description:                                                                                                                                       */
/* vector : The vector                                                                                                                                             */
/* length: The length of the vector                                                                                                                         */
/* value: The value to initialize it with                                                                                                                    */
/****************************************************************************************************/ 

void init_dvector(double *vector, int *length, int value)
{
  int i;

 
  /* Initialize to value */
  for(i=0;i<*length;i++)
    vector[i]=value;
}

/****************************************************************************************************/
/*                                                                sum _vec                                                                                             */
/*  Purpose:  Return the sum  of vector (remove the missing values)                                                                   */
/* Argument description:                                                                                                                                       */
/* vector : The sample vector                                                                                                                                 */
/* length: The length of the vector                                                                                                                         */
/****************************************************************************************************/ 

double  sum_vec(double *vector,int *length)
{
  int i,count=0;
  double sum=0;

  for(i=0;i<*length;i++)
    {
      if(vector[i]!=code_miss)
	{
	  count=count+1;
	  sum=sum+vector[i];
	}
    }
  if (count>0)
    {
      return(sum);
    }
  else
    {
      return(code_miss);
    }
}
