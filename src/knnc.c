#include "util.h"


int is_na(double *array,int* nb_col,int *array_nb); 
void fill_up(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos);
double distance(double *array1,double *array2, int *length);
int comp_na(double *vector, int *length, int *position);

/***************************************************************************************************/
/*                                                    knn                                                                                                                */
/* Purpose: Replace the missing values in a two dimmensional array using a k-th nearest neighboors               */
/* algorithm. The algorithm select a row with missing values, then search the k-th nearest rows using the       */
/* Euclidian distance.  Then the missing values are replaced by the mean of the k-th nearest neighboors           */
/* The code for missing values is -99.                                                                                                                   */
/*                                                                                                                                                                           */
/* Argument description                                                                                                                                        */
/* array_vec: A one array vector that will be coerce into a 2 dimensionnal array                                                  */ 
/* nb_col: The number of columns                                                                                                                         */
/* nb_row : The number of rows                                                                                                                           */
/* k: The number of neighboors                                                                                                                            */
/****************************************************************************************************/



void knnc(double *array_vec,int *nb_col,int *nb_row, int*k)
{
  int missing,i,j,ii;
  int comp_missing;
  int count;
  double value;
  double *temp;
  double *row_nb;
  double ** array;
  int *miss_pos;
  int *niet;

  int min=0;
  int max=*k-1;
  
  array=dmatrix(*nb_row,*nb_col);

/** contain the row numbers of the missing values **/
  miss_pos=ivector(*nb_col);

  /** niet is a ghost vector argument required by the is na funtion **/
  niet=ivector(*nb_col);
  temp=dvector(*k);
  row_nb=dvector(*k);

  /** coerce the vector into a two dimmensional array **/
  vec_mat(array_vec,nb_row,nb_col,array); 
 

  for(i=0;i<*nb_row;i++)
    {
      /** Check for missing values **/
      missing=is_na(array[i],nb_col,miss_pos);
      
      if (missing==1 && miss_pos[*nb_col-1]==0) /**at least one missing value at most nb_col**/
	{
	  count=0; 
	  for(j=0;j<*nb_row;j++) 
	     { 
	       missing=is_na(array[j],nb_col,niet); /** check if the neighboor has missing values **/
	       if(missing==0) 
		 {  
		   /** if no missing values can be use as a neighboor **/
		   value=distance(array[i],array[j],nb_col);  /** compute the distance **/
		   
		   if(value!=-9999999)
		     {
		       if (count<*k)
			 {
			   temp[count]=value;
			   row_nb[count]=j;
			   count++;
			 }
		       else
			 {
			   quicksort2(temp,row_nb,&min,&max); /** sort the neighboors to keep the kth nearest **/
			   if (temp[*k-1]>value)  /** keep it if the distance is shorter **/
			     {
			       temp[*k-1]=value;
			       row_nb[*k-1]=j;
			     }
			 }
		       
		    } 
		 }
	      
	     } 

	   fill_up(array,row_nb,nb_col,k,i,miss_pos); /** fill up the missing values **/
	}
    }
  
  mat_vec(array_vec, nb_row, nb_col,array); /** recoerce the matrix into a vector **/
}

/***************************************************************************************************/
/*                                                    fill_up                                                                                                            */
/* Purpose: Compute the average of the neighboors and replace the missing values                                           */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array : the two dimmensional array                                                                                                                  */ 
/* nb_col: The number of columns                                                                                                                        */
/* nb_row : The number of rows                                                                                                                          */
/* k: The number of neighboors                                                                                                                           */  
/* position : the positions of the neighboors                                                                                                         */
/* miss_pos : the positions of the missing values                                                                                                   */
/****************************************************************************************************/

  
void fill_up(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos)
{
  int i,index,j,count;
  double sum;
  int missing_spot;

  i=0;
  while(miss_pos[i]!=0 || i==0) /* The positions are stored at the beginning then zeros */
    {
      
      missing_spot=miss_pos[i];
      count=0;
      sum=0;
      for (j=0;j<*k;j++)
	{
	  index=row_nb[j];
	  if(array[index][missing_spot]!=-9999999)
	    {
	      sum=sum+array[index][missing_spot];
	      count=count+1;
	    }
	}
      array[position][missing_spot]=sum/count; /* replace the missing value by the average of the neighboors*/
      i++;
    }
  
}

/***************************************************************************************************/
/*                                                    is_na                                                                                                              */
/* Purpose: Check if a vector contains missing values                                                                                          */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array : the one dimmensional array                                                                                                                  */ 
/* nb_col: The length of the array                                                                                                                     */
/* miss_pos: (output)the positions of the missing values                                                                                      */
/****************************************************************************************************/

int is_na(double *array,int* nb_col,int *miss_pos)
{
  int i;
  int count=0;
  
  iintit_vect(miss_pos,*nb_col);

  for(i=0;i<*nb_col;i++)
    {
      if(array[i]==-9999999)
	{
	  miss_pos[count]=i;
	  count=count+1;
	}
    }
  if(count>0)
    return(1); /** one if missing's **/
  else
    return(0); /** 0 if no missing's **/
  
}
/***************************************************************************************************/
/*                                                    distance                                                                                                         */
/* Purpose: Compute the Euclidean distance between 2 vectors omitting the missing values                              */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array1 : the first one dimmensional array                                                                                                         */ 
/* array2 : the second  dimmensional array                                                                                                           */
/* length: the length of the arrays                                                                                                                         */
/****************************************************************************************************/
double distance(double *array1,double *array2, int *length)
{
  int i;
  int count=0;
  double value=0;
  
  for(i=0;i<*length;i++)
    {
 
      
      if(array1[i]!=-9999999) /** omit the missing values **/
	if(array2[i]!=-9999999)
	  {
	    count=count+1;
	    value=value+(array1[i]-array2[i])*(array1[i]-array2[i]);
	  }
    }
  if(count>0)
    {
  return(value);
    }
  else
    {
      return(-9999999);
    }
}
