#include "util.h"



  
int is_na(double *array,int* nb_col,int *array_nb); 
void fill_up(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos, double *temp, double *dist_bound);
double distance(double *array1,double *array2, int *length);
double correlation(double *array1,double *array2, int *length);
int comp_na(double *vector, int *length, int *position);
void neighboors(double **array, int *nb_row, int *nb_col, int *n_position, int *nb_neighboors);
void fill_up_corr(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos, double *temp, double *dist_bound);

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



void knnc(double *array_vec,int *nb_col,int *nb_row, int*k, int *corre_flag, double *dist, double *dist_bound)
{
  int missing,i,j,ii;
  int count;
  double value;
  double *temp;
  double *row_nb;
  double ** array;
  int *miss_pos;
  int index;
  int *n_position;
  int* nb_neighboors;

  int min=0;
  int max=*k-1;
  
  array=dmatrix(*nb_row,*nb_col);

  /** contain the row numbers of the missing values **/
  miss_pos=ivector(*nb_col, code_miss);

  /** contains the distances of the neighboors **/
  temp=dvector(*k,code_miss); 
  /** contains the row numbers of the neighboors **/
  row_nb=dvector(*k,code_miss);
  /** initilize all the distances with the missing codes **/
  init_dvector(dist, nb_row, code_miss);

  n_position=ivector(*nb_row, code_miss); /** positions of potential neighboors **/
  nb_neighboors=ivector(1, code_miss); /** number of neighboors **/

  /** coerce the vector into a two dimmensional array **/
  vec_mat(array_vec,nb_row,nb_col,array); 
 
  neighboors(array, nb_row, nb_col, n_position, nb_neighboors);

  if(*nb_neighboors==0) /** Stop if no neighboors **/ 
    {
      error("No rows without missing values"); 
    }
  else 
    {
      if(*nb_neighboors<*k) /** If less than k neighboors give a warning **/
	warning("Only %d neighboors could be used", *nb_neighboors); 
      
      for(i=0;i<*nb_row;i++)
	{
	  /** Check for missing values **/
	  missing=is_na(array[i],nb_col,miss_pos);
	  
	  if (missing==1 && miss_pos[*nb_col-1]==code_miss) /**at least one missing value at most nb_col**/
	    {
	      if(*corre_flag==1 && miss_pos[*nb_col-2]!=code_miss) /** Give a warning if based on correlation and only one observation **/
		warning("Could not estimate the missing values for the row %d\n One observation is not enough to compute the sample correlation", i+1); 
	      else
		{
		  count=0;  
		  
		  for(j=0;j<*nb_neighboors;j++)  /** loop on the neighboors only **/
		    { 
		      index=n_position[j];
		      
		      if(*corre_flag==0)
			value=distance(array[i],array[index],nb_col);  /** compute the distance **/
		      else
			value=-correlation(array[i],array[index],nb_col);  /** compute the correlation **/
		      
		      if(value!=code_miss)
			{
			  if (count<*k) /** store the first k **/
			    {
			      temp[count]=value;
			      row_nb[count]=index;
			      count++;
			    }
			  else
			    {
			      quicksort2(temp,row_nb,&min,&max); /** sort the neighboors to keep the kth nearest **/
			      if (temp[*k-1]>value)  /** keep it if the distance is shorter **/
				{
				  temp[*k-1]=value;
				  row_nb[*k-1]=index;
				}       
			    } 
			}
		      
		    } 
		  
		  if(*corre_flag==0)
		    {
		      fill_up(array,row_nb,nb_col,k,i,miss_pos,temp, dist_bound); /** fill up the missing values by the averaging the distance**/
		      dist[i]=mean_vec(temp, k); /** Compute the average distances **/
		    }
		  else
		    {
		      fill_up_corr(array,row_nb, nb_col, k,i, miss_pos, temp, dist_bound); /** fill up the missing values based on correlations**/
		      dist[i]=-mean_vec(temp, k); /** Compute the average distances **/
		    }
		  
		  
		  init_dvector(row_nb, k, code_miss);    /** initialize row_nb with missing codes **/
		  init_dvector(temp, k, code_miss);        /** initialize temp with missing codes **/
		}
	    }
	  else if(missing==1 && miss_pos[*nb_col-1]!=code_miss)
	    warning("Could not estimate the missing values for the row %d\n The row only contains missing values", i+1); 
	}
    }
  
  mat_vec(array_vec, nb_row, nb_col,array); /** recoerce the matrix into a vector **/

  /** free the memory **/
  free_dmatrix(array,*nb_row); 
  Free(miss_pos);
  Free(temp);
  Free(row_nb);
  Free(n_position);
  Free(nb_neighboors);
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
/* temp : the distances of the neighboors                                                                                                               */
/* dist_bound : the bound for the distances                                                                                                           */
/****************************************************************************************************/


void fill_up(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos, double *temp, double *dist_bound)
{
  int i,index,j,count;
  double sum;
  int missing_spot;
  int flag=0;
  
  i=0;
  while(miss_pos[i]!=code_miss) /* The positions are stored at the beginning then zeros */
    {
      
      missing_spot=miss_pos[i];
      count=0;
      sum=0;
      for (j=0;j<*k;j++)
	{
	  index=row_nb[j];
	  if(index!=code_miss)
	    {
	      if(temp[j]<*dist_bound | *dist_bound==0)
		{
		  sum=sum+array[index][missing_spot];
		  count=count+1;
		}
	      else
		{
		  row_nb[j]=code_miss;
		  temp[j]=code_miss;
		}
	    }
	}
      if(count>0)
	{
	  array[position][missing_spot]=sum/count; /** replace the missing value by the average of the neighboors **/
	  flag=1; 
	}
      
      if(flag==0)
	 warning("Could not estimate the missing values for the row %d\n  dist.bound is too small", position+1); 

      i++;
    }
  
}


/***************************************************************************************************/
/*                                                    fill_up based on the correlations                                                                   */
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


void fill_up_corr(double **array,double *row_nb,int *nb_col,int *k,int position, int* miss_pos, double *temp, double *dist_bound)
{
  int i,index,j,count,kk;
  double sum;
  int missing_spot;
  int *finite;
  double sdx,meanx;
  double *tmp,*meany,*sdy;
  int flag=0;
  
  tmp=dvector(*nb_col,code_miss);
  meany=dvector(*k,code_miss);
  sdy=dvector(*k,code_miss);

  finite=ivector(1,code_miss);
  sdx=stdd(array[position],nb_col, finite);
  meanx=mean_vec(array[position],nb_col);

  for(kk=0;kk<*k;kk++)
    {
      index=row_nb[kk];
      if(index!=code_miss)
	{
	  for(j=0;j<*nb_col;j++)
	    {
	      if(array[position][j]!=code_miss)
		tmp[j]=array[index][j];
	      else
		tmp[j]=code_miss;
	    }
	  meany[kk]=mean_vec(tmp,nb_col);
	  sdy[kk]=stdd(tmp,nb_col,finite);
	}
    }

  i=0;
  while(miss_pos[i]!=code_miss ) /** The positions are stored at the beginning then code_miss **/ 
    {
      
      missing_spot=miss_pos[i];
      count=0;
      sum=0;
      for (j=0;j<*k;j++) /** Compute the average of the z-scores **/ 
	{
	  index=row_nb[j];
	  
    	  if(index!=code_miss)/** if not enough neighboors, do not use the missing ones **/ 
	    { 
	      if(temp[j]<-(*dist_bound) |*dist_bound==0)
		{
		  sum=sum+(array[index][missing_spot]-meany[j])/sdy[j];
		  count=count+1;
		}
	      else
		{
		  row_nb[j]=code_miss;
		  temp[j]=code_miss;
		}
	    }
	}
      if(count>0)
	{
	  array[position][missing_spot]=meanx+sdx*(sum/count); /** replace the missing value by the average of the standardized neighboors **/ 
	  flag=1;
	}

      i++;
    }
  if(flag==0)
    warning("Could not estimate the missing values for the row %d\n  dist.bound is too large", position+1); 

  /** free the memory **/

  Free(finite); 
  Free(tmp); 
  Free(sdy); 
  Free(meany); 

}


/***************************************************************************************************/
/*                                                    is_na                                                                                                              */
/* Purpose: Check if a vector contains missing values                                                                                          */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array : the one dimmensional array                                                                                                                  */ 
/* nb_col: The length of the array                                                                                                                         */
/* miss_pos: (output)the positions of the missing values                                                                                      */
/****************************************************************************************************/

int is_na(double *array,int* nb_col,int *miss_pos)
{
  int i;
  int count=0;
  
  init_ivector(miss_pos,nb_col,code_miss);
  
  for(i=0;i<*nb_col;i++)
    {
      if(array[i]==code_miss)
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
 
      
      if(array1[i]!=code_miss) /** omit the missing values **/
	if(array2[i]!=code_miss)
	  {
	    count=count+1;
	    value=value+(array1[i]-array2[i])*(array1[i]-array2[i]);
	  }
    }
  if(count>0)
    {
  return(sqrt(value));
    }
  else
    {
      return(code_miss);
    }
}

/***************************************************************************************************/
/*                                                    correlation                                                                                                     */
/* Purpose: Compute the correlation between 2 vectors omitting the missing values                                          */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array1 : the first one dimmensional array                                                                                                         */ 
/* array2 : the second  dimmensional array                                                                                                           */
/* length: the length of the arrays                                                                                                                         */
/****************************************************************************************************/
double correlation(double *array1,double *array2, int *length)
{
  int i;
  int count=0;
  double value=0;
  double meanx,meany;
  double sdx,sdy;
  int *finite1,*finite2;
  double *temp;
  
  finite1=ivector(1,code_miss);
  finite2=ivector(1, code_miss);
  temp=dvector(*length,code_miss);

  for(i=0;i<*length;i++)
    {
      if(array1[i]!=code_miss)
	temp[i]=array2[i];
    }
  

  meanx=mean_vec(array1, length);
  meany=mean_vec(temp, length);

  sdx=stdd(array1,length,finite1);
  sdy=stdd(temp,length,finite2);
  

  if(*finite1 >1 && *finite2>1) /* check if more than two observations for both rows */
    {
      for(i=0;i<*length;i++)
	{
	  if(array1[i]!=code_miss) /** omit the missing values **/
	    if(array2[i]!=code_miss)
	      {
		count=count+1;
		value=value+(array1[i]-meanx)*(array2[i]-meany);
	      }
	}
      
      if(count>1)
	{
	  value=(value/(count-1))/(sdx*sdy);
	}
      else
	{
	  value=code_miss;
	}
    }
  else
    value=code_miss;


  /** Free the memory **/
  Free(temp); 
  Free(finite1);
  Free(finite2);


  return(value);
  
}


/***************************************************************************************************/
/*                                                    neighboors                                                                                                    */
/* Purpose: Look for the rows without missing values                                                                                        */
/*                                                                                                                                                                          */
/* Argument description                                                                                                                                       */
/* array : the data matrix with missing values                                                                                                      */ 
/* nb_row : the number of rows                                                                                                                           */
/* nb_col : the number of col                                                                                                                                */
/* n_position (output) : the indexes of the rows with no missing values                                                              */
/* nb_neighboors : the number of rows with no missing values  (i.e. the potential neighboors)                          */
/****************************************************************************************************/




void neighboors(double **array, int *nb_row, int *nb_col, int *n_position, int *nb_neighboors)
{
  int i;
  int *niet; /* ghost variable */
  int count=0;
  int missing;

  niet=ivector(*nb_col, code_miss);
  
  for(i=0;i<*nb_row;i++)
    {
      missing=is_na(array[i],nb_col,niet); /** check if the neighboor has missing values **/
      if(missing==0) 
	{
	  n_position[count]=i;
	  count++;
	}
    }
  *nb_neighboors=count;

  Free(niet);
}
