////////////////////////////////////////////////////////////////////////////////////
// File       : MIS_sequential.c
// Language   : C/C++ and MPI
// Description: Solving Maximum Independent Set (MIS) problem sequentially
// Author     : Sandesh Prabhakar
// Created    : November 9, 2008
// Modified   : August 19, 2026
// Version    : 1.0
//
// Copyright (c) 2026 Sandesh Prabhakar
// All rights reserved
//////////////////////////////////////////////////////////////////////////////////////

#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

/* Local functions */

static void master(void);
static int MIS_Solver(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c, int **B, int *index);
static int MIS_Solver_worker(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c);

int main(int argc, char **argv)
{
  int myrank;
  double endTime, startTime, time;
  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    	/* start timer */
	startTime = MPI_Wtime();
        master();
        /* stop timer */
	endTime = MPI_Wtime();
	time = endTime - startTime;
       	fprintf(stdout,"\n Total Sequential time (Ts) = %e seconds for n = %d \n",time, 50);fflush(stdout); 
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


static void master(void)
{

        double temp, theta, p;
	int **A,**B, **L, *data_packet, tmp, i, j, n, k, m, z, c, v, *X, *degree, *degree_sort, *vertices_sort, mod_S;
        int *S, *I, subtree_count; 
        int temp_best=0, best=0, depth=1, y=5, y_temp=5;
        int *index;
        int temp_sort;
	n=50, theta = 0.5;
        k=pow(2,y);
        mod_S=n;

        int ntasks, rank;
        MPI_Status status;
        A = new int*[n+1]; 
        for(i=0;i<=n;i++)
	{
        A[i] = new int[n+1];
	}
        X = new int [n+1];
        degree = new int [n+1];
        degree_sort = new int[n+1]; 
        vertices_sort = new int[n+1];
        S = new int [n+1];
        I = new int[n+1];
        data_packet = new int[n+4]; 
        index = new int[1];
        
        B = new int*[k+1];
        L = new int*[k+1];
       for(i=0;i<=k;i++)
	{
          B[i] = new int[n+2];
          L[i] = new int[n+2];
	}
        index[0]=1;
        for(i=1;i<=k;i++)
	{
             for(j=0;j<=n+1;j++)
                {
                  L[i][j]=0;
                  B[i][j]=0;
		}
	}

        for(i=1;i<=n;i++)
	{
        degree[i]=0;
        degree_sort[i]=0;
        vertices_sort[i]=0;
        S[i]=i;
        I[i]=0; 
	}
        fprintf(stdout,"\n Initializing graph.\n");fflush(stdout);
        for(i=1; i<=n; i++)
        {
	  for(j=i; j<=n; j++)
                {
                   //temp = rand()+0.1;
                   //p=(rand()+0.1)/(2*temp);
	           srand48(i);	   
                   temp = drand48()+0.1;
                   p=(drand48()+0.1)/(2*temp);
                   if(i==j)
                   {
	           A[i][j]=0; 
                   }
                   else if(p<theta)
                   {
	           A[i][j]=1; 
                   A[j][i]=1;
                   }
                   else
                   {
                   A[i][j]=0;
                   A[j][i]=0;
                   }
                }
        }
        fprintf(stdout,"\n Random Graph \n");fflush(stdout);
        fprintf(stdout,"\n ------------ \n");fflush(stdout);
        for(i=1;i<=n;i++)
           for(j=1;j<=n;j++)
               {
               fprintf(stdout,"\n A[%d][%d] = %d \n", i, j, A[i][j]);fflush(stdout);
               }

        fprintf(stdout,"\n");fflush(stdout);
      for(i=1;i<=n;i++) 
      for(j=1;j<=n;j++)
      {
      degree[i] = degree[i] + A[i][j];
      degree_sort[i] = degree_sort[i] + A[i][j];
      }
    for ( i = 1; i < n; i++)
        {
         for(j = 1; j <= n-i; j++)
          {

             if (degree_sort[j] > degree_sort[j+1])
             {
                    temp_sort = degree_sort[j];
                    degree_sort[j] = degree_sort[j+1];
                    degree_sort[j+1] = temp_sort;
              }
            }
         }
           c=1;
           for(i=1;i<=n;i++)
            {
             for(j=1;j<=n;j++)
                {
                   if(i>1) 
                   {
                   if(degree_sort[i] == degree_sort[i-1]) 
                   break;
		   }
                if(degree_sort[i] == degree[j])
                   {
                     vertices_sort[c] = j;
                     c++;
                   }
                }
            }
c=1;
while(mod_S!=0)
{
  v=vertices_sort[c]; 

  if(A[v][v]==1)
  {
  c--;
  break; 
  }


  I[c]=v;
  
  if(S[v]!=0)
   {
    S[v]=0;
    mod_S--;
   }
  for(i=1;i<=n;i++) 
  {
   
    if(A[v][i]==1)
    {
       if(S[i]!=0)
        {
        S[i]=0;
        mod_S--;
        }
    }
  }
int flag=0;
  for(j=1;j<c;j++)
    {
     if(A[v][I[j]]==1)
     {
      mod_S=0; 
      I[c]=0;
      c--;
      flag=1;
      break;
     }
    }

if(mod_S!=0 && flag==0)
c++;

}

fprintf(stdout,"\n\n Maximal Independent Set \n");fflush(stdout);
fprintf(stdout,"\n ----------------------- \n");fflush(stdout);
fprintf(stdout,"\n |Maximal IS| = %d \n",c);fflush(stdout);
fprintf(stdout,"\n");fflush(stdout);
fprintf(stdout," [");fflush(stdout);
for(i=1;i<=c;i++)
{
fprintf(stdout," %d ",I[i]);fflush(stdout);
}
fprintf(stdout," ]");fflush(stdout);
fprintf(stdout,"\n\n");fflush(stdout);

best = MIS_Solver(X, depth, best, temp_best, y, n, A, c, B, index);

subtree_count=0;
i=1;
while(B[i][0]!=0)
{
subtree_count++;
i++;
if(i>k)
break;
}
fprintf(stdout,"\n\n Subtree Count = %d \n",subtree_count);fflush(stdout);

fprintf(stdout,"\n\n Independent Sets \n");fflush(stdout);
fprintf(stdout,"\n ---------------- \n");fflush(stdout);
for(i=1;i<=subtree_count;i++)
for(j=0;j<=y_temp;j++)
{
L[i][j]=B[i][j];
B[i][j]=0;
}
y=n;

int *MIS, *MIS_sort;
MIS = new int[subtree_count+1];
MIS_sort = new int[subtree_count+1];

for (j=1;j<=subtree_count;j++)
{
index[0]=1;
depth=y_temp+1;  
temp_best=L[j][0]-1;
for(i=1;i<=y_temp;i++)
X[i]=L[j][i];
MIS[j] = MIS_Solver_worker(X, depth, best, temp_best, y, n, A, c);
MIS_sort[j]=MIS[j];
}


    for ( i = 1; i < subtree_count; i++)
        {
         for(j = 1; j <= subtree_count-i; j++)
          {

             if (MIS_sort[j] > MIS_sort[j+1])
             {
                    temp_sort = MIS_sort[j];
                    MIS_sort[j] = MIS_sort[j+1];
                    MIS_sort[j+1] = temp_sort;
              }
            }
         }

printf("\n\n Maximum Independent Set \n");
printf("\n ----------------------- \n");
printf("\n |MIS| = %d \n\n", MIS_sort[subtree_count]);

fprintf(stdout,"\n MASTER EXITS \n");fflush(stdout);
	return;

}


int MIS_Solver(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c, int **B, int *index)
    {
              bool flag=0;
              int i; 

              X[depth]=0;

              if(depth==n)
              {
              int l;
              for(l=1;l<=n;l++)
              B[index[0]][l]=X[l];
              B[index[0]][0]=temp_best+1;
              index[0]++;
              }

              if(depth<n && (best+(y-depth))>=c)
	      {
              best=MIS_Solver(X, depth+1, best, temp_best, n, y, A, c, B, index);
	      }

              X[depth]=1;

	      //if(A[depth][depth]==1)
              //flag=1;

              for(i=depth-1; i>0; i--) 
              {
               if(X[i]==1)
               {
                   if(A[i][depth]==1)
                   flag=1;
	       }
	      }

              if(flag==0)
	      { 
              if(temp_best>=best) 
                {
                best=temp_best+1;
		}

              if(depth<n && (best+(y-depth))>=c)
              best=MIS_Solver(X, depth+1, best, temp_best+1, n, y, A, c, B, index);

              if(depth==n)
              {
              int l;
              for(l=1;l<=n;l++)
              B[index[0]][l]=X[l];
              B[index[0]][0]=temp_best+2;
              index[0]++;
              }
	      } 
             
              else return best;
return best;

    }


int MIS_Solver_worker(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c)
    {
              bool flag=0;
              int i; 

              X[depth]=0;
             if(depth<n && (best+(y-depth))>=c)
	      {
              best=MIS_Solver_worker(X, depth+1, best, temp_best, n, y, A, c);
	      }

              X[depth]=1;
	      //if(A[depth][depth]==1)
              //flag=1;
              for(i=depth-1; i>0; i--) 
              {
               if(X[i]==1)
               {
                   if(A[i][depth]==1)
                   flag=1;
	       }
	      }
              if(flag==0)
	      { 
              if(temp_best>=best) 
                {
                best=temp_best+1;
                if(best>=c)
                    {
                     printf("\n [");
                     for(i=1;i<=depth;i++)
                       {
                       if(X[i]==1)
                        {
                        fprintf(stdout," %d ",i);fflush(stdout);
                         }
                       }
                      printf("] ");
                      fprintf(stdout,"\n \n");fflush(stdout);
		    }
		}
              if(depth<n && (best+(y-depth))>=c)
	      {
              best=MIS_Solver_worker(X, depth+1, best, temp_best+1, n, y, A, c);
	      }
	      } 
              
              else return best;
        return best;
    }

