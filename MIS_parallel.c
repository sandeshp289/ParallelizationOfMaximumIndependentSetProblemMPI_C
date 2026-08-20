////////////////////////////////////////////////////////////////////////////////////
// File       : MIS_parallel.c
// Language   : C/C++ and MPI
// Description: Parallelization of Maximum Independent Set (MIS) problem
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

#define WORKTAG 1
#define DIETAG 2


/* Local functions */

static void master(void);
static void slave(int);
static int MIS_Solver(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c, int **B, int *index);
static int MIS_Solver_worker(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c, int myrank);

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
       	fprintf(stdout,"\n Total Parallel time = %e \n",time);fflush(stdout); 
  } else {
    slave(myrank);
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
	n=10, theta = 0.5;
        k=pow(2,y);
        mod_S=n;

        int ntasks, rank;
        MPI_Status status;
        A = new int * [n+1]; 
        for(i=0;i<=n;i++)
	{
        A[i] = new int[n+1];
	}
        X = new int [n+1];
        degree = new int[n+1];
        degree_sort = new int[n+1]; 
        vertices_sort = new int [n+1];
        S = new int [n+1];
        I = new int[n+1];
        data_packet = new int [n+4]; 
        index = new int[1];
        
        B = new int *[k+1];
        L = new int* [k+1];
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

int *MIS, *MIS_sort;
MIS = new int[subtree_count+1];
MIS_sort = new int[subtree_count+1];
  /* Find out how many processes there are in the default
     communicator */

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  /* Seed the slaves; send one unit of work to each slave. */
  fprintf(stdout,"\n Number of processors=%d \n",ntasks);fflush(stdout);
  for (rank = 1; rank < ntasks; ++rank) {
    /* Find the next item of work to do */
    temp_best=L[rank][0]-1;
    for(i=1;i<=n;i++)
    {
    X[i]=L[rank][i];
    data_packet[i]=X[i];
    }
    data_packet[n+1]=best;
    data_packet[n+2]=temp_best;
    data_packet[n+3]=c;
    /* Send it to each rank */

    MPI_Send(data_packet,       /* message buffer */
             n+4,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

  //  fprintf(stdout,"\n Master sends data_packet of subtree %d to slave %d \n",rank, rank );fflush(stdout);
//for(i=1;i<=n+3;i++)
//fprintf(stdout,"\n Master sends message data_packet[]=%d to slave %d\n", i,data_packet[i], rank );fflush(stdout);

  }

  /* Loop over getting new work requests until there is no more work
     to be done */
  tmp=ntasks;
  j=1;
  while(tmp<=subtree_count)
  {
    temp_best=L[tmp][0]-1;
    for(i=1;i<=n;i++)
    {
    X[i]=L[tmp][i];
    data_packet[i]=X[i];
    }    
    data_packet[n+1]=best;
    data_packet[n+2]=temp_best;
    data_packet[n+3]=c;


    /* Receive results from a slave */

    MPI_Recv(&MIS[j],           /* message buffer */
             1,                 /* one data item */
             MPI_INT,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

/* Check the tag of the received message. */
    if (status.MPI_TAG == 0) {
    
    MIS_sort[j]=MIS[j];
//fprintf(stdout,"\n Master receives message MIS=%d from slave %d \n", MIS[j], status.MPI_SOURCE);fflush(stdout);
    j++;
    }
    /* Send the slave a new work unit */
    MPI_Send(data_packet,             /* message buffer */
             n+4,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

//fprintf(stdout,"\n Master sends data_packet of subtree %d to slave %d \n",tmp, status.MPI_SOURCE );fflush(stdout);
//for(i=1;i<=n+3;i++)
//fprintf(stdout,"\n Master sends message data_packet[%d]=%d to slave %d \n", i,data_packet[i], status.MPI_SOURCE );fflush(stdout);
    /* Get the next unit of work to be done */
    tmp++;
  }
  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(&MIS[j], 1, MPI_INT, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);

 if (status.MPI_TAG == 0) {
    MIS_sort[j]=MIS[j];
 //   fprintf(stdout,"\n Master receives message MIS[%d]=%d from slave %d\n", j, MIS[j], status.MPI_SOURCE);fflush(stdout);
    j++;
 }
  }

//for(i=1;i<=subtree_count;i++)
//fprintf(stdout,"\n MIS[%d]=%d MIS_sort[%d]=%d \n", i, MIS[i], i, MIS_sort[i]);fflush(stdout);

//fprintf(stdout,"\n subtree_count_last=%d \n", subtree_count);fflush(stdout);

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

fprintf(stdout,"\n\n Maximum Independent Set \n");fflush(stdout);
fprintf(stdout,"\n ----------------------- \n");fflush(stdout);
fprintf(stdout,"\n |MIS| = %d \n\n", MIS_sort[subtree_count]);fflush(stdout);

  /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
//fprintf(stdout,"\n Master sends message 0 and DIETAG to slave %d \n", rank );fflush(stdout);
  }

fprintf(stdout,"\n MASTER EXITS \n");fflush(stdout);
return;
}


static void slave(int myrank)
{

        double temp, theta, p;
	int **A, *data_packet, i, j, n, k, m, z, c, v, *X;
        int temp_best=0, best=0, depth, y, y_temp=5;
	n=10, theta = 0.5;
        y=n;
        depth=y_temp+1;
        MPI_Status status;

        A = new int *[n+1]; 
        for(i=0;i<=n;i++)
	{
        A[i] = new int[n+1];
	}
        X = new int[n+1];
        data_packet = new int[n+4];
        //printf("\n Initializing graph.\n");
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

    while (1) {
    /* Receive a message from the master */

    MPI_Recv(data_packet, n+4, MPI_INT, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

//for(i=1;i<=n+3;i++)
//fprintf(stdout,"\n Slave %d receives message data_packet[%d]=%d from Master\n", myrank, i,data_packet[i] );fflush(stdout);


    /* Check the tag of the received message. */

    if (status.MPI_TAG == DIETAG) {
    fprintf(stdout,"\n SLAVE %d DIES \n", myrank);fflush(stdout);
      return;
    }
    else if(status.MPI_TAG == WORKTAG)
    {
    /* Do the work */
  //  fprintf(stdout,"\n Slave %d receives data_packet %d from Master\n", myrank, j );fflush(stdout);
    for(i=1;i<=n;i++)
    X[i]=data_packet[i];
    best = data_packet[n+1];
    temp_best = data_packet[n+2];
    c = data_packet[n+3];
    best = MIS_Solver_worker(X, depth, best, temp_best, y, n, A, c, myrank);  
    /* Send the result back */
    MPI_Send(&best, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
 //   fprintf(stdout,"\n Slave %d sends message MIS=%d \n to Master", myrank, best );fflush(stdout);
    }
  }
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


int MIS_Solver_worker(int *X, int depth, int best, int temp_best, int n, int y, int **A, int c, int myrank)
    {
              bool flag=0;
              int i; 

              X[depth]=0;
             if(depth<n && (best+(y-depth))>=c)
	      {
              best=MIS_Solver_worker(X, depth+1, best, temp_best, n, y, A, c, myrank);
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
              best=MIS_Solver_worker(X, depth+1, best, temp_best+1, n, y, A, c, myrank);
	      }
	      } 
              
              else return best;
        return best;
    }

