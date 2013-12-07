#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "mpi.h"

#define CHAR char
#define REAL double
#define INT  int

#define OUTPUT stdout       /* output to standard out                     */
#define PLOT_FILE "plots"   /* output files base name                     */
#define INCREMENT 500       /* number of steps between convergence check  */

#define P 1                 /* define processor count for serial codes    */
#define K 0                 /* current thread number for serial code is 0 */
#define MAX_M 4096           /* maximum size of indices of Array u         */
#define M_DEFAULT 512        /* default size of indices of array u */
#define MAXSTEPS 500000      /* Maximum number of iterations               */
#define EPS     0.01        /* Numerical Tolerance */
#define EPS_DEFAULT 0.01    /* default accuracy */
#define PI 3.14159265       /* pi */
#define FILE_OUT 0          /* Output results to the file                  */
/* begin function prototyping  */

REAL **allocate_2D(int m, int n);
REAL my_max(REAL a, REAL b);
void init_array( INT m, INT n, REAL **a);
void bc( INT m, INT n, REAL **a, INT k, INT p );
INT write_file( INT m, INT n, REAL **u, INT k, INT p );
INT update_jacobi( INT m, INT n, REAL **u, REAL **unew, REAL *gdel);
INT update_sor( INT m, INT n, REAL **u, REAL **unew, REAL *gdel, REAL w);
REAL update_w(REAL* w, REAL ro);
INT replicate( INT m, INT n, REAL **u, REAL **ut );
INT transpose( INT m, INT n, REAL **u, REAL **ut );
void neighbors(INT k, INT p, INT UNDEFINED, INT *below, INT *above);
INT update_bc_2( INT mp, INT m, REAL **vt, INT k, INT below, INT above );
/* end function prototyping */


INT main(INT argc, CHAR *argv[]) {
/********** MAIN PROGRAM *********************************
 * Solve Laplace equation using Jacobi iteration method  *
        *
 *********************************************************/
      INT m = M_DEFAULT, mp, k, p, below, above;
      long iter=MAXSTEPS;
      REAL TOL=EPS_DEFAULT, del, gdel,start,finish,time;
      CHAR line[80];
      REAL **v, **vt, **vnew;

      MPI_Init(&argc, &argv);                       /* starts MPI */
      
      MPI_Comm_rank(MPI_COMM_WORLD, &k); /* get current process id */
      MPI_Comm_size(MPI_COMM_WORLD, &p); /* get # procs from env or */

      if(k == 0) {
          fprintf(OUTPUT," Number of args in command line, argc :\n");
          
          fprintf(OUTPUT," argc = %d :\n", argc);
          if (argc >= 2) 
              fprintf(OUTPUT," arg[1]= %s :\n", argv[1]);
          if (argc >= 3) 
              fprintf(OUTPUT," arg[2]= %s :\n", argv[2]);
          
          if (argc >1){
              m=atoi(argv[1]);
              m=MAX_M;
          }  
          
          if (argc >2 ) {
              TOL= atof(argv[2]);
              TOL=EPS;
          }
          
          fprintf(OUTPUT,"Size of interior points, m :\n");
          fprintf(OUTPUT,"m = %d\n",m);
          fprintf(OUTPUT,"Numerical accuracy, eps :\n");
          fprintf(OUTPUT,"eps = %f\n",TOL);

          fprintf(OUTPUT,"Number of processes, p :\n");
          fprintf(OUTPUT,"p = %d\n",p);
          
          start=-MPI_Wtime();
      }
      MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&TOL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      mp = m/p;
      
      v  = allocate_2D(m, mp);  /* allocate mem for 2D array */
      vt = allocate_2D(mp, m);
      vnew = allocate_2D(mp, m);

      gdel = 1.0;
      iter = 0;

      bc(m, mp, v, k, p); /* initialize and define B.C. for v */
      transpose(m, mp, v, vt);  /* solve for vt */
                                /* driven by need of update_bc_2 */

      replicate(mp, m, vt, vnew);            /* vnew = vt */
      neighbors(k, p, MPI_PROC_NULL, &below, &above);   /* domain borders */
      
      REAL ro = 1.0 - pow (PI / (2 * (m + 1)), 2);
      REAL w = 0;
      while (gdel > TOL) {  /* iterate until error below threshold */
        iter++;             /* increment iteration counter */

        if(iter > MAXSTEPS) {
          fprintf(OUTPUT,"Iteration terminated (exceeds %6d", MAXSTEPS);
          fprintf(OUTPUT," )\n");
          return (0);       /* nonconvergent solution */
        }
/* compute new solution according to the Jacobi scheme */
        /*update_jacobi(mp, m, vt, vnew, &del);  [> compute new vt <]*/
        update_sor(mp, m, vt, vnew, &del, w);
        update_w(&w, ro);
        if(iter%INCREMENT == 0) {
            gdel = 0;
          MPI_Allreduce( &del, &gdel, 1, MPI_DOUBLE,
               MPI_SUM, MPI_COMM_WORLD );  /* find global max error */
          if( k == 0) {
            fprintf(OUTPUT,"iter,del,gdel: %6d, %lf %lf\n",iter,del,gdel);
          }
        }
        update_bc_2( mp, m, vt, k, below, above); /* update b.c. */
      }
      if (k == 0) {
	  finish=MPI_Wtime();
        time=start+finish;
        fprintf(OUTPUT,"Stopped at iteration %d\n",iter);
        fprintf(OUTPUT,"The maximum error = %f\n",gdel);
        fprintf(OUTPUT,"Time = %f\n",time);

      }
      if (FILE_OUT) {
	  /* write v to file for use in MATLAB plots */
	  transpose(mp, m, vt, v);
	  write_file( m, mp, v, k, p );
	  
	  MPI_Barrier(MPI_COMM_WORLD);
      }

      free(v); free(vt); free(vnew); /* release allocated arrays  */
      MPI_Finalize();  
      return (0);
}
INT update_bc_2( INT mp, INT m, REAL **vt, INT k, INT below, INT above ) {
    MPI_Status status[6];  

    MPI_Sendrecv( vt[mp  ]+1, m, MPI_DOUBLE, above, 0,
		  vt[0   ]+1, m, MPI_DOUBLE, below, 0,
		  MPI_COMM_WORLD, status );

    MPI_Sendrecv( vt[1   ]+1, m, MPI_DOUBLE, below, 1,
		  vt[mp+1]+1, m, MPI_DOUBLE, above, 1,
		  MPI_COMM_WORLD, status );
    return (0);
}

REAL **allocate_2D(INT m, INT n) {
    INT i;
    REAL **a;

    a = (REAL **) malloc((unsigned) (m+2)*sizeof(REAL*));
    
/* Each pointer array element points to beginning of a row with n+2 entries*/
    for (i = 0; i <=m+1; i++) {
	a[i] = (REAL *) malloc((unsigned) (n+2)*sizeof(REAL));
    }
    
    return a;
}

INT write_file( INT m, INT n, REAL **u, INT k, INT p ) {
/*****************************************************
 * Writes 2D array ut columnwise (i.e. C convention) *
 * m   - size of rows m+2                            *
 * n   - size of columns n+2                         *
 * u   - scratch array                               *
 * k   - 0 <= k < p; = 0 for single thread code      *
 * p   - p >= 0; =1 for single thread code           *
 *****************************************************/
  INT ij, i, j, per_line;
  CHAR filename[50], file[53];
  FILE *fd;
/*
   prints u, 6 per line; used for matlab plots;
   PLOT_FILE contains the array size and number of procs;
   PLOT_FILE.(k+1) contains u pertaining to proc k;
   for serial job, PLOT_FILE.1 contains full u array.
*/

   (void) sprintf(filename, "%s", PLOT_FILE);

   if ( k == 0 ) {
      fd = fopen(filename, "w");
      fprintf(fd, "%5d %5d %5d\n", m+2, n+2, p);
      fclose(fd);
   }

   per_line = 6;                        /* to print 6 per line */
   (void) sprintf(file, "%s.%d", filename, k); /* create output file */
   fd = fopen(file, "w");
   ij = 0;
   for (j = 1; j <=n; j++) {
     for (i = 1; i <=m; i++) {
       fprintf(fd, "%11.4f ", u[i][j]);
     }
     fprintf(fd, "\n");
   }
   fprintf(fd, "\n");
   fclose(fd);
   return (0);
}

void init_array(INT m, INT n, REAL **a) {
/********* Initialize Array **********************
 * Initialize array with nx rows and ny columns  *
 *************************************************/
  INT i, j;

  for (i = 0; i <=m+1; i++) {
    for (j = 0; j <=n+1; j++) {
      a[i][j] = 0.0;          /* initialize all entries to zero */
    }
  }
}

void bc(INT m, INT n, REAL **u, INT k, INT p) {
/*********** Boundary Conditions ****************************************
 *  PDE: Laplacian u = 0;      0<=x<=1;  0<=y<=1                        *
 *  B.C.: u(x,0)=sin(pi*x); u(x,1)=sin(pi*x)*exp(-pi); u(0,y)=u(1,y)=0  *
 *  SOLUTION: u(x,y)=sin(pi*x)*exp(-pi*y)                               *
 ************************************************************************/
    INT i;

    init_array( m, n, u);                         /* initialize u to 0 */

  if (p > 1) {
    if (k == 0) {
      for (i = 0; i <=m+1; i++) {
        u[i][0] = sin(PI*i/(m+1));              /* at y = 0; all x */
      }
    }
    if (k == p-1) {
      for (i = 0; i <=m+1; i++) {
        u[i][n+1] = sin(PI*i/(m+1))*exp(-PI);   /* at y = 1; all x */
      }
    }
  } else if (p == 1) {
      for (i = 0; i <=m+1; i++) {
        u[i][  0] = sin(PI*i/(m+1));    /* at y = 0; all x */
        u[i][n+1] = u[i][0]*exp(-PI);   /* at y = 1; all x */
      }
  } else {
    printf("p is invalid\n");
  }
}

INT update_jacobi( INT m, INT n, REAL **u, REAL **unew, REAL *del) {
/**********************************************************************
 * Updates u according to Jacobi method                               *
 * m        - (INPUT)  size of interior rows                          *
 * n        - (INPUT)  size of interior columns                       *
 * u        - (INPUT)  solution array                                 *
 * unew     - (INPUT)  next solution array                            *
 * del      - (OUTPUT) error norm between 2 solution steps            *
 **********************************************************************/
  INT i, j;
  *del = 0.0;
  for (i = 1; i <=m; i++) {
    for (j = 1; j <=n; j++) {
      unew[i][j] = ( u[i  ][j+1] + u[i+1][j  ] +
                     u[i-1][j  ] + u[i  ][j-1] )*0.25;
      *del += fabs(unew[i][j] - u[i][j]);    /* find local max error */
    }
  }
  for (i = 1; i <=m; i++) {
    for (j = 1; j <=n; j++) {
      u[i][j] = unew[i][j];
    }
  }
      
  return (0);
}


INT replicate( INT m, INT n, REAL **a, REAL **b ) {
/********************************************************
 * Replicates array a into array b                      *
 * m    - (INPUT)  size of interior points in 1st index *
 * n    - (INPUT)  size of interior points in 2st index *
 * a    - (INPUT)  solution at time N                   *
 * b    - (OUTPUT) solution at time N + 1               *
 ********************************************************/
  INT i, j;

  for (i = 0; i <=m+1; i++) {
    for (j = 0; j <=n+1; j++) {
      b[i][j] = a[i][j];
    }
  }
  return (0);
}

INT transpose( INT m, INT n, REAL **a, REAL **at ) {
/**********************************************************
 * Transpose a(0:m+1,0:n+1) into at(0:n+1,0:m+1)          *
 * m    - (INPUT)  size of interior points in 1st index   *
 * n    - (INPUT)  size of interior points in 2st index   *
 * a    - (INPUT)  a  =  a(0:m+1,0:n+1)                   *
 * at   - (OUTPUT) at = at(0:n+1,0:m+1)                   *
 **********************************************************/
  INT i, j;

  for (i = 0; i <=m+1; i++) {
    for (j = 0; j <=n+1; j++) {
      at[j][i] = a[i][j];
    }
  }
  return (0);
}

void neighbors(INT k, INT p, INT UNDEFINED, INT *below, INT *above) {
/****************************************************************
 * determines two adjacent threads                              *
 * k         - (INPUT)  current thread                          *
 * p         - (INPUT)  number of processes (threads)           *
 * UNDEFINED - (INPUT)  code to assign to out-of-bound neighbor *
 * below     - (OUTPUT) neighbor thread below k (usually k-1)   *
 * above     - (OUTPUT) neighbor thread above k (usually k+1)   *
 ****************************************************************/
    *below = k != 0? k - 1 : UNDEFINED;
    *above = k != p - 1? k + 1 : UNDEFINED;
}

REAL update_w (REAL* w, REAL ro)
{
    REAL oldw = *w? *w : 1;

    *w = 1 / ( 1 - pow(ro, 2)*oldw / 4);
}


INT update_sor( INT m, INT n, REAL **u, REAL **unew, REAL *del, REAL w) {
  INT i, j;
  *del = 0.0;
  for (i = 1; i <=m; i++) {
    for (j = 1; j <=n; j++) {
      unew[i][j] = ( u[i  ][j+1] + u[i+1][j  ] +
                     u[i-1][j  ] + u[i  ][j-1] )*0.25;
      unew[i][j] = w * unew[i][j] + (1 - w) * u[i][j];
      *del += fabs(unew[i][j] - u[i][j]);    /* find local max error */
    }
  }
  for (i = 1; i <=m; i++) {
    for (j = 1; j <=n; j++) {
      u[i][j] = unew[i][j];
    }
  }
      
  return (0);
}
