/* *****************************************************************
  PROGRAM: forwardalg.c
AUTHOR:  Andrew Titman
DATE:    February 2009

Forward algorithm

******************************************************************  */

  #include <R.h>
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
  #include <Rmath.h>


int MI(
  int i,
  int j,
  int nrow
)
{
  return j * nrow + i;
}

int MI2(
  int i,
  int j,
  int k,
  int n1,
  int n2
)
{
  return k*(n1*n2) + j*n1 + i;
}

/*Version for 4 dimensional arrays*/
  int MI3(
    int i,
    int j,
    int k,
    int l,
    int n1,
    int n2,
    int n3
  )
{
  return l*(n1*n2*n3) + k*(n1*n2) + j*n1 + i;
  }



  void nhm_forwardalg(
    double *init,
    double *dinit,
    double *e,
    double *de,
    double *alp,
    double *dalp,
    double *pmat,
    double *dpmat,
    int *nob,
    int *nlst,
    int *npar,
    double *prob,
    double *dprob
  )
    /* Here e is a nlst x nob vector and de a nlst x nob x npar vector */
    /* This should allow a) Censoring and b) Exact death times to be passed easily.*/
    /* Also permits covariates and time dependendent covariates*/
    /* Moreover, allows more general hidden Markov model structures */
    {
      int i,j,k,l;
      /*Set initial values to their initial values */
        for (i=0;i<=nlst[0];++i) {
          alp[MI(0,i,nob[0]+1)] = init[i];
          for (j=0;j<=npar[0];++j) {
            dalp[MI2(0,i,j,nob[0]+1,nlst[0]+1)] = dinit[MI(i,j,nlst[0]+1)];
          }
        }
      for (k=1;k<=nob[0];++k) {
        /*Run through each observation */
          /*Set to zero first*/
          for (i=0;i<=nlst[0];++i) {
            alp[MI(k,i,nob[0]+1)]= 0;
            for (l=0;l<=npar[0];++l) {
              dalp[MI2(k,i,l,nob[0]+1,nlst[0]+1)]=0;
            }
            for (j=0;j<=nlst[0];++j) {
              alp[MI(k,i,nob[0]+1)]+= alp[MI(k-1,j,nob[0]+1)]*e[MI(i,k-1,nlst[0]+1)]*pmat[MI2(j,i,k-1,nlst[0]+1,nlst[0]+1)];


              for (l=0;l<=npar[0];++l) {
                dalp[MI2(k,i,l,nob[0]+1,nlst[0]+1)]+= dalp[MI2(k-1,j,l,nob[0]+1,nlst[0]+1)]*e[MI(i,k-1,nlst[0]+1)]*pmat[MI2(j,i,k-1,nlst[0]+1,nlst[0]+1)] + alp[MI(k-1,j,nob[0]+1)]*de[MI2(i,k-1,l,nlst[0]+1,nob[0])]*pmat[MI2(j,i,k-1,nlst[0]+1,nlst[0]+1)] + alp[MI(k-1,j,nob[0]+1)]*e[MI(i,k-1,nlst[0]+1)]*dpmat[MI3(j,i,k-1,l,nlst[0]+1,nlst[0]+1,nob[0])];
              }

            }
          }
      }
      prob[0]=0;
      for (l=0;l<=npar[0];++l) {
        dprob[l]=0;
      }

      for (i=0;i<=nlst[0];++i) {
        prob[0] += alp[MI(nob[0],i,nob[0]+1)];

        for (l=0;l<=npar[0];++l) {
          dprob[l] += dalp[MI2(nob[0],i,l,nob[0]+1,nlst[0]+1)];
        }

      }
    }
