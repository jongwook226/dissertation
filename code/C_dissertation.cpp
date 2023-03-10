#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "C_functions.h"

//spherical distance function
/*
void SphiDist(double* lat1, double* long1, double* lat2, double* long2, double* dist){
  double** mdat = (double**) malloc(sizeof(double*) * *ncol);
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;
  dist* = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(long2 - long1)));
}
*/

void SphDist(int* nrow, int* ncol, double* dat, double* sDmat){
  double** mdat = (double**) malloc(sizeof(double*) * *ncol);
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

  double** msDmat = (double**) malloc(sizeof(double*) * *nrow);
  for(int i = 0; i < *nrow; i++) msDmat[i] = sDmat + i * *nrow;

  for(int i=0; i<*nrow; i++){
    for(int j=i; j<*nrow; j++){
      double dist = acos(sin(mdat[0][i])*sin(mdat[0][j]) + cos(mdat[0][i])*cos(mdat[0][j])*cos(mdat[1][i] - mdat[1][j]));
      if(isnan(dist)){
        msDmat[i][j] = msDmat[j][i] = 0;
      }else{
        msDmat[i][j] = msDmat[j][i] = dist;
      }
    }
  }
  free(mdat);
  free(msDmat);
}


//Function to compute covariance matrix
void Cov_mat(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, double* a1, double* a2, double* a3, int* kappa, double* q, int* ncol_q, double* tauDist, double* bwtauDist, double* R){

  const double pi = 3.141592653589793238462643383279502884L; //define pi
  double sDist = 0; //spherical diatance value
  double tDist = 0; //time lag value

  double** mdat = (double**) malloc(sizeof(double*) * *ncol); //dataset
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

  double** mDmat = (double**) malloc(sizeof(double*) * *nrow_Dmat); //spherical distance matrix
  for(int i = 0; i < *nrow_Dmat; i++) mDmat[i] = sDmat + i * *nrow_Dmat;

  double** mq = (double**) malloc(sizeof(double*) * *ncol_q);  //lower order of spherical harmonics
  for(int i = 0; i < *ncol_q; i++) mq[i] = q + i * *nrow; //nrow of low_spherical is the same as nrow of dat

  double** mtauDist = (double**) malloc(sizeof(double*) * *ncol_q);  //distance between taus and dataset
  for(int i = 0; i < *ncol_q; i++) mtauDist[i] = tauDist + i * *nrow; //nrow of tauDist is the same as nrow of dat

  double** mbwtauDist = (double**) malloc(sizeof(double*) * *ncol_q);  //distance between taus
  for(int i = 0; i < *ncol_q; i++) mbwtauDist[i] = bwtauDist + i * *ncol_q; //ncol_low_spherical is kappa^2 (the number of taus)

  double** mR = (double**) malloc(sizeof(double*) * *nrow); //nrow and ncol of mR(covariance matrix) are the same as nrow of dat
  for(int i = 0; i < *nrow; i++) mR[i] = R + i * *nrow;


  for(int i=0; i<*nrow; i++){
    for(int j=i; j<*nrow; j++){
      //spherical distance of a pair by index (distancce bw P and Q)
      sDist = mDmat[(int) mdat[3][i]][(int) mdat[3][j]]; //mdat[3][] is index of location

      //time lag of a pair
      tDist = abs(mdat[2][i] - mdat[2][j]); //mdat[2][] is time

      if(*kappa==0){
        double Ct = a1[0] * exp(-(a2[0]) * abs(tDist)); //a_l(h)
        mR[i][j] = mR[j][i] = a3[0] * (1 - Ct*Ct)/((1-2*cos(sDist)*Ct+Ct*Ct) * sqrt((1-2*cos(sDist)*Ct+Ct*Ct)));
      }else if(*kappa==1){
        //save pairwise covariances
        double Ct1 = a1[0] * exp(-(a2[0]) * abs(tDist)); //a_l(h)
        double Ct2 = a1[1] * exp(-(a2[1]) * abs(tDist)); //for b_0(h)
        double Ct3 = a1[2] * exp(-(a2[2]) * abs(tDist)); //for b_0^{l,m}(h)

        double icf1 = (1 - Ct1*Ct1)/((1-2*cos(sDist)*Ct1+Ct1*Ct1) * sqrt((1-2*cos(sDist)*Ct1+Ct1*Ct1))); //bw Truncated process
        double icf2 = (1 - Ct2*Ct2)/((1-2*Ct2+Ct2*Ct2) * sqrt((1-2*Ct2+Ct2*Ct2))); //bw Nil space
        double icf3 = (1 - Ct3*Ct3)/((1-2*cos(mtauDist[0][i])*Ct3+Ct3*Ct3) * sqrt((1-2*cos(mtauDist[0][i])*Ct3+Ct3*Ct3))); //mdat[3][] is tauDist
        double icf4 = (1 - Ct3*Ct3)/((1-2*cos(mtauDist[0][j])*Ct3+Ct3*Ct3) * sqrt((1-2*cos(mtauDist[0][j])*Ct3+Ct3*Ct3))); //mdat[3][] is tauDist

        //mR[i][j] = mR[j][i] = a3[0] * icf1 + a3[1] * 1/(4*pi) * icf2 + a3[2] * ( 1/(2*sqrt(pi)) * icf3 + a3[2] * 1/(2*sqrt(pi)) * icf4 - a3[0]*1/(4*pi) - a3[1]*1/(16*pi*pi)) - a3[2]*1/(4*pi*sqrt(pi));
        //mR[i][j] = mR[j][i] = a3[0] * icf1 + a3[1] * 1/(4*pi) * icf2 + a3[2] * ( 1/(2*sqrt(pi)) * icf3 + a3[2] * 1/(2*sqrt(pi)) * icf4 - 1/(4*pi) - 1/(16*pi*pi)) - 1/(4*pi*sqrt(pi));
        //mR[i][j] = mR[j][i] = a3[0]*(icf1 + 1/(4*pi) * icf2 + (1/(2*sqrt(pi)) * icf3 + 1/(2*sqrt(pi)) * icf4 - 1/(4*pi) - 1/(16*pi*pi)) - 1/(4*pi*sqrt(pi)));
        mR[i][j] = mR[j][i] = (a3[0]*icf1) + (a3[1] * (icf2 + 1)) - (a3[2] * icf3) - (a3[2] * icf4);

      }else if(*kappa==2){

        double part1=0; //icf part
        double part2=0; //nil space part
        double part3=0; //b/w nil and truncated part
        double part4=0; //b/w nil and truncated part

        double Ct1 = a1[0] * exp(-(a2[0]) * abs(tDist)); //a_l(h)
        double Ct2 = a1[1] * exp(-(a2[1]) * abs(tDist)); //for b_0(h)
        double Ct3 = a1[2] * exp(-(a2[2]) * abs(tDist)); //for b_0^{l,m}(h)

        //low order spherical harmonics
        double q2=0;
        double q3=0;
        double q4=0;

        //bw Truncated process (homogeneous part)
        double icf1 = (1 - Ct1*Ct1)/((1-2*cos(sDist)*Ct1+Ct1*Ct1) * sqrt((1-2*cos(sDist)*Ct1+Ct1*Ct1)));
        double trunc1 = 1/(4*pi) + 3/(4*pi)*Ct1*sDist; // when ell=0 and 1 for Nil space
        part1 = a3[0]*(icf1 - trunc1);

        ////bw Nil space
        for(int k=0; k<*ncol_q; k++){
          for(int l=0; l<*ncol_q; l++){
            double icf2 = (1 - Ct2*Ct2)/((1-2*cos(mbwtauDist[k][l])*Ct2+Ct2*Ct2) * sqrt((1-2*cos(mbwtauDist[k][l])*Ct2+Ct2*Ct2)));
            double trunc2 = 1/(4*pi) + 3/(4*pi)*Ct2*mbwtauDist[k][l];
            q2 = mq[k][i]*mq[l][j];

            part2 += (a3[1] * (icf2 - trunc2) * q2);
            if(k == l){
              part2 += a3[1] * (mq[k][i] * mq[k][j]); //When q_nu = q_mu
            }
          }
        }

        //between Nil space and truncated process
        for(int k=0; k<*ncol_q; k++){
          double icf3 = (1 - Ct3*Ct3)/((1-2*cos(mtauDist[k][i])*Ct3+Ct3*Ct3) * sqrt((1-2*cos(mtauDist[k][i])*Ct3+Ct3*Ct3))); //mtauDist[3][i] is distance bw P and tau
          double icf4 = (1 - Ct3*Ct3)/((1-2*cos(mtauDist[k][j])*Ct3+Ct3*Ct3) * sqrt((1-2*cos(mtauDist[k][j])*Ct3+Ct3*Ct3))); //mtauDist[3][i] is distance bw Q and tau

          double trunc3 = 1/(4*pi) + 3/(4*pi)*Ct3*mtauDist[k][i];
          double trunc4 = 1/(4*pi) + 3/(4*pi)*Ct3*mtauDist[k][j];

          q3 = mq[k][j];
          q4 = mq[k][i];

          part3 += a3[2] * (icf3 - trunc3) * q3;
          part4 += a3[2] * (icf4 - trunc4) * q4;
        }

        mR[i][j] = mR[j][i] = part1 + part2 - part3 - part4;
      }
    }
  }
  free(mdat);
  free(mDmat);
  free(mq);
  free(mR);
}

//Function for MOM estimator
void G_hat(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, int* length_Hvec, double* Hvec, double* eps, int* t, double* G, double* N){


  double sDist = 0; //spherical diatance value
  double tDist = 0; //time lag value

  double** mdat = (double**) malloc(sizeof(double*) * *ncol);
  for(int i = 0; i < *ncol; i++) mdat[i] = dat + i * *nrow;

  double** mDmat = (double**) malloc(sizeof(double*) * *nrow_Dmat);
  for(int i = 0; i < *nrow_Dmat; i++) mDmat[i] = sDmat + i * *nrow_Dmat;

  double** mG = (double**) malloc(sizeof(double*) * *t);
  for(int i = 0; i < *t; i++) mG[i] = G + i * *length_Hvec;

  double** mN = (double**) malloc(sizeof(double*) * *t);
  for(int i = 0; i < *t; i++) mN[i] = N + i * *length_Hvec;


  for(int i = 0; i < *nrow; i++){
    for(int j = i; j < *nrow; j++){

     //spherical distance of a pair by index
     sDist = mDmat[(int) mdat[3][i]][(int) mdat[3][j]];

     //time lag of a pair
     tDist = abs(mdat[2][i] - mdat[2][j]);

     for(int k = 0; k < *t; k++){
      for(int l = 0; l < *length_Hvec; l++){
        if((Hvec[l] - *eps) < sDist & sDist <= (Hvec[l] + *eps) & tDist == k){
          //#number of obs
          mN[k][l] += 1;
          //#sum of covariances to compute MOM estimate
          mG[k][l] += mdat[4][i] * mdat[4][j];
        }
      }
     }
    }
  }
  free(mdat);
  free(mDmat);
  free(mG);
  free(mN);
}
