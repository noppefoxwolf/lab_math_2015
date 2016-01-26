//
//  file0602.c
//  OpenGLTest
//
//  Created by Tomoya_Hirano on 6/2/15.
//  Copyright (c) 2015 Tomoya_Hirano. All rights reserved.
//

#include "main.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
#include <locale.h>

#define imax 20
#define jmax 20
#define kmax 20

double u[imax+1][jmax+1][kmax+1];
double uu[imax+1][jmax+1][kmax+1];
double v[imax+1][jmax+1][kmax+1];
double vv[imax+1][jmax+1][kmax+1];
double w[imax+1][jmax+1][kmax+1];
double ww[imax+1][jmax+1][kmax+1];

double P[imax+1][jmax+1][kmax+1];
double T[imax+1][jmax+1][kmax+1];
double TT[imax+1][jmax+1][kmax+1];

//include utility.c
void makeDir(char[100]);

void paud(int i,int j,int k,float dx, float dy, float dz, float dt, float Re, float Ra, float Pr);
void taud(int i,int j,int k,float dx, float dy, float dz, float dt, float Re, float Ra, float Pr);

void calcData(){
  struct tm *stm;
  time_t tim;
  char s[100];
  setlocale( LC_ALL, "jpn" );
  time( &tim );
  stm = localtime( &tim );
  strftime( s, 100, "%Y_%m_%d_%H_%M_%S", stm );
  char dirName[100] = {'\0'};
  snprintf(dirName, 100, "%s",s);
  makeDir(dirName);
  
  //レイノルズ数が大きいほど物質はさらさらする
  //風邪が吹くと物体の反対側で渦が生じ
  float Re, Ra, Pr;
  int i, j, k, n, m, N;					// 型の宣言時には、変数演算を用いて値を代入しながら宣言することはできない
  float dx, dy, dz, dt;						// したがってdxを用意した後dtを用意するとき、float dt=0.05*dx;とせず、float dt; dt=0.05*dx;とする
  double A, B;							// もしかしたら単純な数演算を用いて値を代入しながら宣言することはできるかも？つまり、int Re=10;は可能かも？
  float courseTime = 0.0;
  
  Re = 1000;
  Ra = 1000000.0;
  Pr = 5.0;
  dx = 1.0 / (float)imax;
  dy = 1.0 / (float)jmax;
  dz = 1.0 / (float)kmax;
  dt = 0.05 * dx;							// 毎回このdtずつ足していく、その行為がN回行われる、よって最終的なtは(N-1)*dt
  
  N = (int)(10 / dt)*2;						// 今は、結局N=200*imax=20000。ここで、int/(int)floatとしてしまうと分母0でエラー。

  //境界条件
  n = 0;
  for(i = 0; i <= imax; i++){
    for(j = 0; j <= jmax; j++){
      for (k = 0; k <= kmax; k++) {
        u[i][j][k] = 0;
        v[i][j][k] = 0;
        w[i][j][k] = 0;
        P[i][j][k] = 0;
        T[i][j][k] = j < (jmax/2.0) ? 0.5 : 1.5;
        if (j > (jmax/5*4) && i > (imax/5*2) && i < (imax/5*3)  && k > (kmax/5*2) && k < (kmax/5*3)) {
          T[i][j][k] = 0;
        }
      }
    }
  }
  
  
  // u[i][j],v[i][j],P[i][j]の計算と書き込み、時刻n*dt(n=1,2,...,N-1)のとき
  // 課題ノートの、[2]~[4]の繰り返し
  for(n = 1; n < N; n++) {
    time_t start_time;
    start_time = time(NULL);
    
    
    for(i = 0; i <= imax; i++){
      for(j = 0; j <= jmax; j++){
        for (k = 0; k <= kmax; k++) {
          uu[i][j][k] = u[i][j][k];			// 全ての[i][j]で、uをコピーしておく
          vv[i][j][k] = v[i][j][k];
          ww[i][j][k] = w[i][j][k];
          TT[i][j][k] = T[i][j][k];
        }
      }
    }
    
    //３次精度
    for (i=2; i<= imax-2; i++) {
      for (j=2; j<= jmax-2; j++) {
        for (k=2; k<= kmax-2; k++) {
          taud(i, j, k, dx, dy, dz, dt, Re, Ra, Pr);
        }
      }
    }
    
    
    //一次精度
    for (i=1; i< imax; i++) {
      for (j=1; j< jmax; j+=jmax-2) {//上下の部分を計算
        for (k=1; k< kmax; k++) {
          paud(i, j, k, dx, dy, dz, dt, Re, Ra, Pr);
        }
      }
    }
    
    for (i=1;i<imax; i+=imax-2){
      for (j=1; j<jmax; j++) {//横の部分を計算
        for (k=1; k<kmax; k++) {
          paud(i, j, k, dx, dy, dz, dt, Re, Ra, Pr);
        }
      }
    }
    
    for (i=1;i<imax; i++) {
      for (j=1; j<jmax; j++) {
        for (k=1; k<kmax; k+=kmax-2) {
          paud(i, j, k, dx, dy, dz, dt, Re, Ra, Pr);
        }
      }
    }
    
    //ポアソン方程式
    for(m = 0; m < 100; m++){
      for(i = 1; i < imax; i++){
        for(j = 1; j < jmax; j++){//ここにもkを足す
          for(k = 1; k < kmax; k++){//ここにもkを足す
            A =
            - (1.0 * ( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx) * ( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx))
            - (1.0 * ( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy) * ( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy))
            - (1.0 * ( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz) * ( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz))
            
            - (2.0 * ( u[i][j+1][k] - u[i][j-1][k] ) / (2.0*dy)
                * ( v[i+1][j][k] - v[i-1][j][k] ) / (2.0*dx))
            
            - (2.0 * ( u[i][j][k+1] - u[i][j][k-1] ) / (2.0*dz)
                * ( w[i+1][j][k] - w[i-1][j][k] ) / (2.0*dx))
            
            - (2.0 * ( v[i][j][k+1] - v[i][j][k-1] ) / (2.0*dz)
                * ( w[i][j+1][k] - w[i][j-1][k] ) / (2.0*dy))
            
            + (1.0/dt) * (
                  (( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx)) +
                  (( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy)) +
                  (( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz))
                  );
            //問題無さそう
            B = - (( P[i+1][j][k] + P[i-1][j][k] ) / (dx*dx))
              - (( P[i][j+1][k] + P[i][j-1][k] ) / (dy*dy))
              - (( P[i][j][k+1] + P[i][j][k-1] ) / (dz*dz)) + A;
            
            P[i][j][k] = B / ( - (2.0 / (dx*dx)) - (2.0 / (dy*dy)) - (2.0 / (dz*dz)) );
          }
        }
      }
    }
    
    // 境界条件
    for(j = 0; j <= jmax; j++){
      for(k = 0; k <= kmax; k++){
        P[0][j][k] = P[1][j][k];
        P[imax][j][k] = P[imax-1][j][k];
        T[0][j][k] = T[1][j][k];
        T[imax][j][k] = T[imax-1][j][k];
      }
    }
    for(i = 0; i <= imax; i++){
      for(k = 0; k <= kmax; k++){
        P[i][0][k] = P[i][1][k];
        P[i][jmax][k] = P[i][jmax-1][k];
        T[i][0][k] = T[i][1][k];
        T[i][jmax][k] = T[i][jmax-1][k];
      }
    }
    for(i = 0; i <= imax; i++){
      for(j = 0; j <= jmax; j++){
        P[i][j][0] = P[i][j][1];
        P[i][j][kmax] = P[i][j][kmax-1];
        T[i][j][0] = T[i][j][1];
        T[i][j][kmax] = T[i][j][kmax-1];
      }
    }
    
    if (n%20==0) {
      char filename[100] = {'\0'};
      snprintf(filename, 100, "%s/paraview.csv.%d",dirName,n);
      FILE *fp = fopen(filename, "w");
      fprintf(fp, "t,x,y,z,T,\n");
      
      // 時刻n*dt(第n回目時点)の書き込み
      for(i = 0; i <= imax; i++){
        for(j = 0; j <= jmax; j++){
          for (k = 0; k <= kmax; k++) {
            fprintf(fp, "%f,%f,%f,%f,%f,\n", n*dt,i*dx,j*dy,k*dz,T[i][j][k]);
            //          printf("i(%d) : j(%d) : k(%d) :: %f\n",i,j,k,T[i][j][k]);
            if (isnan(T[i][j][k])){
              printf("nan value ...\n");
              exit(0);
            }else if (isinf(T[i][j][k])){
              printf("inf value ...\n");
              exit(0);
            }
          }
        }
      }
      fclose(fp);
    }
    
    //進捗報告
    float diffTime = difftime(time(NULL), start_time);
    courseTime += diffTime;
    int second = (int)((float)(N-(n+1))*(courseTime/(float)n));
    int time   = second / 3600;
    int minute = (second - time * 3600) / 60;
    second = second % 60;
    printf("%f％ 残り%d時間%d分%d秒\n", (float)(n+1)/(float)N,time,minute,second);
  }
  printf("終了\n");
}

//三次精度上流差分
void taud(int i,int j,int k,float dx, float dy, float dz, float dt, float Re, float Ra, float Pr){
  //P183の7.75
  float ududx = uu[i][j][k] * (-uu[i+2][j][k]+(8*(uu[i+1][j][k]-uu[i-1][j][k]))+uu[i-2][j][k])/(12*dx)
    + ((fabs(uu[i][j][k])*1.0/12.0)*(uu[i+2][j][k]-(4*uu[i+1][j][k])+(6*uu[i][j][k])-(4*uu[i-1][j][k])+uu[i-2][j][k])/dx);
  float vdudy = vv[i][j][k] * (-uu[i][j+2][k]+(8*(uu[i][j+1][k]-uu[i][j-1][k]))+uu[i][j-2][k])/(12*dy)
    + ((fabs(vv[i][j][k])*1.0/12.0)*(uu[i][j+2][k]-(4*uu[i][j+1][k])+(6*uu[i][j][k])-(4*uu[i][j-1][k])+uu[i][j-2][k])/dy);
  float wdudz = ww[i][j][k] * (-uu[i][j][k+2]+(8*(uu[i][j][k-1]-uu[i][j][k-1]))+uu[i][j][k-2])/(12*dz)
    + ((fabs(ww[i][j][k])*1.0/12.0)*(uu[i][j][k+2]-(4*uu[i][j][k+1])+(6*uu[i][j][k])-(4*uu[i][j][k-1])+uu[i][j][k-2])/dz);
  float dPdx = (P[i+1][j][k]-P[i-1][j][k])/(2.0*dx);
  float dx2 = dx*dx;
  float dy2 = dy*dy;
  float dz2 = dz*dz;
  
  float d2udx2 = ( uu[i+1][j][k] - (2.0*uu[i][j][k]) + uu[i-1][j][k] ) / dx2;
  float d2udy2 = ( uu[i][j+1][k] - (2.0*uu[i][j][k]) + uu[i][j-1][k] ) / dy2;
  float d2udz2 = ( uu[i][j][k+1] - (2.0*uu[i][j][k]) + uu[i][j][k-1] ) / dz2;
  u[i][j][k] = uu[i][j][k] + (dt*( -ududx -vdudy -wdudz - dPdx + (1.0/Re*(d2udx2+d2udy2+d2udz2)) ));
  
  float udvdx = uu[i][j][k] * (-vv[i+2][j][k]+(8*(vv[i+1][j][k]-vv[i-1][j][k]))+vv[i-2][j][k])/(12*dx)
    + ((fabs(uu[i][j][k])*1.0/12.0)*(vv[i+2][j][k]-(4*vv[i+1][j][k])+(6*vv[i][j][k])-(4*vv[i-1][j][k])+vv[i-2][j][k])/dx);
  float vdvdy = vv[i][j][k] * (-vv[i][j+2][k]+(8*(vv[i][j+1][k]-vv[i][j-1][k]))+vv[i][j-2][k])/(12*dy)
    + ((fabs(vv[i][j][k])*1.0/12.0)*(vv[i][j+2][k]-(4*vv[i][j+1][k])+(6*vv[i][j][k])-(4*vv[i][j-1][k])+vv[i][j-2][k])/dy);
  float wdvdz = ww[i][j][k] * (-vv[i][j][k+2]+(8*(vv[i][j][k+1]-vv[i][j][k-1]))+vv[i][j][k-2])/(12*dz)
    + ((fabs(ww[i][j][k])*1.0/12.0)*(vv[i][j][k+2]-(4*vv[i][j][k+1])+(6*vv[i][j][k])-(4*vv[i][j][k-1])+vv[i][j][k-2])/dz);
  float dPdy = (P[i][j+1][k]-P[i][j-1][k])/(2.0*dy);
  float d2vdx2 = ( vv[i+1][j][k] - (2.0*vv[i][j][k]) + vv[i-1][j][k] ) / dx2;
  float d2vdy2 = ( vv[i][j+1][k] - (2.0*vv[i][j][k]) + vv[i][j-1][k] ) / dy2;
  float d2vdz2 = ( vv[i][j][k+1] - (2.0*vv[i][j][k]) + vv[i][j][k-1] ) / dz2;
  //浮力項
  v[i][j][k] = vv[i][j][k] + (dt*( -udvdx -vdvdy -wdvdz -dPdy + (1.0/Re*(d2vdx2+d2vdy2+d2vdz2)) + ((Ra/(Re*Re*Pr))*T[i][j][k]) ));
  
  float udwdx = uu[i][j][k] * (-ww[i+2][j][k]+(8*(ww[i+1][j][k]-ww[i-1][j][k]))+ww[i-2][j][k])/(12*dx)
    + ((fabs(uu[i][j][k])*1.0/12.0)*(ww[i+2][j][k]-(4*ww[i+1][j][k])+(6*ww[i][j][k])-(4*ww[i-1][j][k])+ww[i-2][j][k])/dx);
  float vdwdy = vv[i][j][k] * (-ww[i][j+2][k]+(8*(ww[i][j+1][k]-ww[i][j-1][k]))+ww[i][j-2][k])/(12*dy)
    + ((fabs(vv[i][j][k])*1.0/12.0)*(ww[i][j+2][k]-(4*ww[i][j+1][k])+(6*ww[i][j][k])-(4*ww[i][j-1][k])+ww[i][j-2][k])/dy);
  float wdwdz = ww[i][j][k] * (-ww[i][j][k+2]+(8*(ww[i][j][k+1]-ww[i][j][k-1]))+ww[i][j][k-2])/(12*dz)
    + ((fabs(ww[i][j][k])*1.0/12.0)*(ww[i][j][k+2]-(4*ww[i][j][k+1])+(6*ww[i][j][k])-(4*ww[i][j][k-1])+ww[i][j][k-2])/dz);
  float dPdz = (P[i][j][k+1]-P[i][j][k-1])/(2.0*dz);
  float d2wdx2 = ( ww[i+1][j][k] - (2.0*ww[i][j][k]) + ww[i-1][j][k] ) / dx2;
  float d2wdy2 = ( ww[i][j+1][k] - (2.0*ww[i][j][k]) + ww[i][j-1][k] ) / dy2;
  float d2wdz2 = ( ww[i][j][k+1] - (2.0*ww[i][j][k]) + ww[i][j][k-1] ) / dz2;
  w[i][j][k] = ww[i][j][k] + (dt*( -udwdx -vdwdy -wdwdz -dPdz + (1.0/Re*(d2wdx2+d2wdy2+d2wdz2)) ));
  
  
  //上流差分
  float udTdx  = uu[i][j][k] * (-TT[i+2][j][k]+(8*(TT[i+1][j][k]-TT[i-1][j][k]))+TT[i-2][j][k])/(12*dx)
    + ((fabs(uu[i][j][k])*1.0/12.0)*(TT[i+2][j][k]-(4*TT[i+1][j][k])+(6*TT[i][j][k])-(4*TT[i-1][j][k])+TT[i-2][j][k])/dt);
  float vdTdy  = vv[i][j][k] * (-TT[i][j+2][k]+(8*(TT[i][j+1][k]-TT[i][j-1][k]))+TT[i][j-2][k])/(12*dy)
    + ((fabs(vv[i][j][k])*1.0/12.0)*(TT[i][j+2][k]-(4*TT[i][j+1][k])+(6*TT[i][j][k])-(4*TT[i][j-1][k])+TT[i][j-2][k])/dy);
  float wdTdy  = ww[i][j][k] * (-TT[i][j][k+2]+(8*(TT[i][j][k+1]-TT[i][j][k-1]))+TT[i][j][k-2])/(12*dz)
    + ((fabs(ww[i][j][k])*1.0/12.0)*(TT[i][j][k+2]-(4*TT[i][j][k+1])+(6*TT[i][j][k])-(4*TT[i][j][k-1])+TT[i][j][k-2])/dz);
  float d2Tdx2 = ( TT[i+1][j][k] - (2.0*TT[i][j][k]) + TT[i-1][j][k] ) / dx2;
  float d2Tdy2 = ( TT[i][j+1][k] - (2.0*TT[i][j][k]) + TT[i][j-1][k] ) / dy2;
  float d2Tdz2 = ( TT[i][j][k+1] - (2.0*TT[i][j][k]) + TT[i][j][k-1] ) / dz2;
  
  T[i][j][k] = TT[i][j][k] + (dt * ( -udTdx - vdTdy -wdTdy + ((1.0/(Re*Pr))*(d2Tdx2 + d2Tdy2 + d2Tdz2)) ));
}

//一次精度上流差分
void paud(int i,int j,int k,float dx, float dy, float dz, float dt, float Re, float Ra, float Pr){
  u[i][j][k] = uu[i][j][k]
  + (dt * ( - uu[i][j][k] * ( uu[i+1][j][k] - uu[i-1][j][k] ) / (2.0*dx)
           + (fabs(uu[i][j][k]) / 2.0 * ( uu[i+1][j][k] - (2.0*uu[i][j][k]) + uu[i-1][j][k] ) / dx)
           - (vv[i][j][k] * ( uu[i][j+1][k] - uu[i][j-1][k] ) / (2.0*dy))
           + (fabs(vv[i][j][k]) / 2.0 * ( uu[i][j+1][k] - (2.0*uu[i][j][k]) + uu[i][j-1][k] ) / dy)
           - (ww[i][j][k] * ( uu[i][j][k+1] - uu[i][j][k-1] ) / (2.0*dz))
           + (fabs(ww[i][j][k]) / 2.0 * ( uu[i][j][k+1] - (2.0*uu[i][j][k]) + uu[i][j][k-1] ) / dz)
           - (( P[i+1][j][k] - P[i-1][j][k] ) / (2.0*dx))
           + ((1.0/Re) * (  (( uu[i+1][j][k] - (2.0*uu[i][j][k]) + uu[i-1][j][k] ) / (dx*dx))
                          + (( uu[i][j+1][k] - (2.0*uu[i][j][k]) + uu[i][j-1][k] ) / (dy*dy))
                          + (( uu[i][j][k+1] - (2.0*uu[i][j][k]) + uu[i][j][k-1] ) / (dz*dz))
                          ))));
  //非線形項の確認
  //上流差分の式が二行に渡っている部分
  v[i][j][k] = vv[i][j][k]
  + (dt * ( - uu[i][j][k] * ( vv[i+1][j][k] - vv[i-1][j][k] ) / (2.0*dx)
           + (fabs(uu[i][j][k]) / 2.0 * ( vv[i+1][j][k] - (2.0*vv[i][j][k]) + vv[i-1][j][k] ) / dx)
           - (vv[i][j][k] * ( vv[i][j+1][k] - vv[i][j-1][k] ) / (2.0*dy))
           + (fabs(vv[i][j][k]) / 2.0 * ( vv[i][j+1][k] - (2.0*vv[i][j][k]) + vv[i][j-1][k] ) / dy)
           - (ww[i][j][k] * ( vv[i][j][k+1] - vv[i][j][k-1] ) / (2.0*dz))
           + (fabs(ww[i][j][k]) / 2.0 * ( vv[i][j][k+1] - (2.0*vv[i][j][k]) + vv[i][j][k-1] ) / dz)
           - (( P[i][j+1][k] - P[i][j-1][k] ) / (2.0*dy))
           + ((1.0/Re) * (  (( vv[i+1][j][k] - (2.0*vv[i][j][k]) + vv[i-1][j][k] ) / (dx*dx))
                          + (( vv[i][j+1][k] - (2.0*vv[i][j][k]) + vv[i][j-1][k] ) / (dy*dy))
                          + (( vv[i][j][k+1] - (2.0*vv[i][j][k]) + vv[i][j][k-1] ) / (dz*dz))
                          ))
           + ( (Ra/(Re*Re*Pr)) * T[i][j][k] )
           ));
  
  w[i][j][k] = ww[i][j][k]
  + (dt * ( - uu[i][j][k] * ( ww[i+1][j][k] - ww[i-1][j][k] ) / (2.0*dx)
           + (fabs(uu[i][j][k]) / 2.0 * ( ww[i+1][j][k] - (2.0*ww[i][j][k]) + ww[i-1][j][k] ) / dx)
           - (vv[i][j][k] * ( ww[i][j+1][k] - ww[i][j-1][k] ) / (2.0*dy))
           + (fabs(vv[i][j][k]) / 2.0 * ( ww[i][j+1][k] - (2.0*ww[i][j][k]) + ww[i][j-1][k] ) / dy)
           - (ww[i][j][k] * ( ww[i][j][k+1] - ww[i][j][k-1] ) / (2.0*dz))
           + (fabs(ww[i][j][k]) / 2.0 * ( ww[i][j][k+1] - (2.0*ww[i][j][k]) + ww[i][j][k-1] ) / dz)
           - (( P[i][j][k+1] - P[i][j][k-1] ) / (2.0*dz))
           + ((1.0/Re) * (  (( ww[i+1][j][k] - (2.0*ww[i][j][k]) + ww[i-1][j][k] ) / (dx*dx))
                          + (( ww[i][j+1][k] - (2.0*ww[i][j][k]) + ww[i][j-1][k] ) / (dy*dy))
                          + (( ww[i][j][k+1] - (2.0*ww[i][j][k]) + ww[i][j][k-1] ) / (dz*dz))
                          ))
           ));
  
  //ここにもTが必要
  T[i][j][k] = TT[i][j][k]
  + (dt * (- uu[i][j][k] * ( TT[i+1][j][k] - TT[i-1][j][k] ) / (2.0*dx)
           + (fabs(uu[i][j][k]) / 2.0 * ( TT[i+1][j][k] - (2.0*TT[i][j][k]) + TT[i-1][j][k] ) / dx)
           - (vv[i][j][k] * ( TT[i][j+1][k] - TT[i][j-1][k] ) / (2.0*dy))
           + (fabs(vv[i][j][k]) / 2.0 * ( TT[i][j+1][k] - (2.0*TT[i][j][k]) + TT[i][j-1][k] ) / dy)
           - (ww[i][j][k] * ( TT[i][j][k+1] - TT[i][j][k-1] ) / (2.0*dz))
           + (fabs(ww[i][j][k]) / 2.0 * ( TT[i][j][k+1] - (2.0*TT[i][j][k]) + TT[i][j][k-1] ) / dz)
           + ((1.0/(Re*Pr)) * (  (( TT[i+1][j][k] - (2.0*TT[i][j][k]) + TT[i-1][j][k] ) / (dx*dx))
                               + (( TT[i][j+1][k] - (2.0*TT[i][j][k]) + TT[i][j-1][k] ) / (dy*dy))
                               + (( TT[i][j][k+1] - (2.0*TT[i][j][k]) + TT[i][j][k-1] ) / (dz*dz))
                               ))
           ));
}

