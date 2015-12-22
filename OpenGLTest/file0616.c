//
//  file0602.c
//  OpenGLTest
//
//  Created by Tomoya_Hirano on 6/2/15.
//  Copyright (c) 2015 Tomoya_Hirano. All rights reserved.
//

//#include <string.h>
//#include <math.h>
#include "myglini.h"
//#define N 8
//#define RHO 0.4

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <dirent.h>
#include <sys/stat.h>

#include <time.h>
#include <locale.h>

#define imax 50
#define jmax 50
#define kmax 50

double u[imax+1][jmax+1][kmax+1];
double uu[imax+1][jmax+1][kmax+1];
double v[imax+1][jmax+1][kmax+1];
double vv[imax+1][jmax+1][kmax+1];
double w[imax+1][jmax+1][kmax+1];
double ww[imax+1][jmax+1][kmax+1];

double P[imax+1][jmax+1][kmax+1];
double T[imax+1][jmax+1][kmax+1];
double TT[imax+1][jmax+1][kmax+1];

void makeDir(char name[]){
    struct stat buf;
    int ret;
    char dir[256];
    char mkdir[512];
    
    snprintf(dir,256,"%s",name);
    snprintf(mkdir,512,"mkdir %s",dir);
    
    ret=stat(dir, &buf);
    
    if(ret!=0){
        
        ret=system("ls");
        
        if(ret==0){
            
            ret=system(mkdir);
            
            if(ret==0){
                
                printf("\n\n");
                printf("%sフォルダ作成成功! \n ",dir);
                printf("\n\n ");
                
                ret=system("ls");
                
                if(ret!=0){
                    printf("dirコマンド失敗! \n ");
                }
                
            }else{
                printf("%sフォルダ作成失敗! \n ",dir);
            }
            
        }else{
            printf("dirコマンド失敗! \n ");
        }
    }else{
        printf("%sフォルダが存在します \n ",dir);
    }
}

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
    
    //fp = fopen("oneper.txt", "w");			// この文は変数宣言の後に書くのが普通(もしかしたらFILE *fp;の直後でも大丈夫かもしれない)
    // wはwrite、rならread、rwも可
    // このファイルをpov-rayのフォルダに手動コピーしてあげよう
    // プログラム上で日本語名フォルダを指定して移動させることはできないから
    
    Re = 1000;
    Ra = 1000000.0;
    Pr = 5.0;
    dx = 1.0 / (float)imax;
    dy = 1.0 / (float)jmax;
    dz = 1.0 / (float)kmax;
    dt = 0.05 * dx;							// 毎回このdtずつ足していく、その行為がN回行われる、よって最終的なtは(N-1)*dt
    
    N = (int)(10 / dt);						// 今は、結局N=200*imax=20000。ここで、int/(int)floatとしてしまうと分母0でエラー。
    

    //境界条件
    n = 0;
    for(i = 0; i <= imax; i++){
        for(j = 0; j <= jmax; j++){
            for (k = 0; k <= kmax; k++) {
                u[i][j][k] = 0;
                v[i][j][k] = 0;
                w[i][j][k] = 0;
                P[i][j][k] = 0;
                T[i][j][k] = i < j ? 1 : 0;
            }
        }
    }
    
    
    // u[i][j],v[i][j],P[i][j]の計算と書き込み、時刻n*dt(n=1,2,...,N-1)のとき
    // 課題ノートの、[2]~[4]の繰り返し
    for(n = 1; n < N; n++){
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

        //kを入れた場合のfor文
//        for (k=1;k<=kmax-1;k+=kmax-2){
//            for (i=1; i<=imax-1; i++) {
//                for (j=1; j<=; j) {
//
//                }
//            }
//        }
        
        for (i=2; i<= imax-2; i++) {
            for (j=2; j<= jmax-2; j++) {
                for (k=2; k<= kmax-2; k++) {
                    float ududx = uu[i][j][k] * (-uu[i+2][j][k]+8*(uu[i+1][j][k]-uu[i-1][j][k])+uu[i-2][j][k])/(12*dx) + (fabs(uu[i][j][k])*powf(dx, 3)/12.0)*(uu[i+2][j][k]-4*uu[i+1][j][k]+6*uu[i][j][k]-4*uu[i-1][j][k]+uu[i-2][j][k])/powf(dx, 4);//P183の7.75
                    float vdudy = vv[i][j][k] * (-uu[i][j+2][k]+8*(uu[i][j+1][k]-uu[i][j-1][k])+uu[i][j-2][k])/(12*dy) + (fabs(vv[i][j][k])*powf(dy, 3)/12.0)*(uu[i][j+2][k]-4*uu[i][j+1][k]+6*uu[i][j][k]-4*uu[i][j-1][k]+uu[i][j-2][k])/powf(dy, 4);
                    float wdudz = ww[i][j][k] * (-uu[i][j][k+2]+8*(uu[i][j][k-1]-uu[i][j][k-1])+uu[i][j][k-2])/(12*dz) + (fabs(ww[i][j][k])*powf(dz, 3)/12.0)*(uu[i][j][k+2]-4*uu[i][j][k+1]+6*uu[i][j][k]-4*uu[i][j][k-1]+uu[i][j][k-2])/powf(dz, 4);
                    float dPdx = (P[i+1][j][k]-P[i-1][j][k])/(2.0*dx);
                    float dx2 = dx*dx;
                    float dy2 = dy*dy;
                    float dz2 = dz*dz;
                    
                    float d2udx2 = ( uu[i+1][j][k] - 2.0*uu[i][j][k] + uu[i-1][j][k] ) / dx2;
                    float d2udy2 = ( uu[i][j+1][k] - 2.0*uu[i][j][k] + uu[i][j-1][k] ) / dy2;
                    float d2udz2 = ( uu[i][j][k+1] - 2.0*uu[i][j][k] + uu[i][j][k-1] ) / dz2;
                    u[i][j][k] = uu[i][j][k] + dt*( -ududx -vdudy -wdudz - dPdx + 1.0/Re*(d2udx2+d2udy2+d2udz2) );
                    
                    float udvdx = uu[i][j][k] * (-vv[i+2][j][k]+8*(vv[i+1][j][k]-vv[i-1][j][k])+vv[i-2][j][k])/(12*dx) + (fabs(uu[i][j][k])*powf(dx, 3)/12.0)*(vv[i+2][j][k]-4*vv[i+1][j][k]+6*vv[i][j][k]-4*vv[i-1][j][k]+vv[i-2][j][k])/powf(dx, 4);//P183の7.75
                    float vdvdy = vv[i][j][k] * (-vv[i][j+2][k]+8*(vv[i][j+1][k]-vv[i][j-1][k])+vv[i][j-2][k])/(12*dy) + (fabs(vv[i][j][k])*powf(dy, 3)/12.0)*(vv[i][j+2][k]-4*vv[i][j+1][k]+6*vv[i][j][k]-4*vv[i][j-1][k]+vv[i][j-2][k])/powf(dy, 4);//P183の7.75
                    float wdvdz = ww[i][j][k] * (-vv[i][j][k+2]+8*(vv[i][j][k+1]-vv[i][j][k-1])+vv[i][j][k-2])/(12*dz) + (fabs(ww[i][j][k])*powf(dz, 3)/12.0)*(vv[i][j][k+2]-4*vv[i][j][k+1]+6*vv[i][j][k]-4*vv[i][j][k-1]+vv[i][j][k-2])/powf(dz, 4);
                    float dPdy = (P[i][j+1][k]-P[i][j-1][k])/(2.0*dy);
                    float d2vdx2 = ( vv[i+1][j][k] - 2.0*vv[i][j][k] + vv[i-1][j][k] ) / dx2;
                    float d2vdy2 = ( vv[i][j+1][k] - 2.0*vv[i][j][k] + vv[i][j-1][k] ) / dy2;
                    float d2vdz2 = ( vv[i][j][k+1] - 2.0*vv[i][j][k] + vv[i][j][k-1] ) / dz2;
                    v[i][j][k] = vv[i][j][k] + dt*( -udvdx -vdvdy -wdvdz -dPdy + 1.0/Re*(d2vdx2+d2vdy2+d2vdz2) );
                    v[i][j][k] += (Ra/(Re*Re*Pr))*T[i][j][k];  //浮力項はwがないのでyにつける
                    
                    float udwdx = uu[i][j][k] * (-ww[i+2][j][k]+8*(ww[i+1][j][k]-ww[i-1][j][k])+ww[i-2][j][k])/(12*dx) + (fabs(uu[i][j][k])*powf(dx, 3)/12.0)*(ww[i+2][j][k]-4*ww[i+1][j][k]+6*ww[i][j][k]-4*ww[i-1][j][k]+ww[i-2][j][k])/powf(dx, 4);//P183の7.75
                    float vdwdy = vv[i][j][k] * (-ww[i][j+2][k]+8*(ww[i][j+1][k]-ww[i][j-1][k])+ww[i][j-2][k])/(12*dy) + (fabs(vv[i][j][k])*powf(dy, 3)/12.0)*(ww[i][j+2][k]-4*ww[i][j+1][k]+6*ww[i][j][k]-4*ww[i][j-1][k]+ww[i][j-2][k])/powf(dy, 4);//P183の7.75
                    float wdwdz = ww[i][j][k] * (-ww[i][j][k+2]+8*(ww[i][j][k+1]-ww[i][j][k-1])+ww[i][j][k-2])/(12*dz) + (fabs(ww[i][j][k])*powf(dz, 3)/12.0)*(ww[i][j][k+2]-4*ww[i][j][k+1]+6*ww[i][j][k]-4*ww[i][j][k-1]+ww[i][j][k-2])/powf(dz, 4);
                    float dPdz = (P[i][j][k+1]-P[i][j][k-1])/(2.0*dz);
                    float d2wdx2 = ( ww[i+1][j][k] - 2.0*ww[i][j][k] + ww[i-1][j][k] ) / dx2;
                    float d2wdy2 = ( ww[i][j+1][k] - 2.0*ww[i][j][k] + ww[i][j-1][k] ) / dy2;
                    float d2wdz2 = ( ww[i][j][k+1] - 2.0*ww[i][j][k] + ww[i][j][k-1] ) / dz2;
                    w[i][j][k] = ww[i][j][k] + dt*( -udwdx -vdwdy -wdwdz -dPdz + 1.0/Re*(d2wdx2+d2wdy2+d2wdz2) );
                    
                    
                    //上流差分
                    float udTdx  = uu[i][j][k] * (-TT[i+2][j][k]+8*(TT[i+1][j][k]-TT[i-1][j][k])+TT[i-2][j][k])/(12*dx) + (fabs(uu[i][j][k])*powf(dx, 3)/12.0)*(TT[i+2][j][k]-4*TT[i+1][j][k]+6*TT[i][j][k]-4*TT[i-1][j][k]+TT[i-2][j][k])/powf(dt, 4);//P183の7.75;
                    float vdTdy  = vv[i][j][k] * (-TT[i][j+2][k]+8*(TT[i][j+1][k]-TT[i][j-1][k])+TT[i][j-2][k])/(12*dy) + (fabs(vv[i][j][k])*powf(dy, 3)/12.0)*(TT[i][j+2][k]-4*TT[i][j+1][k]+6*TT[i][j][k]-4*TT[i][j-1][k]+TT[i][j-2][k])/powf(dy, 4);//P183の7.75;
                    float wdTdy  = ww[i][j][k] * (-TT[i][j][k+2]+8*(TT[i][j][k+1]-TT[i][j][k-1])+TT[i][j][k-2])/(12*dz) + (fabs(ww[i][j][k])*powf(dz, 3)/12.0)*(TT[i][j][k+2]-4*TT[i][j][k+1]+6*TT[i][j][k]-4*TT[i][j][k-1]+TT[i][j][k-2])/powf(dz, 4);//P183の7.75;
                    float d2Tdx2 = ( TT[i+1][j][k] - 2.0*TT[i][j][k] + TT[i-1][j][k] ) / dx2;
                    float d2Tdy2 = ( TT[i][j+1][k] - 2.0*TT[i][j][k] + TT[i][j-1][k] ) / dy2;
                    float d2Tdz2 = ( TT[i][j][k+1] - 2.0*TT[i][j][k] + TT[i][j][k-1] ) / dz2;
                    
                    T[i][j][k] = TT[i][j][k] + dt * ( -udTdx - vdTdy -wdTdy + 1.0/(Re*Pr)*(d2Tdx2 + d2Tdy2 + d2Tdz2));
                }
            }
        }
        
        for (i=1; i< imax; i++) {
            for (j=1; j< jmax; j+=jmax-2) {//上下の部分を計算
                for (k=1; k< kmax; k++) {
                    u[i][j][k] = uu[i][j][k]
                    + dt * ( - uu[i][j][k] * ( uu[i+1][j][k] - uu[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( uu[i+1][j][k] - 2.0*uu[i][j][k] + uu[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( uu[i][j+1][k] - uu[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( uu[i][j+1][k] - 2.0*uu[i][j][k] + uu[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( uu[i][j][k+1] - uu[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( uu[i][j][k+1] - 2.0*uu[i][j][k] + uu[i][j][k-1] ) / dz
                            - ( P[i+1][j][k] - P[i-1][j][k] ) / (2.0*dx)
                            + 1.0/Re * (  ( uu[i+1][j][k] - 2.0*uu[i][j][k] + uu[i-1][j][k] ) / (dx*dx)
                                        + ( uu[i][j+1][k] - 2.0*uu[i][j][k] + uu[i][j-1][k] ) / (dy*dy)
                                        + ( uu[i][j][k+1] - 2.0*uu[i][j][k] + uu[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    //非線形項の確認
                    //上流差分の式が二行に渡っている部分
                    v[i][j][k] = vv[i][j][k]
                    + dt * ( - uu[i][j][k] * ( vv[i+1][j][k] - vv[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( vv[i+1][j][k] - 2.0*vv[i][j][k] + vv[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( vv[i][j+1][k] - vv[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( vv[i][j+1][k] - 2.0*vv[i][j][k] + vv[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( vv[i][j][k+1] - vv[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( vv[i][j][k+1] - 2.0*vv[i][j][k] + vv[i][j][k-1] ) / dz
                            - ( P[i][j+1][k] - P[i][j-1][k] ) / (2.0*dy)
                            + 1.0/Re * (  ( vv[i+1][j][k] - 2.0*vv[i][j][k] + vv[i-1][j][k] ) / (dx*dx)
                                        + ( vv[i][j+1][k] - 2.0*vv[i][j][k] + vv[i][j-1][k] ) / (dy*dy)
                                        + ( vv[i][j][k+1] - 2.0*vv[i][j][k] + vv[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    
                    w[i][j][k] = ww[i][j][k]
                    + dt * ( - uu[i][j][k] * ( ww[i+1][j][k] - ww[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( ww[i+1][j][k] - 2.0*ww[i][j][k] + ww[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( ww[i][j+1][k] - ww[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( ww[i][j+1][k] - 2.0*ww[i][j][k] + ww[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( ww[i][j][k+1] - ww[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( ww[i][j][k+1] - 2.0*ww[i][j][k] + ww[i][j][k-1] ) / dz
                            - ( P[i][j+1][k] - P[i][j-1][k] ) / (2.0*dy)
                            + 1.0/Re * (  ( ww[i+1][j][k] - 2.0*ww[i][j][k] + ww[i-1][j][k] ) / (dx*dx)
                                        + ( ww[i][j+1][k] - 2.0*ww[i][j][k] + ww[i][j-1][k] ) / (dy*dy)
                                        + ( ww[i][j][k+1] - 2.0*ww[i][j][k] + ww[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    
                    //ここにもTが必要
                    T[i][j][k] = TT[i][j][k]
                    + dt * (- uu[i][j][k] * ( TT[i+1][j][k] - TT[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( TT[i+1][j][k] - 2.0*TT[i][j][k] + TT[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( TT[i][j+1][k] - TT[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( TT[i][j+1][k] - 2.0*TT[i][j][k] + TT[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( TT[i][j][k+1] - TT[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( TT[i][j][k+1] - 2.0*TT[i][j][k] + TT[i][j][k-1] ) / dz
                            + 1.0/(Re*Pr) * (  ( TT[i+1][j][k] - 2.0*TT[i][j][k] + TT[i-1][j][k] ) / (dx*dx)
                                             + ( TT[i][j+1][k] - 2.0*TT[i][j][k] + TT[i][j-1][k] ) / (dy*dy)
                                             + ( TT[i][j][k+1] - 2.0*TT[i][j][k] + TT[i][j][k-1] ) / (dz*dz)
                                             )
                            );
                }
            }
        }
        
        for (i=1;i<imax; i+=imax-2){
            for (j=1; j<jmax; j++) {//横の部分を計算
                for (k=1; k<kmax; k++) {
                    u[i][j][k] = uu[i][j][k]
                    + dt * ( - uu[i][j][k] * ( uu[i+1][j][k] - uu[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( uu[i+1][j][k] - 2.0*uu[i][j][k] + uu[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( uu[i][j+1][k] - uu[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( uu[i][j+1][k] - 2.0*uu[i][j][k] + uu[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( uu[i][j][k+1] - uu[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( uu[i][j][k+1] - 2.0*uu[i][j][k] + uu[i][j][k-1] ) / dz
                            - ( P[i+1][j][k] - P[i-1][j][k] ) / (2.0*dx)
                            + 1.0/Re * (  ( uu[i+1][j][k] - 2.0*uu[i][j][k] + uu[i-1][j][k] ) / (dx*dx)
                                        + ( uu[i][j+1][k] - 2.0*uu[i][j][k] + uu[i][j-1][k] ) / (dy*dy)
                                        + ( uu[i][j][k+1] - 2.0*uu[i][j][k] + uu[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    //非線形項の確認
                    //上流差分の式が二行に渡っている部分
                    v[i][j][k] = vv[i][j][k]
                    + dt * ( - uu[i][j][k] * ( vv[i+1][j][k] - vv[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( vv[i+1][j][k] - 2.0*vv[i][j][k] + vv[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( vv[i][j+1][k] - vv[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( vv[i][j+1][k] - 2.0*vv[i][j][k] + vv[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( vv[i][j][k+1] - vv[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( vv[i][j][k+1] - 2.0*vv[i][j][k] + vv[i][j][k-1] ) / dz
                            - ( P[i][j+1][k] - P[i][j-1][k] ) / (2.0*dy)
                            + 1.0/Re * (  ( vv[i+1][j][k] - 2.0*vv[i][j][k] + vv[i-1][j][k] ) / (dx*dx)
                                        + ( vv[i][j+1][k] - 2.0*vv[i][j][k] + vv[i][j-1][k] ) / (dy*dy)
                                        + ( vv[i][j][k+1] - 2.0*vv[i][j][k] + vv[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    
                    w[i][j][k] = ww[i][j][k]
                    + dt * ( - uu[i][j][k] * ( ww[i+1][j][k] - ww[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( ww[i+1][j][k] - 2.0*ww[i][j][k] + ww[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( ww[i][j+1][k] - ww[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( ww[i][j+1][k] - 2.0*ww[i][j][k] + ww[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( ww[i][j][k+1] - ww[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( ww[i][j][k+1] - 2.0*ww[i][j][k] + ww[i][j][k-1] ) / dz
                            - ( P[i][j+1][k] - P[i][j-1][k] ) / (2.0*dy)
                            + 1.0/Re * (  ( ww[i+1][j][k] - 2.0*ww[i][j][k] + ww[i-1][j][k] ) / (dx*dx)
                                        + ( ww[i][j+1][k] - 2.0*ww[i][j][k] + ww[i][j-1][k] ) / (dy*dy)
                                        + ( ww[i][j][k+1] - 2.0*ww[i][j][k] + ww[i][j][k-1] ) / (dz*dz)
                                        )
                            );
                    
                    //ここにもTが必要
                    T[i][j][k] = TT[i][j][k]
                    + dt * (- uu[i][j][k] * ( TT[i+1][j][k] - TT[i-1][j][k] ) / (2.0*dx)
                            + fabs(uu[i][j][k]) / 2.0 * ( TT[i+1][j][k] - 2.0*TT[i][j][k] + TT[i-1][j][k] ) / dx
                            - vv[i][j][k] * ( TT[i][j+1][k] - TT[i][j-1][k] ) / (2.0*dy)
                            + fabs(vv[i][j][k]) / 2.0 * ( TT[i][j+1][k] - 2.0*TT[i][j][k] + TT[i][j-1][k] ) / dy
                            - ww[i][j][k] * ( TT[i][j][k+1] - TT[i][j][k-1] ) / (2.0*dz)
                            + fabs(ww[i][j][k]) / 2.0 * ( TT[i][j][k+1] - 2.0*TT[i][j][k] + TT[i][j][k-1] ) / dz
                            + 1.0/(Re*Pr) * (  ( TT[i+1][j][k] - 2.0*TT[i][j][k] + TT[i-1][j][k] ) / (dx*dx)
                                             + ( TT[i][j+1][k] - 2.0*TT[i][j][k] + TT[i][j-1][k] ) / (dy*dy)
                                             + ( TT[i][j][k+1] - 2.0*TT[i][j][k] + TT[i][j][k-1] ) / (dz*dz)
                                             )
                            );
                }
            }
        }
        
        for(m = 0; m < 100; m++){
            for(i = 1; i < imax; i++){
                for(j = 1; j < jmax; j++){//ここにもkを足す
                    for(k = 1; k < kmax; k++){//ここにもkを足す
                        A =
                        - 1.0 * ( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx) * ( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx)
                        - 1.0 * ( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy) * ( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy)
                        - 1.0 * ( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz) * ( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz)
                        
                        - 2.0 * ( u[i][j+1][k] - u[i][j-1][k] ) / (2.0*dy)
                              * ( v[i+1][j][k] - v[i-1][j][k] ) / (2.0*dx)
                        - 2.0 * ( u[i][j][k+1] - u[i][j][k-1] ) / (2.0*dz)
                              * ( w[i+1][j][k] - w[i-1][j][k] ) / (2.0*dx)
                        - 2.0 * ( v[i][j][k+1] - v[i][j][k-1] ) / (2.0*dz)
                              * ( w[i][j+1][k] - w[i][j-1][k] ) / (2.0*dy)
                        + 1.0/dt * (
                                    ( u[i+1][j][k] - u[i-1][j][k] ) / (2.0*dx) +
                                    ( v[i][j+1][k] - v[i][j-1][k] ) / (2.0*dy) +
                                    ( w[i][j][k+1] - w[i][j][k-1] ) / (2.0*dz)
                                    );//TODO
                        
                        B = - ( P[i+1][j][k] + P[i-1][j][k] ) / (dx*dx)
                            - ( P[i][j+1][k] + P[i][j-1][k] ) / (dy*dy)
                            - ( P[i][j][k+1] + P[i][j][k-1] ) / (dz*dz) + A;
                        P[i][j][k] = B / ( - 2.0 / (dx*dx) - 2.0 / (dy*dy) - 2.0 / (dz*dz) );
                    }
                }
            }
        }
        
        // 境界条件
        for(j = 0; j <= jmax; j++){
            for(k = 0; k <= kmax; k++){
                P[0][j][k] = P[1][j][k];
                P[imax][j][k] = P[imax-1][j][k];
            }
        }
        for(i = 0; i <= imax; i++){
            for(k = 0; k <= kmax; k++){
                P[i][0][k] = P[i][1][k];
                P[i][jmax][k] = P[i][jmax-1][k];
            }
        }
        for(i = 0; i <= imax; i++){
            for(j = 0; j <= jmax; j++){
                P[i][j][0] = P[i][j][1];
                P[i][j][kmax] = P[i][j][kmax-1];
            }
        }
        
        char filename[100] = {'\0'};
        snprintf(filename, 100, "%s/paraview.csv.%d",dirName,n);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, "t,x,y,z,u,v,w,T,\n");
        // 時刻n*dt(第n回目時点)の書き込み
        int flag = 0;
        for(i = 0; i <= imax; i++){
            for(j = 0; j <= jmax; j++){
                for (k = 0; k <= kmax; k++) {
                    if (flag == 0){
                        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,\n", n*dt,i*dx,j*dy,k*dz,u[i][j][k],v[i][j][k],w[i][j][k], T[i][j][k]);
                        flag = 1;
                    }if (flag == 1){
                        flag = 2;
                    }else if (flag == 2){
                        flag = 0;
                    }
                }
            }
        }
        fclose(fp);
        // 計算が長くなるため、コマンドプロンプトに進捗状況を表示させる
        // 今回は全部で20000回
        printf("%d/10000\n", n+1);
        // printf("終了");
    }
    printf("終了");
}

void mydraw(){
    //枠
//    myline(0, 0, 0, 1, 1, 1, 1, 1);
//    myline(0, 1, 1, 1, 1, 1, 1, 1);
//    myline(1, 1, 1, 0, 1, 1, 1, 1);
//    myline(1, 0, 0, 0, 1, 1, 1, 1);
    calcData();
//    showData(t);
}

//
//void mydraw_xx(){
//    double u[N+1], uu[N+1];				// N分割するから点はN+1個。uはn+1、uuはnのときのもの。
//    float dx, dt;						// dx, dtをdoubleにしてしまうと、n>=1のときのu[j]の値がなぜかおかしくなる。
//    int j, n;
//    dx = 1.0 / (float)N / 2.0;				// 「x」は0～1。
//    dt = dx * dx / 3.0 / 2.0;				// ρ=Δt/(Δx)^2=1/6にしたい
//
//    // 初期条件(すなわちn=0,t=0のとき。ここより、最初両端ではuは0になる。これをいじらなければずっと両端は0で出力される)
//    for(j = 0; j <= N; j++){
//        if(0 <= j && j <= N/2){
//            u[j] = j * dx;			// 半分より左側
//        }
//        if(N/2 < j && j <= N ){
//            u[j] = 1.0 - j * dx;		// 半分より右側
//        }
//        mypoint(j*dx, u[j], 1, 0, 0, 3);
//        myline((j-1)*dx, u[j-1], (j)*dx, u[j], 1, 0, 0, 1);
//    }
//    
//    // ここ以降で出力されるu[j]の値が、dx,dtがdoubleだとおかしくなってしまう。
//    for(n = 1; n <= 1.0/dt; n++){		// dtは0未満の小数だから、1.0/dtは1より大きい
//        for(j = 1; j < N; j++){
//            uu[j] = u[j];				// u[j]の値をコピーしておく
//        }
//        for(j = 1; j < N; j++){			// 両端j=0,Nのときuは初期値0、これらはいじらなければ常に0
//            u[j] = uu[j] + dt/(dx*dx) * (uu[j+1] - 2*uu[j] + uu[j-1]);		// dt/(dx*dx)は常に1.0/6.0だから、そう書いてもOK
//            printf("%f\n",dt/(dx*dx));
//            //u[j] = uu[j] + RHO * (uu[j+1] - 2*uu[j] + uu[j-1]);		// dt/(dx*dx)は常に1.0/6.0だから、そう書いてもOK
//        }
//        if(n%10 == 0){					// 出力される図がベタ塗りにならないように間引く
//            for(j = 0; j <= N; j++){		// 両端j=0,Nのときも出力するため、出力のためのforを別個に作る
//                printf("%f\n",u[j]);
//                mypoint(j*dx, u[j], 1, 0, 0, 3);
////                myline(j*dx, u[j], (j+1)*dx, u[j+1], 1, 0, 0, 1);
//                myline((j-1)*dx, u[j-1], (j)*dx, u[j], 1, 0, 0, 1);
//            }
//        }
//    }
//}

/**draw関数*/
void showData(int fileNo){
    char filename[10] = {'\0'};
    snprintf(filename, 10, "%.d.txt", fileNo);
    FILE *file;
    if((file=fopen(filename,"r"))!=NULL){
        double buf[6];
        //ファイルが終わるまで読み込む
        int index = 0;
        double dt = 0.2;
        while( fscanf(file,"%lf,%lf,%lf,%lf,%lf,%lf,",&buf[0],&buf[1],&buf[2],&buf[3],&buf[4],&buf[5]) != EOF ){
            index++;
            //if (index%4==0){//間引き
            //myline(buf[1],buf[2], buf[1]+(dt*buf[3]),buf[2]+(dt*buf[4]),  1, 0, 0, 1);
            powLine(buf[1],//x0
                    buf[2],//y0
                    buf[1]+(dt*buf[3]),//x1
                    buf[2]+(dt*buf[4]));//y1
            //}
            for(int i=0;i<10 ;i++){
                char disp[12] = {'\0'};
                snprintf(disp, 12, "t = %.03f", buf[0]);
                glRasterPos3f(0.05*i - 0.5,-0.6,0);
                glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, disp[i]);
            }
        }
        fclose(file);
    }
}
