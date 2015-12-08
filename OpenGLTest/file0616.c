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

#define imax 50
#define jmax 50

double u[imax+1][jmax+1];
double uu[imax+1][jmax+1];
double v[imax+1][jmax+1];
double vv[imax+1][jmax+1];

double P[imax+1][jmax+1];
double T[imax+1][jmax+1];
double TT[imax+1][jmax+1];

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


void calcData(){
    //レイノルズ数が大きいほど物質はさらさらする
    //風邪が吹くと物体の反対側で渦が生じ
    float Re, Ra, Pr;
    int i, j, n, m, N;					// 型の宣言時には、変数演算を用いて値を代入しながら宣言することはできない
    float dx, dy, dt;						// したがってdxを用意した後dtを用意するとき、float dt=0.05*dx;とせず、float dt; dt=0.05*dx;とする
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
    dt = 0.05 * dx;							// 毎回このdtずつ足していく、その行為がN回行われる、よって最終的なtは(N-1)*dt
    
    N = (int)(10 / dt);						// 今は、結局N=200*imax=20000。ここで、int/(int)floatとしてしまうと分母0でエラー。
    

    //境界条件
    n = 0;
    for(i = 0; i <= imax; i++){
        for(j = 0; j <= jmax; j++){
            v[i][j] = 0;
            P[i][j] = 0;
            T[i][j] = i < j ? 1 : 0;
        }
    }
    
    
    // u[i][j],v[i][j],P[i][j]の計算と書き込み、時刻n*dt(n=1,2,...,N-1)のとき
    // 課題ノートの、[2]~[4]の繰り返し
    for(n = 1; n < N; n++){
        for(i = 0; i <= imax; i++){
            for(j = 0; j <= jmax; j++){
                uu[i][j] = u[i][j];			// 全ての[i][j]で、uをコピーしておく
                vv[i][j] = v[i][j];
                TT[i][j] = T[i][j];
            }
        }

        //kを入れた場合のfor文
//        for (k=1;k<=kmax-1;k+=kmax-2){
//            for (i=1; i<=imax-1; i++) {
//                for (j=1; j<=; j) {
//                    <#statements#>
//                }
//            }
//        }
        
        for (i=2; i<= imax-2; i++) {
            for (j=2; j<= jmax-2; j++) {
                float ududx = uu[i][j] * (-uu[i+2][j]+8*(uu[i+1][j]-uu[i-1][j])+uu[i-2][j])/(12*dx) + (fabs(uu[i][j])*powf(dx, 3)/12.0)*(uu[i+2][j]-4*uu[i+1][j]+6*uu[i][j]-4*uu[i-1][j]+uu[i-2][j])/powf(dx, 4);//P183の7.75
                float vdudy = vv[i][j] * (-uu[i][j+2]+8*(uu[i][j+1]-uu[i][j-1])+uu[i][j-2])/(12*dy) + (fabs(vv[i][j])*powf(dy, 3)/12.0)*(uu[i][j+2]-4*uu[i][j+1]+6*uu[i][j]-4*uu[i][j-1]+uu[i][j-2])/powf(dy, 4);
                float dPdx = (P[i+1][j]-P[i-1][j])/(2.0*dx);
                float dx2 = dx*dx;
                float dy2 = dy*dy;
                float d2udx2 = ( uu[i+1][j] - 2.0*uu[i][j] + uu[i-1][j] ) / dx2;
                float d2udy2 = ( uu[i][j+1] - 2.0*uu[i][j] + uu[i][j-1] ) / dy2;
                u[i][j] = uu[i][j] + dt*( -ududx -vdudy - dPdx + 1.0/Re*(d2udx2+d2udy2) );
                
                float udvdx = uu[i][j] * (-vv[i+2][j]+8*(vv[i+1][j]-vv[i-1][j])+vv[i-2][j])/(12*dx) + (fabs(uu[i][j])*powf(dx, 3)/12.0)*(vv[i+2][j]-4*vv[i+1][j]+6*vv[i][j]-4*vv[i-1][j]+vv[i-2][j])/powf(dx, 4);//P183の7.75
                float vdvdy = vv[i][j] * (-vv[i][j+2]+8*(vv[i][j+1]-vv[i][j-1])+vv[i][j-2])/(12*dy) + (fabs(vv[i][j])*powf(dy, 3)/12.0)*(vv[i][j+2]-4*vv[i][j+1]+6*vv[i][j]-4*vv[i][j-1]+vv[i][j-2])/powf(dy, 4);//P183の7.75
                float dPdy = (P[i][j+1]-P[i][j-1])/(2.0*dy);
                float d2vdx2 = ( vv[i+1][j] - 2.0*vv[i][j] + vv[i-1][j] ) / dx2;
                float d2vdy2 = ( vv[i][j+1] - 2.0*vv[i][j] + vv[i][j-1] ) / dy2;
                v[i][j] = vv[i][j] + dt*( -udvdx -vdvdy -dPdy + 1.0/Re*(d2vdx2+d2vdy2) );
                v[i][j] += (Ra/(Re*Re*Pr))*T[i][j];  //浮力項はwがないのでyにつける
                
                //上流差分
                float udTdx  = uu[i][j] * (-TT[i+2][j]+8*(TT[i+1][j]-TT[i-1][j])+TT[i-2][j])/(12*dx) + (fabs(uu[i][j])*powf(dx, 3)/12.0)*(TT[i+2][j]-4*TT[i+1][j]+6*TT[i][j]-4*TT[i-1][j]+TT[i-2][j])/powf(dt, 4);//P183の7.75;
                float vdTdy  = vv[i][j] * (-TT[i][j+2]+8*(TT[i][j+1]-TT[i][j-1])+TT[i][j-2])/(12*dy) + (fabs(vv[i][j])*powf(dy, 3)/12.0)*(TT[i][j+2]-4*TT[i][j+1]+6*TT[i][j]-4*TT[i][j-1]+TT[i][j-2])/powf(dy, 4);//P183の7.75;
                float d2Tdx2 = ( TT[i+1][j] - 2.0*TT[i][j] + TT[i-1][j] ) / dx2;
                float d2Tdy2 = ( TT[i][j+1] - 2.0*TT[i][j] + TT[i][j-1] ) / dy2;
                
                T[i][j] = TT[i][j] + dt * ( -udTdx - vdTdy + 1.0/(Re*Pr)*(d2Tdx2 + d2Tdy2));
            }
        }
        
        for (i=1; i< imax; i++) {
            for (j=1; j< jmax; j+=jmax-2) {//上下の部分を計算
                u[i][j] = uu[i][j]
                + dt * ( - uu[i][j] * ( uu[i+1][j] - uu[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( uu[i+1][j] - 2.0*uu[i][j] + uu[i-1][j] ) / dx
                        - vv[i][j] * ( uu[i][j+1] - uu[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( uu[i][j+1] - 2.0*uu[i][j] + uu[i][j-1] ) / dy
                        - ( P[i+1][j] - P[i-1][j] ) / (2.0*dx)
                        + 1.0/Re * ( ( uu[i+1][j] - 2.0*uu[i][j] + uu[i-1][j] ) / (dx*dx)
                                    + ( uu[i][j+1] - 2.0*uu[i][j] + uu[i][j-1] ) / (dy*dy)
                                    )
                        );
                //非線形項の確認
                //上流差分の式が二行に渡っている部分
                v[i][j] = vv[i][j]
                + dt * ( - uu[i][j] * ( vv[i+1][j] - vv[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( vv[i+1][j] - 2.0*vv[i][j] + vv[i-1][j] ) / dx
                        - vv[i][j] * ( vv[i][j+1] - vv[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( vv[i][j+1] - 2.0*vv[i][j] + vv[i][j-1] ) / dy
                        - ( P[i][j+1] - P[i][j-1] ) / (2.0*dy)
                        + 1.0/Re * ( ( vv[i+1][j] - 2.0*vv[i][j] + vv[i-1][j] ) / (dx*dx)
                                    + ( vv[i][j+1] - 2.0*vv[i][j] + vv[i][j-1] ) / (dy*dy)
                                    )
                        );
                
                //ここにもTが必要
                T[i][j] = TT[i][j]
                + dt * ( - uu[i][j] * ( TT[i+1][j] - TT[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( TT[i+1][j] - 2.0*TT[i][j] + TT[i-1][j] ) / dx
                        - vv[i][j] * ( TT[i][j+1] - TT[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( TT[i][j+1] - 2.0*TT[i][j] + TT[i][j-1] ) / dy
                        + 1.0/(Re*Pr) * ( ( TT[i+1][j] - 2.0*TT[i][j] + TT[i-1][j] ) / (dx*dx)
                                    + ( TT[i][j+1] - 2.0*TT[i][j] + TT[i][j-1] ) / (dy*dy)
                                    )
                        );
            }
        }
        
        for (i=1;i<imax; i+=imax-2){
            for (j=1; j<jmax; j++) {//横の部分を計算
                u[i][j] = uu[i][j]
                + dt * ( - uu[i][j] * ( uu[i+1][j] - uu[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( uu[i+1][j] - 2.0*uu[i][j] + uu[i-1][j] ) / dx
                        - vv[i][j] * ( uu[i][j+1] - uu[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( uu[i][j+1] - 2.0*uu[i][j] + uu[i][j-1] ) / dy
                        - ( P[i+1][j] - P[i-1][j] ) / (2.0*dx)
                        + 1.0/Re * ( ( uu[i+1][j] - 2.0*uu[i][j] + uu[i-1][j] ) / (dx*dx)
                                    + ( uu[i][j+1] - 2.0*uu[i][j] + uu[i][j-1] ) / (dy*dy)
                                    )
                        );
                v[i][j] = vv[i][j]
                + dt * ( - uu[i][j] * ( vv[i+1][j] - vv[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( vv[i+1][j] - 2.0*vv[i][j] + vv[i-1][j] ) / dx
                        - vv[i][j] * ( vv[i][j+1] - vv[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( vv[i][j+1] - 2.0*vv[i][j] + vv[i][j-1] ) / dy
                        - ( P[i][j+1] - P[i][j-1] ) / (2.0*dy)
                        + 1.0/Re * ( ( vv[i+1][j] - 2.0*vv[i][j] + vv[i-1][j] ) / (dx*dx)
                                    + ( vv[i][j+1] - 2.0*vv[i][j] + vv[i][j-1] ) / (dy*dy)
                                    )
                        );
                //ここにもTが必要
                
                //ここにもTが必要
                T[i][j] = TT[i][j]
                + dt * ( - uu[i][j] * ( TT[i+1][j] - TT[i-1][j] ) / (2.0*dx)
                        + fabs(uu[i][j]) / 2.0 * ( TT[i+1][j] - 2.0*TT[i][j] + TT[i-1][j] ) / dx
                        - vv[i][j] * ( TT[i][j+1] - TT[i][j-1] ) / (2.0*dy)
                        + fabs(vv[i][j]) / 2.0 * ( TT[i][j+1] - 2.0*TT[i][j] + TT[i][j-1] ) / dy
                        + 1.0/(Re*Pr) * ( ( TT[i+1][j] - 2.0*TT[i][j] + TT[i-1][j] ) / (dx*dx)
                                         + ( TT[i][j+1] - 2.0*TT[i][j] + TT[i][j-1] ) / (dy*dy)
                                         )
                        );
            }
        }
        
        for(m = 0; m < 100; m++){
            for(i = 1; i < imax; i++){
                for(j = 1; j < jmax; j++){//ここにもkを足す
                    A = - 1.0 * ( u[i+1][j] - u[i-1][j] ) / (2.0*dx) * ( u[i+1][j] - u[i-1][j] ) / (2.0*dx)
                    - 1.0 * ( v[i][j+1] - v[i][j-1] ) / (2.0*dy) * ( v[i][j+1] - v[i][j-1] ) / (2.0*dy)
                    - 2.0 * ( u[i][j+1] - u[i][j-1] ) / (2.0*dy) * ( v[i+1][j] - v[i-1][j] ) / (2.0*dx)
                    + 1.0/dt * ( ( u[i+1][j] - u[i-1][j] ) / (2.0*dx) + ( v[i][j+1] - v[i][j-1] ) / (2.0*dy) );
                    B = - ( P[i+1][j] + P[i-1][j] ) / (dx*dx) - ( P[i][j+1] + P[i][j-1] ) / (dy*dy) + A;
                    P[i][j] = B / ( - 2.0 / (dx*dx) - 2.0 / (dy*dy) );
                }
            }
        }
        
        // 課題ノートの[4]
        for(j = 0; j <= jmax; j++){
            P[0][j] = P[1][j];
            P[imax][j] = P[imax-1][j];
        }
        for(i = 0; i <= imax; i++){
            P[i][0] = P[i][1];
            P[i][jmax] = P[i][jmax-1];
        }
        
        char filename[20] = {'\0'};
        snprintf(filename, 20, "paraview_%d.csv",n);
        FILE *fp = fopen(filename, "w");
        fprintf(fp, "t,x,y,u,v,P,\n");
        // 時刻n*dt(第n回目時点)の書き込み
        int flag = 0;
        for(i = 0; i <= imax; i++){
            for(j = 0; j <= jmax; j++){
                if (flag == 0){
                    fprintf(fp, "%f,%f,%f,%f,%f,%f,\n", n*dt,i*dx, j*dy, u[i][j], v[i][j], P[i][j]);
                    flag = 1;
                }if (flag == 1){
                    flag = 2;
                }else if (flag == 2){
                    flag = 0;
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