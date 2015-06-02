#include <math.h>
#include "myglini.h"

/* x3d[],y3d[],z3d[] : 物体の拡張展の３次元座標*/
/*f[][]:               物体を構成する各面と頂点との関係*/
/*nv[]:                物体の各面に含まれる頂点の数(<=10)*/
/*n_vertex:             物体に含まれる頂点の数(<=20)*/
/*n_face:               物体に含まれる面の数(<=12)*/

float x3d[20],y3d[20],z3d[20];
int f[12][10],nv[12],n_vertex,n_face;

void modeling (void);/*モデリング*/
void calc_modeling (int n);/*モデリング(計算で求める)*/

void mydraw(){
    float x[20],y[20],z0 = 1.0f; /*投影面z0 = 1.0*/
    float R,G,B;
    int i,k,i0;
    /*物体の位相と幾何の定義（モデリング）*/
    //modeling();
    calc_modeling(10);
    /*物体の色*/
    R = 1.0f;
    G = 1.0f;
    B = 1.0f;
    /*物体の描画*/
    for (k=0; k<n_face; k++) {
        /*を引く（透視投影）*/
        for (i=0; i<nv[k]; i=i+1) {/*nv[k]:面F_kに含まれる頂点数*/
            i0 = f[k][i];/*面f_kのi番目の頂点番号をi0に代入*/
            x[i] = z0*x3d[i0]/z3d[i0];/*透視投影による投影変換*/
            y[i] = z0*y3d[i0]/z3d[i0];/*透視投影による投影変換*/
        }
        mylineloop(nv[k], x, y, R, G, B, 1.0);/*(nv[k])角形の描画*/
    }
}

//calc_modelingは与えられたn角形のモデルデータ（表）を生成する
void calc_modeling (int n){
    
    float pi;
    /*πの値をpiに代入*/
    pi = 2.0 * acos(0.0);
    
    /*n角形錐台（n>=6）の場合について各頂点の座標を指定する*/
    n_vertex = 2 * n;
    for (int i=0; i<n_vertex; i++) {
        if (i<n) {
            x3d[i] = cos(i * 2.0f * pi / n);
            y3d[i] = -1.0f;
            z3d[i] = sinf(i * 2.0f * pi / n) + 3.0f;
        }else{
            x3d[i] = 0.5f * cos((i - n) * 2.0f * pi / n);
            y3d[i] = 1.0f;
            z3d[i] = 0.5f * sin((i - n) * 2.0f * pi / n) + 3.0f;
        }
    }
    
    /*物体に含まれる各面f[k]と頂点の関係*/
    /*各面f[k]の頂点f[k][i]を物体の外側から見て反時計まわりにf[k][0,1,2,...]で順序付けて定義する*/
    /*先ほど作成したfとnvの表をそのまま書けば良い*/
    n_face = n + 2;
    
    for (int i=0; i<n_face; i++) {
        for (int j=0; j<n; j++) {
            if (i==0) {//底面
                f[i][j] = j;
            }else if (i==n_face-1){//上面
                f[i][j] = n + j;
            }else{//側面
                int tmp;
                if (j==0) {//左下
                    tmp = i-1;//i-1を基準とする
                }else if (j==1){
                    tmp = i;//右下
                    if (tmp>n-1) {//底面に接するのでn-1より大きくなることはない
                        tmp -= n;
                    }
                }else if (j==3){//左上
                    //j=0の真上になる
                    tmp = i-1+n;
                }else if (j==2){//右上
                    tmp = i;//右下
                    if (tmp>n-1) {//底面に接するのでn-1より大きくなることはない
                        tmp -= n;
                    }
                    tmp += n;
                }
                f[i][j] = tmp;
            }
        }
        
        if (i==0) {//底面
            nv[i] = n;
        }else if (i==n_face-1){//上面
            nv[i] = n;
        }else{//側面
            nv[i] = 4;
        }
    }
}

void modeling (void){
    float pi;
    /*πの値をpiに代入*/
    pi = 2.0 * acos(0.0);
    
    /*n角形錐台（n>=6）の場合について各頂点の座標を指定する*/
    int n;
    n = 5;
    n_vertex = 2 * n;
    
    x3d[0] = cos(0.0f * 2.0f * pi / 5.0f);
    x3d[1] = cos(1.0f * 2.0f * pi / 5.0f);
    x3d[2] = cos(2.0f * 2.0f * pi / 5.0f);
    x3d[3] = cos(3.0f * 2.0f * pi / 5.0f);
    x3d[4] = cos(4.0f * 2.0f * pi / 5.0f);
    x3d[5] = 0.5f * cos(0.0f * 2.0f * pi / 5.0f);
    x3d[6] = 0.5f * cos(1.0f * 2.0f * pi / 5.0f);
    x3d[7] = 0.5f * cos(2.0f * 2.0f * pi / 5.0f);
    x3d[8] = 0.5f * cos(3.0f * 2.0f * pi / 5.0f);
    x3d[9] = 0.5f * cos(4.0f * 2.0f * pi / 5.0f);
    
    y3d[0] = -1.0f;
    y3d[1] = -1.0f;
    y3d[2] = -1.0f;
    y3d[3] = -1.0f;
    y3d[4] = -1.0f;
    y3d[5] = 1.0f;
    y3d[6] = 1.0f;
    y3d[7] = 1.0f;
    y3d[8] = 1.0f;
    y3d[9] = 1.0f;
    
    z3d[0] = sin(0.0 * 2 * pi / 5) + 3.0f;
    z3d[1] = sin(1.0 * 2 * pi / 5) + 3.0f;
    z3d[2] = sin(2.0 * 2 * pi / 5) + 3.0f;
    z3d[3] = sin(3.0 * 2 * pi / 5) + 3.0f;
    z3d[4] = sin(4.0 * 2 * pi / 5) + 3.0f;
    z3d[5] = 0.5f * sin(0.0 * 2 * pi / 5) + 3.0f;
    z3d[6] = 0.5f * sin(1.0 * 2 * pi / 5) + 3.0f;
    z3d[7] = 0.5f * sin(2.0 * 2 * pi / 5) + 3.0f;
    z3d[8] = 0.5f * sin(3.0 * 2 * pi / 5) + 3.0f;
    z3d[9] = 0.5f * sin(4.0 * 2 * pi / 5) + 3.0f;
    
    nv[0] = 5;
    nv[1] = 4;
    nv[2] = 4;
    nv[3] = 4;
    nv[4] = 4;
    nv[5] = 4;
    nv[6] = 5;
    
    /*物体に含まれる各面f[k]と頂点の関係*/
    /*各面f[k]の頂点f[k][i]を物体の外側から見て反時計まわりにf[k][0,1,2,...]で順序付けて定義する*/
    /*先ほど作成したfとnvの表をそのまま書けば良い*/
    n_face = n + 2;

    f[0][0] = 0;f[0][1] = 1;f[0][2] = 2;f[0][3] = 3;f[0][4] = 4;
    f[1][0] = 0;f[1][1] = 1;f[1][2] = 6;f[1][3] = 5;
    f[2][0] = 1;f[2][1] = 2;f[2][2] = 7;f[2][3] = 6;
    f[3][0] = 2;f[3][1] = 3;f[3][2] = 8;f[3][3] = 7;
    f[4][0] = 3;f[4][1] = 4;f[4][2] = 9;f[4][3] = 8;
    f[5][0] = 4;f[5][1] = 0;f[5][2] = 5;f[5][3] = 9;
    f[6][0] = 5;f[6][1] = 6;f[6][2] = 7;f[6][3] = 8;f[6][4] = 9;
}

void mydraw_002(){
//void mydraw_0428(){
    float xmin,xmax,ymin,ymax,x[4],y[4];
    float R0,G0,B0,R1,G1,B1;
    float R,G,B,h;
    int n;
    
    //(186, 85, 211)
    R0 = 0.9;
    G0 = 0.9;
    B0 = 0;
    
    //852BA9
    R1 = 0.6;
    G1 = 0;
    B1 = 0.6;
    
    /*長方形の数*/
    n = 256;
    /*xの最小値*/
    xmin = -1.0;
    /*xの最大値*/
    xmax = 1.0;
    /*yの最小値*/
    ymin = -1.0;
    /*yの最大値*/
    ymax = 1.0;
    /*長方形の横幅*/
    h = (ymax - ymin) / n;
    
    for (int i=0; i<n; i++) {
        /*i番目の長方形の４頂点の座標*/
        x[0] = xmin;
        y[0] = ymin + h * i;
        x[1] = xmax;
        y[1] = y[0];
        x[2] = x[1];
        y[2] = y[1] + h;
        x[3] = x[0];
        y[3] = y[2];
        
        /*i番目の長方形の色*/
        R = R0 + (float)i*(R1 - R0)/(n-1);
        G = G0 + (float)i*(G1 - G0)/(n-1);
        B = B0 + (float)i*(B1 - B0)/(n-1);
        
        /*i番目の長方形の描画*/
        mypolygon(4, x, y, R, G, B);
    }
}

void mydraw_001(){//課題１日目
    for (int i=0; i<100; i++) {
        mypoint(i/100.0, 0, 1, 1, 1, 2);
        mypoint(0, i/100.0, 1, 1, 1, 2);
    }
}

void mydraw_sample()
{
    
    float x0,y0,x1,y1,the,dthe,pi;
    float x[100],y[100];
    
    pi=2.0*acos(0.0);
    dthe=pi/20.0;
    
    /* Line */
    myline(-1.0,-1.0,1.0,1.0, 0.5,1.0,0.0, 2.0);
    
    /* Square */
    x[0]=-0.5; y[0]=-0.5;
    x[1]=-0.5; y[1]= 0.5;
    x[2]= 0.5; y[2]= 0.5;
    x[3]= 0.5; y[3]=-0.5;
    mylineloop(4, x,y, 1.0,1.0,1.0, 1.0);
    
    /* Circle */
    for(the=0.0; the<= 2.0*pi; the=the+dthe){
        x0=0.5*cos(the);      y0=0.5*sin(the);
        x1=0.5*cos(the+dthe); y1=0.5*sin(the+dthe);
        myline(x0,y0,x1,y1, 1.0,0.0,0.0, 5.0);
    }
    
    /* Triangle */
    x[0]= 0.1; y[0]= 0.1;
    x[1]= 0.0; y[1]=-0.1;
    x[2]=-0.1; y[2]= 0.1;
    mypolygon(3, x,y, 0.8,0.0,1.0);
    
    /* Point */
    mypoint(0.0,0.0, 0.0,0.0,1.0, 3.0);
    
    
}
