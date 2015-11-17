
//#include <GL/glut.h>
//GLの互換性のあるGLUTを利用する
#include <GLUT/glut.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#pragma clang diagnostic ignored "-Wdeprecated-declarations"//deplicateワーニングを抑制

void mydraw(void);
void showData(int);
GLint t = 0;

void mypoint(float x, float y, float r, float g, float b, float s)
{
    x-=1;
    y-=1;
    if(r > 1.0) r=1.0; if(g > 1.0) g=1.0; if(b > 1.0) b=1.0;
    if(r < 0.0) r=0.0; if(g < 0.0) g=0.0; if(b < 0.0) b=0.0;
    glColor3f(r,g,b);
    glPointSize(s);
    glBegin(GL_POINTS);
    glVertex2f(x,y);
    glEnd();
}

void myline(float x0, float y0, float x1, float y1,
            float r, float g, float b, float s)
{
    x0-=1;
    y0-=1;
    x1-=1;
    y1-=1;
    if(r > 1.0) r=1.0; if(g > 1.0) g=1.0; if(b > 1.0) b=1.0;
    if(r < 0.0) r=0.0; if(g < 0.0) g=0.0; if(b < 0.0) b=0.0;
    glColor3f(r,g,b);
    glLineWidth(s);
    glBegin(GL_LINES);
    glVertex2f(x0,y0);
    glVertex2f(x1,y1);
    glEnd();
}

//
void powLine(float x0, float y0, float x1, float y1){

    float width = x1 - x0;
    float height = y1 - y0;
    float  value = fabs(sqrt(width*width + height*height))*20.0;

    int r, g, b;    // RGB値
    double  tmp_val = cos( 4 * M_PI * value );
    int     col_val = (int)( ( -tmp_val / 2 + 0.5 ) * 255 );
    if ( value >= ( 4.0 / 4.0 ) ) { r = 255;     g = 0;       b = 0;       }   // 赤
    else if ( value >= ( 3.0 / 4.0 ) ) { r = 255;     g = col_val; b = 0;       }   // 黄～赤
    else if ( value >= ( 2.0 / 4.0 ) ) { r = col_val; g = 255;     b = 0;       }   // 緑～黄
    else if ( value >= ( 1.0 / 4.0 ) ) { r = 0;       g = 255;     b = col_val; }   // 水～緑
    else if ( value >= ( 0.0 / 4.0 ) ) { r = 0;       g = col_val; b = 255;     }   // 青～水
    else {                               r = 0;       g = 0;       b = 255;     }   // 青
    
    
    int isLine = 1;
    if (isLine) {
        myline(x0*2.0,
               y0*2.0,
               x1*2.0,
               y1*2.0,
               r/255.0, g/255.0, b/255.0, 1.7);
        mypoint(x1*2.0, y1*2.0, r/255.0, g/255.0, b/255.0, 5.0);
    }else{
//        mypoint(x0*2.0, y0*2.0, value, 0, 1.0-value, 20.0);
        mypoint(x0*2.0, y0*2.0, r/255.0, g/255.0, b/255.0, 20.0);
    }
}

void mypolygon(int n, float *x0, float *y0, float r, float g, float b)
{
    int i;
    float x,y;
    
    if(r > 1.0) r=1.0; if(g > 1.0) g=1.0; if(b > 1.0) b=1.0;
    if(r < 0.0) r=0.0; if(g < 0.0) g=0.0; if(b < 0.0) b=0.0;
    glColor3f(r,g,b);
    glBegin(GL_POLYGON);
    for(i=0;i<n;++i) {
        x=*(x0+i); y=*(y0+i);
        glVertex2f(x,y);
    }
    glEnd();
}

void mylineloop(int n, float *x0, float *y0,
                float r, float g, float b, float s)
{
    int i;
    float x,y;
    
    if(r > 1.0) r=1.0; if(g > 1.0) g=1.0; if(b > 1.0) b=1.0;
    if(r < 0.0) r=0.0; if(g < 0.0) g=0.0; if(b < 0.0) b=0.0;
    glColor3f(r,g,b);
    glLineWidth(s);
    glBegin(GL_LINE_LOOP);
    for(i=0;i<n;++i) {
        x=*(x0+i); y=*(y0+i);
        glVertex2f(x,y);
    }
    glEnd();
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT);
    mydraw();
    glFlush();
}

void init(){
    glClearColor(0.0, 0.0, 0.0, 0.0);
}

/* 100ミリ秒ごとに実行される関数 */
void timer(int value) {
    /* 正方形のサイズを増加 */
    t += 50;
    if (t>5000){
        t = 0;
    }
    //showData(t);
        
    /* 画面を再描写 */
    glutPostRedisplay();
    /* 100ミリ秒後に再実行 */
    glutTimerFunc(100, timer, 0);
}

int main(int argc, char** argv)
{
    mydraw();
    //    glutInit(&argc, argv);
//    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
//    glutInitWindowSize(800, 800);
//    glutInitWindowPosition(0, 100);
//    glutCreateWindow("[-1,1] X [-1,1]");
//    glutDisplayFunc(display);
//    
//    //glutTimerFunc(50, timer, 0);
//    init();
//    glutMainLoop();
    
}
