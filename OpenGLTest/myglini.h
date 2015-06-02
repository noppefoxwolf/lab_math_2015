
//#include <GL/glut.h>
//GLの互換性のあるGLUTを利用する
#include <GLUT/glut.h>
#pragma clang diagnostic ignored "-Wdeprecated-declarations"//deplicateワーニングを抑制

void mydraw(void);

void mypoint(float x, float y, float r, float g, float b, float s)
{
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
    if(r > 1.0) r=1.0; if(g > 1.0) g=1.0; if(b > 1.0) b=1.0;
    if(r < 0.0) r=0.0; if(g < 0.0) g=0.0; if(b < 0.0) b=0.0;
    glColor3f(r,g,b);
    glLineWidth(s);
    glBegin(GL_LINES);
    glVertex2f(x0,y0);
    glVertex2f(x1,y1);
    glEnd();
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

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    mydraw();
    glFlush();
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("[-1,1] X [-1,1]     (500 X 500 pixel)");
    glutDisplayFunc(display);
    init();
    glutMainLoop();
}
