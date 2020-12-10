#include "openglwidget.h"
#include "ui_openglwidget.h"
#include <iostream>
#include "src/Mesh.h"
#include "src/GLProjector.h"
#include "src/Deformer.h"
#include "src/Parameterization.h"
#include <commdlg.h>
#include <string>
#include <opencv2/opencv.hpp>


openglwidget::openglwidget(QWidget *parent) :
    QOpenGLWidget(parent),
    ui(new Ui::openglwidget)
{
    ui->setupUi(this);
}

openglwidget::~openglwidget()
{
    delete ui;
}

void openglwidget::initializeGL()
{
    glClearColor(0,0,0,1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
}

void openglwidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_TRIANGLES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(-0.5, -0.5, 0);
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f( 0.5, -0.5, 0);
        glColor3f(0.0, 0.0, 1.0);
        glVertex3f( 0.0,  0.5, 0);
    glEnd();
}

void openglwidget::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (float)w/h, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,0,0,0,0,1,0);
}
