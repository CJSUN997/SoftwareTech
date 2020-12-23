#include "openglwidget.h"
#include "ui_openglwidget.h"
#include<QDebug>
#include<QFileDialog>
#include <QMouseEvent>

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

void openglwidget::Changeimportpath()
{
    QString path=QFileDialog::getOpenFileName(this,"open","../models/","*.obj");
    String s = path.toStdString();
    paintGL(&s);

}

void openglwidget::initializeGL()
{

    glClearColor(1,1,1,1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    //String objpath = "../models/cube.obj";
    String objpath = "E:/Courses/SoftwareTech/build-Cube-Qt_5_14_1_msvc2017-Debug/models/cube.obj";
    const char *charpath = objpath.data();
    mesh.LoadObjFile(charpath);

    mesh.SetIteration(6);
    GetIteration(0);
//    timer.start(12, this);

}
void openglwidget::DisplayFunc()
{

}
// Wireframe render function
void openglwidget::DrawWireframe() {
    HEdgeList heList = mesh.Edges();
    HEdgeList bheList = mesh.BoundaryEdges();
    FaceList fList = mesh.Faces();
    glColor3f(1.0f, 0.0f, 0.0f);
//    glRotatef(rotation.scalar(),rotation.x(),rotation.y(),rotation.z());

    glBegin(GL_LINES);

    for (size_t i = showstart; i < showend; i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
        glVertex3dv(pos3.data());
        glVertex3dv(pos1.data());
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);
}

// Hidden Line render function
void openglwidget::DrawHiddenLine() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glColor3f(0.5, 0.5, 0);
//    glRotatef(rotation.scalar(),rotation.x(),rotation.y(),rotation.z());

    glBegin(GL_TRIANGLES);
    for (size_t i = showstart; i < showend; i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());

    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    DrawWireframe();
}

void openglwidget::DrawFlatShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
//    glRotatef(rotation.scalar(),rotation.x(),rotation.y(),rotation.z());

    glBegin(GL_TRIANGLES);

    for (size_t i = showstart; i < showend; i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        Eigen::Vector3d normal = (pos2 - pos1).cross(pos3 - pos1);
        normal.normalize();
        glNormal3dv(normal.data());
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}


void openglwidget::DrawColorSmoothShaded()
{
    FaceList fList = mesh.Faces();
       glShadeModel(GL_SMOOTH);
       glEnable(GL_LIGHTING);
       glColor3f(0.5f, 0.5f, 0.5f);
       //qDebug() <<"mouse event:" << rotation.scalar() <<" " << rotation.x() ;

       //glRotatef(rotation.scalar(),rotation.x(),rotation.y(),rotation.z());
       glBegin(GL_TRIANGLES);

    for (size_t i = showstart; i < showend; i++) {
           Face* f = fList[i];
           Vertex* v1 = f->HalfEdge()->Start();
           Vertex* v2 = f->HalfEdge()->End();
           Vertex* v3 = f->HalfEdge()->Next()->End();
           glNormal3dv(v1->Normal().data());

           //glColor3dv(v1->Color().data());
           glVertex3dv(v1->Position().data());
           glNormal3dv(v2->Normal().data());
           //glColor3dv(v2->Color().data());
           glVertex3dv(v2->Position().data());
           glNormal3dv(v3->Normal().data());
           //glColor3dv(v3->Color().data());
           glVertex3dv(v3->Position().data());
       }
       glEnd();
       glDisable(GL_LIGHTING);
}

void openglwidget::paintGL()
{
    // Calculate model view transformation
    QMatrix4x4 matrix;
    matrix.translate(0.0, 0.0, -5.0);
    matrix.rotate(rotation);
    glRotatef(stheta/128.0, 1.0, 0.0, 0.0);
    glRotatef(sphi/128.0, 0.0, 1.0, 0.0);
    switch(displaymode){
        case 0:
            DrawColorSmoothShaded();
            break;
        case 1:
            DrawFlatShaded();
            break;
        case 2:
            DrawHiddenLine();
            break;
        case 3:
            DrawWireframe();
            break;

    }

}

void openglwidget::paintGL(String *path2)
{
    const char *charpath = path2->data();
    mesh.LoadObjFile(charpath);
    mesh.SetIteration(6);
    GetIteration(0);
    switch(displaymode){
        case 0:
            DrawColorSmoothShaded();
            break;
        case 1:
            DrawFlatShaded();
            break;
        case 2:
            DrawHiddenLine();
            break;
        case 3:
            DrawWireframe();
            break;

    }
    update();

}

void openglwidget::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //gluPerspective(45, (float)w/h, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //gluLookAt(0,0,5,0,0,0,0,1,0);
}

void openglwidget::mousePressEvent(QMouseEvent *e)
{
    // Save mouse press position
    mousePressPosition = QVector2D(e->localPos());
}
//
void openglwidget::mouseReleaseEvent(QMouseEvent *e)
{
    // Mouse release position - mouse press position
//    QVector2D diff = QVector2D(e->localPos()) - mousePressPosition ;

//    sphi += (double)diff.x() / 4.0;
//    stheta += (double)diff.y() / 4.0;
//    // Rotation axis is perpendicular to the mouse position difference
//    // vector
//    QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();

//    // Accelerate angular speed relative to the length of the mouse sweep
//    qreal acc = diff.length() / 100.0;

//    // Calculate new rotation axis as weighted sum
//    rotationAxis = (rotationAxis * angularSpeed + n * acc).normalized();

//    // Increase angular speed
//    angularSpeed += acc;
//    update();
}
void openglwidget::mouseMoveEvent(QMouseEvent *event)
{
//    QVector2D diff = QVector2D(event->localPos()) - mousePressPosition ;
//    qDebug()<<diff.x() ;
    sphi += (double)(event->localPos().x() - lastx) ;
    stheta += (double)(event->localPos().y() - lasty);
    lastx = event->localPos().x();
     lasty = event->localPos().y();
    update();
}
//! [0]

//! [1]
void openglwidget::timerEvent(QTimerEvent *)
{
    // Decrease angular speed (friction)
    angularSpeed *= 0.99;

    // Stop rotation when speed goes below threshold
    if (angularSpeed < 0.01) {
        angularSpeed = 0.0;
    } else {
        // Update rotation
        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angularSpeed) * rotation;

        // Request an update
        update();

    }
}
void openglwidget::changeintocolorsmoothshaded()
{
    displaymode=0;
    update();
}
void openglwidget::changeintohiddenline()
{
    displaymode=2;
    update();
}
void openglwidget::changeintoflatshaded()
{
    displaymode=1;
    update();
}
void openglwidget::changeintowireframe()
{
    displaymode=3;
    update();
}
//当前显示的代数
void openglwidget::GetIteration(int iter){
    viewiter=iter;
    showstart = showend =0;
    for(int j = 0; j < viewiter;j++){
        showstart += mesh.origialFaceNum*pow(4,j);

    }
    for(int j = 0; j < viewiter+1;j++){

        showend  += mesh.origialFaceNum*pow(4,j);
    }
    qDebug() <<showstart <<" " << showend;
    update();
}
