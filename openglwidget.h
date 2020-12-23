#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include <QWidget>
#include <QOpenGLWidget>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/glut.h>
#include <iostream>
#include "src/Mesh.h"
#include "src/GLProjector.h"
#include "src/Deformer.h"
#include "src/Parameterization.h"
#include <commdlg.h>
#include <string>
#include <opencv2/opencv.hpp>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

using namespace cv;
namespace Ui {
class openglwidget;
}

class openglwidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    explicit openglwidget(QWidget *parent = 0);
    ~openglwidget();
    void Changeimportpath();
    void GetIteration(int iter);
    void changeintocolorsmoothshaded();
    void changeintohiddenline();
    void changeintoflatshaded();
    void changeintowireframe();
private:
    Ui::openglwidget *ui;

    QBasicTimer timer;
    QMatrix4x4 projection;

    QVector2D diff;
    QVector2D mousePressPosition;
    QVector3D rotationAxis;
    qreal angularSpeed;
    QQuaternion rotation;
    double sphi = 90.0, stheta = 45.0, sdepth = 10;     // for simple trackball

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void paintGL(String *path2);

    void DisplayFunc();
    void DrawColorSmoothShaded();
    void DrawHiddenLine();
    void DrawFlatShaded();
    void DrawWireframe();

    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void timerEvent(QTimerEvent *e) override;
    Mesh mesh;
    int displaymode=0;
    int lastx=0,lasty=0;
    int viewiter=0, showstart = 0, showend = 0;

};

#endif // OPENGLWIDGET_H
