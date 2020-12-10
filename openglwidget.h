#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include <QWidget>
#include <QOpenGLWidget>
#include <GL/gl.h>
#include <GL/glu.h>
#include <opencv2/opencv.hpp>
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

private:
    Ui::openglwidget *ui;

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
};

#endif // OPENGLWIDGET_H
