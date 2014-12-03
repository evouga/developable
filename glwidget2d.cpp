#include "glwidget2d.h"
#include "controller.h"
#include <iostream>

using namespace std;

GLWidget2D::GLWidget2D(QWidget *parent) :
    QGLWidget(QGLFormat(QGL::SampleBuffers),parent), cont_(NULL), center_(0,0), radius_(1.0)
{
}

void GLWidget2D::setController(Controller &cont)
{
    cont_ = &cont;
    centerCamera();
}

void GLWidget2D::centerCamera()
{
    assert(cont_);
    cont_->getMaterialBounds(center_, radius_);
}


void GLWidget2D::initializeGL()
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glDisable(GL_DEPTH_TEST);
}

void GLWidget2D::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
}

void GLWidget2D::paintGL()
{
    assert(cont_);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(center_[0]-radius_, center_[0]+radius_, center_[1]-radius_, center_[1]+radius_, 0, 1);
    glMatrixMode(GL_MODELVIEW);

    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear (GL_COLOR_BUFFER_BIT);
    glColor3f (0.0, 0.0, 0.0);

    cont_->renderMaterial();
}
