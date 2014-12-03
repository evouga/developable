#include "glwidget.h"
#include <iostream>
#include <QtGui>
#include <Eigen/Core>
#include "yimage.h"
#include "controller.h"
#include <iostream>
#include <GL/glu.h>

using namespace std;
using namespace Eigen;

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(QGLFormat(QGL::SampleBuffers),parent), cont_(NULL), c_(), translator_(c_, 1.0), rotator_(c_), zoomer_(c_, 10.0), takeScreenshot_(false)
{
}

void GLWidget::setController(Controller &cont)
{
    cont_ = &cont;
    centerCamera();
}

void GLWidget::initializeGL()
{
    glShadeModel(GL_SMOOTH);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClearDepth(1.0);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
}

void GLWidget::resizeGL(int w, int h)
{
    c_.setPerpective(60.0, 1.0);
    c_.setViewport(w, h);
}

void GLWidget::paintGL()
{
    assert(cont_);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3f (0.0, 0.0, 0.0);

    c_.applyViewport();
    c_.applyProjection();
    c_.applyLookAt();

    static GLfloat lightPosition[4] = { 0.0, -1.0, 0.0, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

    cont_->renderMesh();

    if(takeScreenshot_)
    {
        takeScreenshot_ = false;
        GLint view[4] ;
        glGetIntegerv( GL_VIEWPORT, view ) ;

        YImage img ;
        img.resize( view[2], view[3] ) ;
        //glReadBuffer( GL_FRONT );
        glReadPixels( view[0], view[1], view[2], view[3], GL_RGBA, GL_UNSIGNED_BYTE, img.data() ) ;
        //glReadBuffer( GL_BACK );

        img.flip() ;
        img.save( ssFilename_.c_str() ) ;
    }
}

void GLWidget::scaleMousePos(int x, int y, double &scaledx, double &scaledy) const
{
    int w, h;
    c_.getViewport(w,h);
    scaledx = 2.0 * x / double(w-1) - 1.0;
    scaledy = 2.0 * (h - y -1) / double(h-1) - 1.0;
}

GLWidget::MouseAction GLWidget::deduceAction(QMouseEvent *event)
{
    if(event->buttons() & Qt::LeftButton)
    {
        if(event->modifiers() & Qt::ShiftModifier)
            return MA_TRANSLATE;
        return MA_ROTATE;
    }
    else if(event->buttons() & Qt::RightButton)
        return MA_ZOOM;
    return MA_NONE;
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{

    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);

    MouseAction ma = deduceAction(event);

    switch(ma)
    {
        case MA_TRANSLATE:
        {
            translator_.startTranslation(pos);
            break;
        }
        case MA_ROTATE:
        {
            rotator_.startRotation(pos);
            break;
        }
        case MA_ZOOM:
        {
            zoomer_.startZoom(pos);
            break;
        }
        default:
            break;
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);
    translator_.updateTranslation(pos);
    rotator_.updateRotation(pos);
    zoomer_.updateZoom(pos);
    updateGL();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    int x = event->pos().x();
    int y = event->pos().y();
    Vector2d pos;
    scaleMousePos(x,y,pos[0],pos[1]);
    translator_.stopTranslation();
    rotator_.stopRotation();
    zoomer_.stopZoom();
}

void GLWidget::centerCamera()
{
    assert(cont_);
    Vector3d centroid;
    double radius;
    cont_->getSceneBounds(centroid, radius);

    c_.setCenter(centroid);
    translator_.setScale(2*radius);
    zoomer_.setScale(2*radius);

    c_.setDefault3D(2*radius);
    lightPos_ = c_.getEye();
}

void GLWidget::saveScreenshot(const std::string &filename)
{
    takeScreenshot_ = true;
    ssFilename_ = filename;
}
