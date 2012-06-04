#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "rotator.h"
#include "zoomer.h"
#include "translator.h"
#include "camera.h"
#include <QGLWidget>
#include <Eigen/Core>

class Controller;

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(QWidget *parent = 0);

    void setController(Controller &cont);
    void centerCamera();
    void saveScreenshot(const std::string &filename);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    enum MouseAction { MA_NONE, MA_TRANSLATE, MA_ROTATE, MA_ZOOM};

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void scaleMousePos(int x, int y, double &scaledx, double &scaledy) const;
    MouseAction deduceAction(QMouseEvent *event);

private:
    Controller *cont_;

    Camera c_;
    Translator translator_;
    Rotator rotator_;
    Zoomer zoomer_;
    Eigen::Vector3d lightPos_;

    bool takeScreenshot_;
    std::string ssFilename_;

signals:

public slots:
};

#endif
