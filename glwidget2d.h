#ifndef GLWIDGET2D_H
#define GLWIDGET2D_H

#include <QGLWidget>
#include <Eigen/Core>

class Controller;

class GLWidget2D : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget2D(QWidget *parent = 0);

    void setController(Controller &cont);
    void centerCamera();

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

private:
    Controller *cont_;
    Eigen::Vector2d center_;
    double radius_;

signals:

public slots:



};

#endif // GLWIDGET2D_H
