#ifndef ZOOMER_H
#define ZOOMER_H

#include <Eigen/Core>

class Camera;

class Zoomer
{
public:
    Zoomer(Camera &c, double scale);

    void setScale(double scale);
    void startZoom(const Eigen::Vector2d &pos);
    void updateZoom(const Eigen::Vector2d &pos);
    void stopZoom();

private:
    Camera &c_;
    bool zooming_;
    Eigen::Vector2d startPos_;
    double scale_;
};

#endif // ZOOMER_H
