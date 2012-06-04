#include "zoomer.h"
#include "camera.h"

using namespace Eigen;

Zoomer::Zoomer(Camera &c, double scale) : c_(c), zooming_(false), scale_(scale)
{
    startPos_.setZero();
}

void Zoomer::setScale(double scale)
{
    scale_ = scale;
}

void Zoomer::startZoom(const Vector2d &pos)
{
    zooming_ = true;
    startPos_ = pos;
}

void Zoomer::updateZoom(const Vector2d &pos)
{
    if(!zooming_)
        return;

    Vector3d right, up, center;
    c_.getSpanningSet(right, up, center);
    Vector3d translation = center*scale_*(pos[1]-startPos_[1]);
    c_.translateEye(translation);
    startPos_ = pos;
}

void Zoomer::stopZoom()
{
    zooming_ = false;
}
