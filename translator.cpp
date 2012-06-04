#include "translator.h"
#include "camera.h"

using namespace Eigen;

Translator::Translator(Camera &c, double scale) : c_(c), translating_(false), scale_(scale)
{
    translation_.setZero();
}

void Translator::setScale(double scale)
{
    scale_ = scale;
}


void Translator::startTranslation(const Vector2d &pos)
{
    translating_ = true;
    startpos_ = pos;
}

void Translator::updateTranslation(const Vector2d &pos)
{
    if(!translating_)
        return;
    Vector3d right, up, center;
    c_.getSpanningSet(right, up, center);
    Vector2d v = (pos-startpos_)*scale_;
    translation_ = v[0] * right + v[1] *  up;
    c_.translateEye(-translation_);
    c_.translateCenter(-translation_);
    startpos_ = pos;
}

void Translator::stopTranslation()
{
    translating_ = false;
}

const Vector3d &Translator::getTranslation() const
{
    return translation_;
}
