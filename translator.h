#ifndef TRANSLATOR_H
#define TRANSLATOR_H

#include <Eigen/Core>

class Camera;

class Translator
{
public:
    Translator(Camera &c, double scale);

    void setScale(double scale);

    void startTranslation(const Eigen::Vector2d &pos);
    void updateTranslation(const Eigen::Vector2d &pos);
    void stopTranslation();

    const Eigen::Vector3d &getTranslation() const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Camera &c_;
    bool translating_;
    double scale_;

    Eigen::Vector2d startpos_;
    Eigen::Vector3d translation_;
};

#endif //TRANSLATOR_H
