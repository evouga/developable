#ifndef AUTODIFFTEMPLATES_H
#define AUTODIFFTEMPLATES_H

template<typename T> void diff(const T *vec1, const T *vec2, T *difference)
{
    for(int i=0; i<3; i++)
        difference[i] = vec1[i]-vec2[i];
}

template<typename T> void diff2D(const T *vec1, const T *vec2, T *difference)
{
    for(int i=0; i<2; i++)
        difference[i] = vec1[i]-vec2[i];
}

template<typename T> T norm(const T *vec)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec[i]*vec[i];
    return sqrt(result);
}

template<typename T> T norm2D(const T *vec)
{
    T result = 0;
    for(int i=0; i<2; i++)
        result += vec[i]*vec[i];
    return sqrt(result);
}

template<typename T> void cross(const T *vec1, const T *vec2, T *result)
{
    result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

template<typename T> void normalize(T *vec)
{
    T n = norm(vec);
    for(int i=0; i<3; i++)
        vec[i] /= n;
}

template<typename T> T dot(const T *vec1, const T *vec2)
{
    T result = 0;
    for(int i=0; i<3; i++)
        result += vec1[i]*vec2[i];
    return result;
}

template<typename T> T dot2D(const T *vec1, const T *vec2)
{
    T result = 0;
    for(int i=0; i<2; i++)
        result += vec1[i]*vec2[i];
    return result;
}

#endif // AUTODIFFTEMPLATES_H
