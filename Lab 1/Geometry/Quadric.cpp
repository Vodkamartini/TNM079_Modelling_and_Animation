/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Quadric.h"

Quadric::Quadric(const Matrix4x4<float> &q) { this->mQuadric = q; }

Quadric::~Quadric() {}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const {

    // Calculate p prime
    Vector4<float> p(x, y, z, 1);
    Vector4<float> p_prime = mWorld2Obj * p;

    // Calculate Q prime
    // Matrix4x4<float> Q_prime = mWorld2Obj.Inverse() * (mQuadric * mWorld2Obj.Inverse());

    Matrix4x4<float> Q_prime = GetTransform() * mQuadric;

    // Calculate value
    return p * (Q_prime * p);
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */
Vector3<float> Quadric::GetGradient(float x, float y, float z) const {

    Vector4<float> p(x, y, z, 1);

    float gradientMatrix[4][4] = {
        {2 * mQuadric(0, 0), 2 * mQuadric(0, 1), 2 * mQuadric(0, 2), 2 * mQuadric(0, 3)},
        {2 * mQuadric(1, 0), 2 * mQuadric(1, 1), 2 * mQuadric(1, 2), 2 * mQuadric(1, 3)},
        {2 * mQuadric(2, 0), 2 * mQuadric(2, 1), 2 * mQuadric(2, 2), 2 * mQuadric(2, 3)},
    };

    Matrix4x4<float> Qsub(gradientMatrix);
    Vector4<float> product = Qsub * p;

    return Vector3<float>(product[1], product[2], product[3]);
}
