
#include "UniformCubicSplineSubdivisionCurve.h"

UniformCubicSplineSubdivisionCurve::UniformCubicSplineSubdivisionCurve(
    const std::vector<Vector3<float>> &joints, Vector3<float> lineColor, float lineWidth)
    : mCoefficients(joints), mControlPolygon(joints) {
    this->mLineColor = lineColor;
    this->mLineWidth = lineWidth;
}

void UniformCubicSplineSubdivisionCurve::Subdivide() {
    // Allocate space for new coefficients
    std::vector<Vector3<float>> newc;

    assert(mCoefficients.size() > 4 && "Need at least 5 points to subdivide");

    /* Set the first and last control points to be the same, don't forget the padding. 
	 * (For a more detailed explaination, see Lab 3 script p.10 equations (30) and (31),
	 * and figure 14 and 15).*/

	/* Define two temporary control points according to equation (30)
	 * cold is the old control points
	 * cnew is the new control points */
    Vector3<float> cold, cnew;
	
	// Insert the padding points
	newc.push_back(mCoefficients.at(0));												// -2 (see figure 15)
	cnew = (1.0f / 8.0f) * (4.0f * mCoefficients.at(0) + 4.0f * mCoefficients.at(1));
    newc.push_back(cnew);																// -1 (see figure 15)

    // Insert rest of points
    for (int i = 1; i < mCoefficients.size() - 1; i++) {

        // Push back old control point ci
        cold =(1.0f / 8.0f) * (mCoefficients.at(i - 1) + 6.0f * mCoefficients.at(i) + mCoefficients.at(i + 1));
        newc.push_back(cold);

        // Push back new control point c(i+1/2)
        cnew = (1.0f / 8.0f) * (4.0f * mCoefficients.at(i) + 4.0f * mCoefficients.at(i + 1));
        newc.push_back(cnew);
    }
	
	// Insert last control point, the other padding point, 1, is added in the loop
    newc.push_back(mCoefficients.at(mCoefficients.size() - 1));						   // 2 (see figure 15)
	
	// Control check if the number of new control points equals 2 * (number of old points) -1
    assert(newc.size() == ((2 * mCoefficients.size()) -1) && "Incorrect number of new coefficients!");

	// Overwrite the old control points with the new ones
    mCoefficients = newc;
}

void UniformCubicSplineSubdivisionCurve::Render() {
    // Apply transform
    glPushMatrix();  // Push modelview matrix onto stack

    // Convert transform-matrix to format matching GL matrix format
    // Load transform into modelview matrix
    glMultMatrixf(mTransform.ToGLMatrix().GetArrayPtr());

    mControlPolygon.Render();

    // save line point and color states
    glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

    // draw segments
    glLineWidth(mLineWidth);
    glColor3fv(mLineColor.GetArrayPtr());
    glBegin(GL_LINE_STRIP);
    // just draw the spline as a series of connected linear segments
    for (size_t i = 0; i < mCoefficients.size(); i++) {
        glVertex3fv(mCoefficients.at(i).GetArrayPtr());
    }
    glEnd();

    // restore attribs
    glPopAttrib();

    glPopMatrix();

    GLObject::Render();
}
