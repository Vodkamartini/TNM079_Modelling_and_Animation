/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces =
    NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize() {
    // Allocate memory for the quadric array
    size_t numVerts = mVerts.size();
    mQuadrics.reserve(numVerts);
    std::streamsize width = std::cerr.precision();  // store stream precision
    for (size_t i = 0; i < numVerts; i++) {

        // Compute quadric for vertex i here
        mQuadrics.push_back(createQuadricForVert(i));

        // Calculate initial error, should be numerically close to 0

        Vector3<float> v0 = mVerts[i].pos;
        Vector4<float> v(v0[0], v0[1], v0[2], 1);
        Matrix4x4<float> m = mQuadrics.back();

        float error = v * (m * v);
        // std::cerr << std::scientific << std::setprecision(2) << error << " ";
    }
    std::cerr << std::setprecision(width) << std::fixed;  // reset stream precision

    // Run the initialize for the parent class to initialize the edge collapses
    DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute,
 * DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
    // Compute collapse->position and collapse->cost here
    // based on the quadrics at the edge endpoints

    // Get the two vertices
    size_t indx1 = e(collapse->halfEdge).vert;
    size_t indx2 = e(e(collapse->halfEdge).pair).vert;

    // Get the quadric for each vertex
    Matrix4x4<float> Q1 = mQuadrics.at(indx1);
    Matrix4x4<float> Q2 = mQuadrics.at(indx2);

    Matrix4x4<float> Q = Q1 + Q2;

    /* We want to find the vertex newPos that gives
     * us the minimal possible cost. This is solved
     * with the equation Q*newPos = 0. To solve the
     * equation we use the inverse of Q and multiply
     * it with the zero vector. */

	/* Create a matrix Qeq from Q (to be used in the equation)
     * from deriving v over x, y and z
     * and set the fourth row to [0 0 0 1] (w = 1) */
    float floatQ[4][4] = {
		{Q(0, 0), Q(0, 1), Q(0, 2), Q(0, 3)}, 
		{Q(1, 0), Q(1, 1), Q(1, 2), Q(1, 3)},
		{Q(2, 0), Q(2, 1), Q(2, 2), Q(2, 3)}, 
		{0, 0, 0, 1}
	};

    Matrix4x4<float> Qeq(floatQ);

    /* Finding the best vertex newPos is only possible
     * if Qeq is not singular */
    if (!Qeq.IsSingular()) {
        // Create zero-vector used for equation and the vector newPos
        Vector4<float> zeroVec(0, 0, 0, 1);
        Vector4<float> newPos = Qeq.Inverse() * zeroVec;

        // Compute position for collapse
        collapse->position = Vector3<float>(newPos[0], newPos[1], newPos[2]);

        // Compute cost for collapse at collapse-position
        //Vector4<float> v = {newPos[0], newPos[1], newPos[2], 1};
        collapse->cost = newPos * (Q * newPos);
    }
    /* In the unfortunate case where Qeq is singular (not inversible)
     * we decide the position of the collapse by comparing the cost of
     * moving v1 to v2, v2 to v1 or using an intermediate vertex
     * v = (v1 + v2 )/ 2 */
    else {
        // Case 1, we use v1->v2
        Vector3<float> v1 = v(indx1).pos;
        Vector4<float> vec1 = {v1[0], v1[1], v1[2], 1};
        float cost1 = vec1 * (Q * vec1);

        // Case 2, we use v2->v1
        Vector3<float> v2 = v(indx2).pos;
        Vector4<float> vec2 = {v2[0], v2[1], v2[2], 1};
        float cost2 = vec2 * (Q * vec2);

        // Case 3, we use v3 = (v1 + v2) / 2
        Vector3<float> v3 = ((v(indx1).pos + v(indx2).pos) / 2);
        Vector4<float> vec3 = {v3[0], v3[1], v3[2], 1};
        float cost3 = vec3 * (Q * vec3);

        // Assign the lowest cost to collapse->cost
        collapse->cost = std::min(std::min(cost1, cost2), cost3);  // std::min only takes two arguments

        // Update collapse->position according to which cost was the smallest
        if (collapse->cost == cost1)
            collapse->position = v1;
        else if (collapse->cost == cost2)
            collapse->position = v2;
        else
            collapse->position = v3;
    }

    // std::cerr << "computeCollapse in QuadricDecimationMesh not implemented.\n";
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(size_t ind) {
    DecimationMesh::updateVertexProperties(ind);
    mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForVert(size_t indx) const {

    // Initialize empty quadric
    float q[4][4] = {
		{0, 0, 0, 0}, 
		{0, 0, 0, 0}, 
		{0, 0, 0, 0}, 
		{0, 0, 0, 0}};
    Matrix4x4<float> Q(q);

    // Get the one-ring of the vertex in order to be able to calculate the quadric for the adjacent
    // faces
    std::vector<size_t> oneRing = HalfEdgeMesh::FindNeighborFaces(indx);

    // The quadric for a vertex is the sum of all the quadrics for the adjacent
    // faces Tip: Matrix4x4 has an operator +=
    for (int i = 0; i < oneRing.size(); i++) {

        Q += createQuadricForFace(oneRing.at(i));
    }

    return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForFace(size_t indx) const {

    // Calculate the quadric (outer product of plane parameters) for a face
    // here using the formula from Garland and Heckbert

    // v0 = (x,y,z) is a point in the plane
    Vector3<float> position = v(e(f(indx).edge).vert).pos;

    // n = (a,b,c) is the normal of the plane
    Vector3<float> normal = f(indx).normal;

    /* d is calculated with dot product of position and normal,
     * v0 * n = -d
     * and multiplied with -1.0 according to ax + by + cz = -d	*/
    float d = (position * normal) * -1.0f;

    /*
     * Kp = | (a*a) (a*b) (a*c) (a*d) |
     *		| (b*a) (b*b) (b*c) (b*d) |
     *		| (c*a) (c*b) (c*c) (c*d) |
     *		| (d*d) (a*d) (a*d) (d*d) |
     */
    float Kp[4][4] = {
        {normal[0] * normal[0], normal[0] * normal[1], normal[0] * normal[2], normal[0] * d},
        {normal[1] * normal[0], normal[1] * normal[1], normal[1] * normal[2], normal[1] * d},
        {normal[2] * normal[0], normal[2] * normal[1], normal[2] * normal[2], normal[2] * d},
        {normal[0] * d, normal[1] * d, normal[2] * d, d * d}};

    // Return Kp to be summed with the other Kp-matrices to create Q
    return Matrix4x4<float>(Kp);
}

void QuadricDecimationMesh::Render() {
    DecimationMesh::Render();

    glEnable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);

    if (mVisualizationMode == QuadricIsoSurfaces) {
        // Apply transform
        glPushMatrix();  // Push modelview matrix onto stack

        // Implement the quadric visualization here
        std::cout << "Quadric visualization not implemented" << std::endl;

        // Restore modelview matrix
        glPopMatrix();
    }
}
