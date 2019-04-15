/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Söderström (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/

#include "LoopSubdivisionMesh.h"
#include <cassert>

/*! Subdivides the mesh uniformly one step
 */
void LoopSubdivisionMesh::Subdivide() {
  // Create new mesh and copy all the attributes
  HalfEdgeMesh subDivMesh;
  subDivMesh.SetTransform(GetTransform());
  subDivMesh.SetName(GetName());
  subDivMesh.SetColorMap(GetColorMap());
  subDivMesh.SetWireframe(GetWireframe());
  subDivMesh.SetShowNormals(GetShowNormals());
  subDivMesh.SetOpacity(GetOpacity());
  if (IsHovering())
    subDivMesh.Hover();
  if (IsSelected())
    subDivMesh.Select();
  subDivMesh.mMinCMap = mMinCMap;
  subDivMesh.mMaxCMap = mMaxCMap;
  subDivMesh.mAutoMinMax = mAutoMinMax;

  // loop over each face and create 4 new ones
  for (size_t i = 0; i < mFaces.size(); i++) {
    // subdivide face
    std::vector<std::vector<Vector3<float>>> faces = Subdivide(i);

    // add new faces to subDivMesh
    for (size_t j = 0; j < faces.size(); j++) {
      subDivMesh.AddFace(faces.at(j));
    }
  }

  // Assigns the new mesh
  *this = LoopSubdivisionMesh(subDivMesh, ++mNumSubDivs);
  Update();
}

/*! Subdivides the face at faceindex into a vector of faces
 */
std::vector<std::vector<Vector3<float>>>
LoopSubdivisionMesh::Subdivide(size_t faceIndex) {
  std::vector<std::vector<Vector3<float>>> faces;
  EdgeIterator eit = GetEdgeIterator(f(faceIndex).edge);

  // get the inner halfedges
  size_t e0, e1, e2;
  // and their vertex indices
  size_t v0, v1, v2;

  e0 = eit.GetEdgeIndex();
  v0 = eit.GetEdgeVertexIndex();
  eit.Next();
  e1 = eit.GetEdgeIndex();
  v1 = eit.GetEdgeVertexIndex();
  eit.Next();
  e2 = eit.GetEdgeIndex();
  v2 = eit.GetEdgeVertexIndex();

  // Compute positions of the vertices
  Vector3<float> pn0 = VertexRule(v0);
  Vector3<float> pn1 = VertexRule(v1);
  Vector3<float> pn2 = VertexRule(v2);

  // Compute positions of the edge vertices
  Vector3<float> pn3 = EdgeRule(e0);
  Vector3<float> pn4 = EdgeRule(e1);
  Vector3<float> pn5 = EdgeRule(e2);

  // add the four new triangles to new mesh
  std::vector<Vector3<float>> verts;
  verts.push_back(pn0);
  verts.push_back(pn3);
  verts.push_back(pn5);
  faces.push_back(verts);
  verts.clear();
  verts.push_back(pn3);
  verts.push_back(pn4);
  verts.push_back(pn5);
  faces.push_back(verts);
  verts.clear();
  verts.push_back(pn3);
  verts.push_back(pn1);
  verts.push_back(pn4);
  faces.push_back(verts);
  verts.clear();
  verts.push_back(pn5);
  verts.push_back(pn4);
  verts.push_back(pn2);
  faces.push_back(verts);
  return faces;
}

/*! Computes a new vertex, replacing a vertex in the old mesh
 */
Vector3<float> LoopSubdivisionMesh::VertexRule(size_t vertexIndex) {
  // Get the current vertex
  Vector3<float> vtx = v(vertexIndex).pos;
  // Get the valence from size of 1-ring
  std::vector<size_t> oneRing = HalfEdgeMesh::FindNeighborVertices(vertexIndex);
  float beta = Beta(oneRing.size());	// Calculate beta

  // Update vtx with the weight (fig 17. 1-k*beta)
  vtx *= (1 - oneRing.size() * beta);

  // Update weight for the neighbouring vertices
  for (int i = 0; i < oneRing.size(); ++i)
  {
	  vtx += v(oneRing.at(i)).pos * beta;
  }
  return vtx;
}

/*! Computes a new vertex, placed along an edge in the old mesh
 */
Vector3<float> LoopSubdivisionMesh::EdgeRule(size_t edgeIndex) {

  // Place the edge vertex halfway along the edge
  HalfEdge &e0 = e(edgeIndex);
  HalfEdge &e1 = e(e0.pair);
  Vector3<float> &v0 = v(e0.vert).pos;			// Left
  Vector3<float> &v1 = v(e1.vert).pos;			// Right
  Vector3<float> &v2 = v(e(e0.prev).vert).pos;	// Up
  Vector3<float> &v3 = v(e(e1.prev).vert).pos;	// Down 
  return (3.0f / 8.0f) * (v0 + v1) + (1.0f / 8.0f) * (v2 + v3);	// Split up fractions for proper precision
  
}

//! Return weights for interior verts
float LoopSubdivisionMesh::Beta(size_t valence) {
  if (valence == 6) {
    return 1.0f / 16.0f;
  } else if (valence == 3) {
    return 3.0f / 16.0f;
  } else {
    return 3.0f / (8.0f * valence);
  }
}
