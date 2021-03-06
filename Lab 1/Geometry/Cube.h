#ifndef __CUBE_H__
#define __CUBE_H__

#include "Geometry/Quadric.h"
#include <vector>

class Cube : public Implicit {
public:
  Cube();
  ~Cube();

  virtual float GetValue(float x, float y, float z) const;

private:
  std::vector<Quadric *> mPlanes;
  bool mEuclideanDistance;
};

#endif
