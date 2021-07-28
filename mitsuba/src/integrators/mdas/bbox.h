#ifndef __BBOX_H
#define __BBOX_H

#include "global.h"
#include "vector.h"

#include <vector>

template <int DIMENSION>
class CBBox
{
  public:
    CVector<DIMENSION> Min;
    CVector<DIMENSION> Max;

	  void Reset();
		void Fit(const CVector<DIMENSION>& Vector);
		int GetLargestAxis() const;
		bool Check(const CVector<DIMENSION>& Vector) const;
    float ComputeVolume() const;
		bool isProximal(const CBBox& BBox) const;

		int GetLargestAxis2D() const;
    int GetLargestAxisID() const;
    float ComputeVolume2D() const;

  private:

};

#include "sbbox.cpp"

#endif
