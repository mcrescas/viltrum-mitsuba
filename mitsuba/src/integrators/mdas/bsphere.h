#ifndef __BSPHERE_H
#define __BSPHERE_H

#include "global.h"
#include "vector.h"

#include <vector>

template <int DIMENSION>
class CBSphere
{
  public:
    CVector<DIMENSION> Center;
		float SquaredRadius;

	  void Combine(const CBSphere<DIMENSION>& BSphere);
		void Expand(const CVector<DIMENSION>& Vector);
		bool Check(const CVector<DIMENSION>& Vector) const;

  private:

};

#include "sbsphere.cpp"


#endif
