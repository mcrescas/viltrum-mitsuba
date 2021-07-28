#ifndef __VECTOR_H
#define __VECTOR_H

#include "global.h"
#include <math.h>
#include "mt.h"




template <int DIMENSION>
class CVector
{
  public:
    float Elements[DIMENSION];
 
    CVector operator+(const CVector& Vector) const;
    inline CVector operator-(const CVector& Vector) const;
    CVector operator*(float Scalar) const;

		CVector()
		{
			for (int i = 0; i < DIMENSION; i++)
			{
				this->Elements[i] = 0.0f;
			}
		}

    inline float DotProduct(const CVector& Vector) const;
    float Distance(const CVector& Vector) const;
    float SquaredDistance(const CVector& Vector) const;
		void Normalize();
    void UniformRandomHyperSpehere(MTRand& mtrand);

  private:

};


 


template <int DIMENSION>
inline CVector<DIMENSION> operator*(float Scalar, const CVector<DIMENSION>& Vector)
{
  CVector<DIMENSION> Result;
  for (int i = 0; i < DIMENSION; i++)
  {
    Result.Elements[i] = Vector.Elements[i] * Scalar;
  }

  return Result;
}


template <int DIMENSION>
inline CVector<DIMENSION> CVector<DIMENSION>::operator-(const CVector<DIMENSION>& Vector) const 
{
  CVector Result;
  for (int i = 0; i < DIMENSION; i++)
  {
    Result.Elements[i] = this->Elements[i] - Vector.Elements[i];
  }

  return Result;
}




template <int DIMENSION>
inline float CVector<DIMENSION>::DotProduct(const CVector<DIMENSION>& Vector) const
{
  float Result = 0.0f;

  for (int i = 0; i < DIMENSION; i++)
  {
    Result = Result + this->Elements[i] * Vector.Elements[i];
  }

  return Result;
}



template <int DIMENSION>
class CMatrix
{
  public:
    CVector<DIMENSION> Rows[DIMENSION];

    CMatrix()
    {
      // empty
    }
};

#include "svector.cpp"

#endif
