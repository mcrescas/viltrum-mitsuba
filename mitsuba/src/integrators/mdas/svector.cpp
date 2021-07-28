#include "vector.h"




template <int DIMENSION>
void CVector<DIMENSION>::UniformRandomHyperSpehere(MTRand& mtrand)
{
  do
  {
    for (int i = 0; i < DIMENSION; i++)
    {
      this->Elements[i] = 2.0f * (float)mtrand.rand() - 1.0f;
    }
  }
  while (this->DotProduct((*this)) > 1.0f);
}




template <int DIMENSION>
void CVector<DIMENSION>::Normalize()
{
	float Length = sqrtf(this->DotProduct(*this));
  float Scale = 1.0f / (Length + 1E-20f);

  for (int i = 0; i < DIMENSION; i++)
  {
    this->Elements[i] = this->Elements[i] * Scale;
  }
}




template <int DIMENSION>
CVector<DIMENSION> CVector<DIMENSION>::operator*(float Scalar) const
{
  CVector Result;
  for (int i = 0; i < DIMENSION; i++)
  {
    Result.Elements[i] = this->Elements[i] * Scalar;
  }

  return Result;
}




template <int DIMENSION>
CVector<DIMENSION> CVector<DIMENSION>::operator+(const CVector<DIMENSION>& Vector) const 
{
  CVector<DIMENSION> Result;
  for (int i = 0; i < DIMENSION; i++)
  {
    Result.Elements[i] = this->Elements[i] + Vector.Elements[i];
  }

  return Result;
}







 template <int DIMENSION>
float CVector<DIMENSION>::SquaredDistance(const CVector<DIMENSION>& Vector) const
{
  CVector<DIMENSION> Difference = (*this) - Vector;
  return Difference.DotProduct(Difference);
}




template <int DIMENSION>
float CVector<DIMENSION>::Distance(const CVector<DIMENSION>& Vector) const
{
  return sqrtf(this->SquaredDistance(Vector));
}
