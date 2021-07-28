#include "bsphere.h"




template <int DIMENSION>
void CBSphere<DIMENSION>::Combine(const CBSphere<DIMENSION>& BSphere)
{
	CVector<DIMENSION> Offset = BSphere.Center - this->Center;
	Offset.Normalize();

	bool isAllZero = true;
	for (int i = 0; i < DIMENSION; ++i)
	{
		if (Offset.Elements[i] != 0.0f)
		{
			isAllZero = false;
			break;
		}
	}
	if (isAllZero)
	{
		Offset.Elements[0] = 1.0f;
	}

	CVector<DIMENSION> Vector0 = BSphere.Center + Offset * sqrtf(BSphere.SquaredRadius);
	float SquaredDistance = Vector0.SquaredDistance(this->Center);

	if (SquaredDistance > this->SquaredRadius)
	{
		this->SquaredRadius = SquaredDistance;
	}
}




template <int DIMENSION>
bool CBSphere<DIMENSION>::Check(const CVector<DIMENSION>& Vector) const
{
	float SquaredDistance = Vector.SquaredDistance(this->Center);

	if (SquaredDistance > this->SquaredRadius)
	{
		return false;
	}
	else
	{
		return true;
	}
}




template <int DIMENSION>
void CBSphere<DIMENSION>::Expand(const CVector<DIMENSION>& Vector)
{
	float SquaredDistance = Vector.SquaredDistance(this->Center);

	if (SquaredDistance > this->SquaredRadius)
	{
		this->SquaredRadius = SquaredDistance;
	}
}
