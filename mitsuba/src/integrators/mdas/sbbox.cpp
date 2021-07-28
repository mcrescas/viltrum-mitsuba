#include "bbox.h"




template <int DIMENSION>
bool CBBox<DIMENSION>::isProximal(const CBBox<DIMENSION>& BBox) const
{
	// get the axis that coincides
	int ProximalAxis = -1;
	for (int i = 0; i < DIMENSION; ++i)
	{
		if ((this->Min.Elements[i] == BBox.Max.Elements[i]) || (this->Max.Elements[i] == BBox.Min.Elements[i]))
		{
			ProximalAxis = i;
			break;
		}
	}

	if (ProximalAxis == -1)
	{
		return false;
	}
 
	for (int i = 0; i < DIMENSION; ++i)
	{

		if (i == ProximalAxis)
		{
			continue;
		}

		if ((this->Min.Elements[i] > BBox.Max.Elements[i]) || (BBox.Min.Elements[i] > this->Max.Elements[i]))
		{
			return false;
		}

		//if ((this->Min.Elements[i] - BBox.Max.Elements[i]) > 1E-3)
		//{
		//	return false;
		//}
		//	
		//if ((BBox.Min.Elements[i] - this->Max.Elements[i]) > 1E-3)
		//{
		//	return false;
		//}


	}

	return true;
}


template <int DIMENSION>
float CBBox<DIMENSION>::ComputeVolume() const
{
  float Volume = 1.0f;

  for (int i = 0; i < DIMENSION; ++i)
  {
    Volume *= (this->Max.Elements[i] - this->Min.Elements[i]);
  }

  return Volume;
}
 



template <int DIMENSION>
void CBBox<DIMENSION>::Reset()
{
	for (int i = 0; i < DIMENSION; ++i)
	{
		this->Min.Elements[i] = 1E20f;
		this->Max.Elements[i] = -1E20f;
	}
}




template <int DIMENSION>
void CBBox<DIMENSION>::Fit(const CVector<DIMENSION>& Vector)
{
	for (int i = 0; i < DIMENSION; ++i)
	{
		if (this->Min.Elements[i] > Vector.Elements[i]) 
		{
			this->Min.Elements[i] = Vector.Elements[i];
		}
		if (this->Max.Elements[i] < Vector.Elements[i]) 
		{
			this->Max.Elements[i] = Vector.Elements[i];
		}
	}
}




template <int DIMENSION>
bool CBBox<DIMENSION>::Check(const CVector<DIMENSION>& Vector) const
{
	for (int i = 0; i < DIMENSION; ++i)
	{
		if (this->Min.Elements[i] > Vector.Elements[i]) 
		{
			return false;
		}
		if (this->Max.Elements[i] < Vector.Elements[i]) 
		{
			return false;
		}
	}

	return true;
}




template <int DIMENSION>
int CBBox<DIMENSION>::GetLargestAxis() const
{
	int LargestAxis = 0;
	float MaxLength = 0.0f;
	for (int i = 0; i < DIMENSION; ++i)
	{
		float Length = this->Max.Elements[i] - this->Min.Elements[i];
		if (Length > MaxLength)
		{
			LargestAxis = i;
			MaxLength = Length;
		}
	}

	return LargestAxis;
}



template <int DIMENSION>
float CBBox<DIMENSION>::ComputeVolume2D() const
{
  float Volume = 1.0f;

  for (int i = 0; i < 2; ++i)
  {
    Volume *= (this->Max.Elements[i] - this->Min.Elements[i]);
  }

  return Volume;
}


template <int DIMENSION>
int CBBox<DIMENSION>::GetLargestAxisID() const
{
	int LargestAxis = 0;
	float MaxLength = 0.0f;
	for (int i = 2; i < DIMENSION; ++i)
	{
		float Length = this->Max.Elements[i] - this->Min.Elements[i];
		if (Length > MaxLength)
		{
			LargestAxis = i;
			MaxLength = Length;
		}
	}

	return LargestAxis;
}

template <int DIMENSION>
int CBBox<DIMENSION>::GetLargestAxis2D() const
{
	int LargestAxis = 0;
	float MaxLength = 0.0f;
	for (int i = 0; i < 2; ++i)
	{
		float Length = this->Max.Elements[i] - this->Min.Elements[i];
		if (Length > MaxLength)
		{
			LargestAxis = i;
			MaxLength = Length;
		}
	}

	return LargestAxis;
}
