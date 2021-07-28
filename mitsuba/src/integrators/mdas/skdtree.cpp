extern "C"
{
	#include "f2c.h"
	#include "clapack.h"
  #undef abs
}

#include "kdtree.h"


#ifndef CELLSPLIT
inline double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (double**)malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = (double*)malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}
#endif



template <int DIMENSION>
CVector<DIMENSION> CKdTree<DIMENSION>::GenerateNewSamplePosition(MTRand& mtrand, bool Is2D)
{
  CVector<DIMENSION> Result;
  float MaxDistance;

#ifdef CELLSPLIT
	// get a parent node
  int ParentNodeIndex = this->Heap.Items[0].Index;

  MaxDistance = 0.0f;
  for (int m = 0; m < this->CandidatesNum; ++m)
  //for (int m = 0; m < 3; ++m)
  {

    CVector<DIMENSION> Candidate;
 
    float Radius = 0.55f * sqrtf(this->Nodes[ParentNodeIndex].BBox.Max.SquaredDistance(this->Nodes[ParentNodeIndex].BBox.Min));
    CVector<DIMENSION> ParentPosition = 0.5f * (this->Nodes[ParentNodeIndex].BBox.Max + this->Nodes[ParentNodeIndex].BBox.Min);

		//CKNNs KNNs;
    //this->KNNSearch(24, ParentPosition, &KNNs);
    //Radius = sqrtf(KNNs.MaxSquaredDistance);
 
    // generate a candidate uniformly within a hypersphere
		bool IsValidCandidate;
    do
		{
			// generate a candidate position
			IsValidCandidate = true;

			CVector<DIMENSION> UniformRandomVector; 
			UniformRandomVector.UniformRandomHyperSpehere(mtrand);
			Candidate = UniformRandomVector * Radius + ParentPosition;

			//// generate a new position within the bounding box
			//for (int i = 0; i < DIMENSION; ++i)
			//{
			//	Candidate.Elements[i] = (this->Nodes[ParentNodeIndex].BBox.Max.Elements[i] - this->Nodes[ParentNodeIndex].BBox.Min.Elements[i]) * mtrand.rand() * (1.0f + 0.2f);
			//}
			//Candidate = Candidate + this->Nodes[ParentNodeIndex].BBox.Min - (this->Nodes[ParentNodeIndex].BBox.Max - this->Nodes[ParentNodeIndex].BBox.Min) * 0.1f;

			// make sure that a candidate position is within the extents of the function
			for (int i = 0; i < DIMENSION; ++i)
			{
				// rejection sampling approach
				if ((Candidate.Elements[i] < 0.0f) || (Candidate.Elements[i] > this->SampleExtent.Elements[i]))
				{
					IsValidCandidate = false;
					break;
				}
			}
    }
		while (!IsValidCandidate);
    //Result = Candidate;


    // compute the minimum distance to the neighbors of the parent sample
    float MinDistance = 1E20f;
    //CKNNs KNNs;

    //this->KNNSearch(1, Candidate, &KNNs);
    int NearestIndex = this->NearestSearch(Candidate);
    MinDistance = this->Samples[NearestIndex].Position.SquaredDistance(Candidate);//KNNs.MaxSquaredDistance;
 
     //take this candidate if it has a larger minimum distance to the existing samples
    if (MinDistance > MaxDistance)
    {
      MaxDistance = MinDistance;
      Result = Candidate;
    }
  }

#else

	// get a parent sample
  int ParentIndex = this->Heap.Items[0].Index;
  const CVector& ParentPosition = this->Samples[ParentIndex].Position;
  float Radius = sqrtf(this->Samples[ParentIndex].KNNs.MaxSquaredDistance);


	// use the best candidate approach to generate a new sample position
  MaxDistance = 0.0f;
  for (int m = 0; m < this->CandidatesNum; ++m)
  {
    // generate a candidate uniformly within a hypersphere
    CVector Candidate;
		bool IsValidCandidate;
    do
		{
			// generate a candidate position
			IsValidCandidate = true;

			CVector UniformRandomVector;
			UniformRandomVector.UniformRandomHyperSpehere(mtrand);
			Candidate = UniformRandomVector * Radius + ParentPosition;


			// make sure that a candidate position is within the extents of the function
			for (int i = 0; i < DIMENSION; ++i)
			{
				// rejection sampling approach
				if ((Candidate.Elements[i] < 0.0f) || (Candidate.Elements[i] > this->SampleExtent.Elements[i]))
				{
					IsValidCandidate = false;
					break;
				}

				// reflection approach (
				//if (Candidate.Elements[i] < 0.0f) Candidate.Elements[i] = -Candidate.Elements[i];
				//if (Candidate.Elements[i] > this->SampleExtent.Elements[i]) Candidate.Elements[i] = this->SampleExtent.Elements[i] + (this->SampleExtent.Elements[i] - Candidate.Elements[i]);
			}
    }
		while (!IsValidCandidate);


    // compute the minimum distance to the neighbors of the parent sample
    float MinDistance = 1E20f;
    CKNNs KNNs;

#ifndef ISOTROPIC
    this->KNNSearch(this->KNNsNum, Candidate, &KNNs);
		for (int i = 0; i < this->KNNsNum; ++i)
		{ 
			//float Distance = this->ComputeDistance(Candidate, KNNs.SamplesIndices[i], KNNs, ParentIndex);
			float Distance = this->ComputeDistance(Candidate, KNNs.SamplesIndices[i], this->Samples[ParentIndex].KNNs, ParentIndex);
      if (Distance < MinDistance)
      {
        MinDistance = Distance;
      }
			//MinDistance = MinDistance + Distance;
    }
#else
		// Euclidean version 
    this->KNNSearch(1, Candidate, &KNNs);
    MinDistance = KNNs.MaxSquaredDistance;
#endif
 
    // take this candidate if it has a larger minimum distance to the existing samples
    if (MinDistance > MaxDistance)
    {
      MaxDistance = MinDistance;
      Result = Candidate;
    }
  }
#endif

  return Result;
}


#ifndef CELLSPLIT
void CKdTree::ComputeSampleKNNs(const int SampleIndex)
{
	// perform KNNs query
	int n = SampleIndex;

	// do KNN query
	const CVector& QueryPosition = this->Samples[n].Position;
	this->KNNSearch(this->KNNsNum + 1, QueryPosition, &(this->Samples[n].KNNs));


	// remove the query sample itself
	int k = 0;
	for (int i = 0; i < (this->KNNsNum + 1); i++)
	{
		if (this->Samples[n].KNNs.SamplesIndices[i] == n)
		{
			k = i;
			break;
		}
	}
	std::vector<int>::iterator it0 = this->Samples[n].KNNs.SamplesIndices.begin();
	for (int j = 0; j < k; j++)
	{
		it0++;
	}
	this->Samples[n].KNNs.SamplesIndices.erase(it0);
	std::vector<float>::iterator it1 = this->Samples[n].KNNs.SamplesDistances.begin();
	for (int j = 0; j < k; j++)
	{
		it1++;
	}
	this->Samples[n].KNNs.SamplesDistances.erase(it1);


	// reconstruct the heap of kNNs
	//this->Samples[n].KNNs.Heap.Clear();
  this->Samples[n].KNNs.Heap.CurrentNumItems = 0;
	for (int i = 0; i < this->KNNsNum; i++)
	{
		this->Samples[n].KNNs.Heap.Put(i, this->Samples[n].KNNs.SamplesDistances[i]);	
	}
}
#endif



#ifndef CELLSPLIT
void CKdTree::ComputeKNNs()
{
	// perform KNNs query for all samples
	int SamplesNum = this->CurrentNumSamples;
	for (int n = 0; n < SamplesNum; n++)
	{
		this->ComputeSampleKNNs(n);
	}
}
#endif




template <int DIMENSION>
void CKdTree<DIMENSION>::SortSamples(std::vector<int>& SamplesIndices, int Axis, int li, int ri)
{
	if ((li - ri) < 50)
	{
		// insertion sort
		for (int i = li + 1; i <= ri; ++i) 
		{
			for (int j = i; j >= (li + 1); --j)
			{
				if (this->Samples[SamplesIndices[j]].Position.Elements[Axis] < this->Samples[SamplesIndices[j - 1]].Position.Elements[Axis])
				{
					int temp = SamplesIndices[j - 1];
					SamplesIndices[j - 1] = SamplesIndices[j];
					SamplesIndices[j] = temp;
				}
			}
		}
	}
	else
	{
		// quick sort
		int i = li;
		int j = ri;

		float pivot = this->Samples[SamplesIndices[(li + ri) / 2]].Position.Elements[Axis];
		for (;;)
		{
			while (this->Samples[SamplesIndices[i]].Position.Elements[Axis] < pivot)
			{
				++i;
			}

			while (this->Samples[SamplesIndices[j]].Position.Elements[Axis] > pivot)
			{
				--j;
			}

			if (i >= j) break;

			int temp = SamplesIndices[i];
			SamplesIndices[i] = SamplesIndices[j];
			SamplesIndices[j] = temp;

			++i;
			--j;
		}
		
		if (li < (i - 1)) 
		{
			this->SortSamples(SamplesIndices, Axis, li, i - 1);
		}
		if ((j + 1) <  ri) 
		{
			this->SortSamples(SamplesIndices, Axis, j + 1, ri);
		}
	}
}



template <int DIMENSION>
void CKdTree<DIMENSION>::Build2D()
{
  //this->Nodes.clear();
  //this->Samples.clear();
  this->CurrentNum2DNodes = 0;

	int SamplesNum = this->CurrentNumSamples;
  std::vector<int> SamplesIndices(SamplesNum);
  for (int i = 0; i < SamplesNum; ++i)
  {
    SamplesIndices[i] = i;
  }

	CBBox<DIMENSION> BBox;
	for (int i = 0; i < DIMENSION; i++)
	{
		BBox.Min.Elements[i] = 0.0f;
		BBox.Max.Elements[i] = this->SampleExtent.Elements[i];
	}

  this->MakeNode2D(SamplesIndices, 0, BBox);
}






template <int DIMENSION>
void CKdTree<DIMENSION>::Build(const std::vector<CSample<DIMENSION>> SamplesList)
{
  //this->Nodes.clear();
  //this->Samples.clear();
  this->CurrentNumNodes = 0;
  this->CurrentNumSamples = 0;

	int SamplesNum = (int)SamplesList.size();
  std::vector<int> SamplesIndices(SamplesNum);
  for (int i = 0; i < SamplesNum; ++i)
  {
    this->Samples[i] = SamplesList[i];
    SamplesIndices[i] = i;
  }
  this->CurrentNumSamples = SamplesNum;

	CBBox<DIMENSION> BBox;
	for (int i = 0; i < DIMENSION; i++)
	{
		BBox.Min.Elements[i] = 0.0f;
		BBox.Max.Elements[i] = this->SampleExtent.Elements[i];
	}

  this->MakeNode(SamplesIndices, 0, BBox);
//	this->Nodes[0].Parent = -1;



/*

  this->CurrentNum2DNodes = 0;
  for (int i = 0; i < SamplesNum; ++i)
  {
    SamplesIndices[i] = i;
  }

	for (int i = 0; i < DIMENSION; i++)
	{
		BBox.Min.Elements[i] = 0.0f;
		BBox.Max.Elements[i] = this->SampleExtent.Elements[i];
	}

  this->MakeNode2D(SamplesIndices, 0, BBox);
*/
}







template <int DIMENSION>
int CKdTree<DIMENSION>::MakeNodeRecon(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox)
{
  int CurrentIndex = this->CurrentNumNodes;
  this->CurrentNumNodes++;

	int SampleIndicesNum = (int)SamplesIndices.size();

  if (SampleIndicesNum <= this->BucketSize)
  {
    // we have a space in the bucket
		this->Nodes[CurrentIndex].SamplesIndices.reserve(this->BucketSize);
    this->Nodes[CurrentIndex].NodeType = NodeTypeLeaf;
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      this->Nodes[CurrentIndex].SamplesIndices.push_back(SamplesIndices[k]);
    }
    this->Nodes[CurrentIndex].Axis = Axis; 
		this->Nodes[CurrentIndex].BBox = BBox;
  }
  else
  {
		// set bounding box
		this->Nodes[CurrentIndex].BBox = BBox;
		int NewAxis = BBox.GetLargestAxisID();

    // we need to split
    this->SortSamples(SamplesIndices, NewAxis, 0, SampleIndicesNum - 1);
    int Median = SampleIndicesNum / 2;

    // split by the median
    std::vector<int> LeftSamplesIndices;
    std::vector<int> RightSamplesIndices;
		LeftSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
		RightSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      if (k < Median)
      {
        LeftSamplesIndices.push_back(SamplesIndices[k]);
      }
      else
      {
        RightSamplesIndices.push_back(SamplesIndices[k]);
      }
    }

    // this node is a splitting node
    this->Nodes[CurrentIndex].NodeType = NodeTypeSplit;

		// store axis and threshold
    this->Nodes[CurrentIndex].Axis = NewAxis;
    this->Nodes[CurrentIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[NewAxis] + this->Samples[SamplesIndices[Median]].Position.Elements[NewAxis]);

		CBBox<DIMENSION> LeftBBox;
		CBBox<DIMENSION> RightBBox;
		LeftBBox = BBox;
		RightBBox = BBox;
		LeftBBox.Max.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;
		RightBBox.Min.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;

		// make child nodes
    this->Nodes[CurrentIndex].LeftChild = this->MakeNodeRecon(LeftSamplesIndices, (Axis + 1) % DIMENSION, LeftBBox);
    this->Nodes[CurrentIndex].RightChild = this->MakeNodeRecon(RightSamplesIndices, (Axis + 1) % DIMENSION, RightBBox);

		// store the index to the parent node
//		this->Nodes[this->Nodes[CurrentIndex].LeftChild].Parent = CurrentIndex;
//		this->Nodes[this->Nodes[CurrentIndex].RightChild].Parent = CurrentIndex;
	}

  return CurrentIndex;
}





template <int DIMENSION>
int CKdTree<DIMENSION>::MakeNode(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox)
{
  int CurrentIndex = this->CurrentNumNodes;
  this->CurrentNumNodes++;

	int SampleIndicesNum = (int)SamplesIndices.size();

  if (SampleIndicesNum <= this->BucketSize)
  {
    // we have a space in the bucket
		this->Nodes[CurrentIndex].SamplesIndices.reserve(this->BucketSize);
    this->Nodes[CurrentIndex].NodeType = NodeTypeLeaf;
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      this->Nodes[CurrentIndex].SamplesIndices.push_back(SamplesIndices[k]);
    }
    this->Nodes[CurrentIndex].Axis = Axis; 
		this->Nodes[CurrentIndex].BBox = BBox;
  }
  else
  {
		// set bounding box
		this->Nodes[CurrentIndex].BBox = BBox;
		int NewAxis = BBox.GetLargestAxis();

    // we need to split
    this->SortSamples(SamplesIndices, NewAxis, 0, SampleIndicesNum - 1);
    int Median = SampleIndicesNum / 2;

    // split by the median
    std::vector<int> LeftSamplesIndices;
    std::vector<int> RightSamplesIndices;
		LeftSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
		RightSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      if (k < Median)
      {
        LeftSamplesIndices.push_back(SamplesIndices[k]);
      }
      else
      {
        RightSamplesIndices.push_back(SamplesIndices[k]);
      }
    }

    // this node is a splitting node
    this->Nodes[CurrentIndex].NodeType = NodeTypeSplit;

		// store axis and threshold
    this->Nodes[CurrentIndex].Axis = NewAxis;
    this->Nodes[CurrentIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[NewAxis] + this->Samples[SamplesIndices[Median]].Position.Elements[NewAxis]);

		CBBox<DIMENSION> LeftBBox;
		CBBox<DIMENSION> RightBBox;
		LeftBBox = BBox;
		RightBBox = BBox;
		LeftBBox.Max.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;
		RightBBox.Min.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;

		// make child nodes
    this->Nodes[CurrentIndex].LeftChild = this->MakeNode(LeftSamplesIndices, (Axis + 1) % DIMENSION, LeftBBox);
    this->Nodes[CurrentIndex].RightChild = this->MakeNode(RightSamplesIndices, (Axis + 1) % DIMENSION, RightBBox);

		// store the index to the parent node
//		this->Nodes[this->Nodes[CurrentIndex].LeftChild].Parent = CurrentIndex;
//		this->Nodes[this->Nodes[CurrentIndex].RightChild].Parent = CurrentIndex;
	}

  return CurrentIndex;
}






template <int DIMENSION>
int CKdTree<DIMENSION>::MakeNode2D(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox)
{
  int CurrentIndex = this->CurrentNum2DNodes;
  this->CurrentNum2DNodes++;

	int SampleIndicesNum = (int)SamplesIndices.size();

  if (SampleIndicesNum <= this->BucketSize)
  {
    // we have a space in the bucket
		this->Nodes2D[CurrentIndex].SamplesIndices.reserve(this->BucketSize);
    this->Nodes2D[CurrentIndex].NodeType = NodeTypeLeaf;
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      this->Nodes2D[CurrentIndex].SamplesIndices.push_back(SamplesIndices[k]);
    }
    this->Nodes2D[CurrentIndex].Axis = Axis; 
		this->Nodes2D[CurrentIndex].BBox = BBox;
  }
  else
  {
		// set bounding box
		this->Nodes2D[CurrentIndex].BBox = BBox;
		int NewAxis = BBox.GetLargestAxis2D();

    // we need to split
    this->SortSamples(SamplesIndices, NewAxis, 0, SampleIndicesNum - 1);
    int Median = SampleIndicesNum / 2;

    // split by the median
    std::vector<int> LeftSamplesIndices;
    std::vector<int> RightSamplesIndices;
		LeftSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
		RightSamplesIndices.reserve(SampleIndicesNum / 2 + 1);
    for (int k = 0; k < SampleIndicesNum; k++)
    {
      if (k < Median)
      {
        LeftSamplesIndices.push_back(SamplesIndices[k]);
      }
      else
      {
        RightSamplesIndices.push_back(SamplesIndices[k]);
      }
    }

    // this node is a splitting node
    this->Nodes2D[CurrentIndex].NodeType = NodeTypeSplit;

		// store axis and threshold
    this->Nodes2D[CurrentIndex].Axis = NewAxis;
    this->Nodes2D[CurrentIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[NewAxis] + this->Samples[SamplesIndices[Median]].Position.Elements[NewAxis]);

		CBBox<DIMENSION> LeftBBox;
		CBBox<DIMENSION> RightBBox;
		LeftBBox = BBox;
		RightBBox = BBox;
		LeftBBox.Max.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;
		RightBBox.Min.Elements[NewAxis] = this->Nodes[CurrentIndex].Threshold;

		// make child nodes
    this->Nodes2D[CurrentIndex].LeftChild = this->MakeNode2D(LeftSamplesIndices, (Axis + 1) % DIMENSION, LeftBBox);
    this->Nodes2D[CurrentIndex].RightChild = this->MakeNode2D(RightSamplesIndices, (Axis + 1) % DIMENSION, RightBBox);
	}

  return CurrentIndex;
}


 

template <int DIMENSION>
void CKdTree<DIMENSION>::Insert(CSample<DIMENSION>& Sample, bool Is2D)
{
	// add a new sample
	//this->Samples.push_back(Sample);
  this->Samples[this->CurrentNumSamples] = Sample;
  this->CurrentNumSamples++;
	int NewSampleIndex = this->CurrentNumSamples - 1;


#ifndef CELLSPLIT
	// do KNNs query for a new sample
	this->KNNSearch(this->KNNsNum, Sample.Position, &(this->Samples[NewSampleIndex].KNNs));

	// computes sample's bounding sphere
	this->Samples[NewSampleIndex].BSphere.Center = Sample.Position;
	this->Samples[NewSampleIndex].BSphere.SquaredRadius = this->Samples[NewSampleIndex].KNNs.MaxSquaredDistance;

	// compute its gradient
	//this->ComputeSampleGradient(NewSampleIndex);


	// get affected samples
	CKNNs ReverseKNNs;
	this->ReverseKNNSearch(Sample.Position, &ReverseKNNs);

	// update KNNs and gradients
	// bounding spheres always shrink, so the hierarchy does not need to be updated
	int ReverseKNNsNum = (int)ReverseKNNs.SamplesIndices.size();
	for (int i = 0; i < ReverseKNNsNum; ++i)
	{
		int ReverseKNNIndex = ReverseKNNs.SamplesIndices[i];
		this->Samples[ReverseKNNIndex].KNNs.Updated = false;
	}
	for (int i = 0; i < ReverseKNNsNum; ++i)
	{
		int ReverseKNNIndex = ReverseKNNs.SamplesIndices[i];
		this->SearchTest(NewSampleIndex, this->Samples[ReverseKNNIndex].Position, &(this->Samples[ReverseKNNIndex].KNNs));
	}
	for (int i = 0; i < ReverseKNNsNum; ++i)
	{
		int ReverseKNNIndex = ReverseKNNs.SamplesIndices[i];
    if (this->Samples[ReverseKNNIndex].KNNs.Updated)
    {
      this->ComputeSampleGradient(ReverseKNNIndex);
		  this->Samples[ReverseKNNIndex].BSphere.SquaredRadius = this->Samples[ReverseKNNIndex].KNNs.MaxSquaredDistance;
    }
	}
#endif

  // find a leaf node that contains a new sample
  int NodeIndex = 0;
  while (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
    if (Sample.Position.Elements[this->Nodes[NodeIndex].Axis] < this->Nodes[NodeIndex].Threshold)
    {
      NodeIndex = this->Nodes[NodeIndex].LeftChild;
    }
    else
    {
      NodeIndex = this->Nodes[NodeIndex].RightChild;
    }
  }



  if ((int)this->Nodes[NodeIndex].SamplesIndices.size() == this->BucketSize)
  {
    // the leaf node is already full, split 

		// gather indices
    std::vector<int> SamplesIndices;
		SamplesIndices.clear();
    for (int k = 0; k < this->BucketSize; k++)
    {
      SamplesIndices.push_back(this->Nodes[NodeIndex].SamplesIndices[k]);
    }
    SamplesIndices.push_back(NewSampleIndex);


		// split into two leaves
		int SampleIndicesNum = this->BucketSize + 1;
		int LeftChildNode = this->CurrentNumNodes;
		int RightChildNode = this->CurrentNumNodes + 1;
		//this->Nodes.resize(this->Nodes.size() + 2);
    this->CurrentNumNodes += 2;

    // split by the median
		int Axis = this->Nodes[NodeIndex].Axis;
		Axis = this->Nodes[NodeIndex].BBox.GetLargestAxis();
		//if (Is2D)
		//{
		//	Axis = this->Nodes[NodeIndex].BBox.GetLargestAxis2D();
		//}
    this->Nodes[NodeIndex].Axis = Axis;

    int Median = SampleIndicesNum / 2;
		this->SortSamples(SamplesIndices, Axis, 0, SampleIndicesNum - 1);

		//// for BucketNum = 4, explicit sort
		//if (this->Samples[SamplesIndices[0]].Position.Elements[Axis] > this->Samples[SamplesIndices[1]].Position.Elements[Axis]) std::swap(SamplesIndices[0], SamplesIndices[1]);
		//if (this->Samples[SamplesIndices[1]].Position.Elements[Axis] > this->Samples[SamplesIndices[2]].Position.Elements[Axis]) std::swap(SamplesIndices[1], SamplesIndices[2]);
		//if (this->Samples[SamplesIndices[2]].Position.Elements[Axis] > this->Samples[SamplesIndices[3]].Position.Elements[Axis]) std::swap(SamplesIndices[2], SamplesIndices[3]);
		//if (this->Samples[SamplesIndices[3]].Position.Elements[Axis] > this->Samples[SamplesIndices[4]].Position.Elements[Axis]) std::swap(SamplesIndices[3], SamplesIndices[4]);
		//if (this->Samples[SamplesIndices[0]].Position.Elements[Axis] > this->Samples[SamplesIndices[1]].Position.Elements[Axis]) std::swap(SamplesIndices[0], SamplesIndices[1]);
		//if (this->Samples[SamplesIndices[1]].Position.Elements[Axis] > this->Samples[SamplesIndices[2]].Position.Elements[Axis]) std::swap(SamplesIndices[1], SamplesIndices[2]);
		//if (this->Samples[SamplesIndices[2]].Position.Elements[Axis] > this->Samples[SamplesIndices[3]].Position.Elements[Axis]) std::swap(SamplesIndices[2], SamplesIndices[3]);
		//if (this->Samples[SamplesIndices[0]].Position.Elements[Axis] > this->Samples[SamplesIndices[1]].Position.Elements[Axis]) std::swap(SamplesIndices[0], SamplesIndices[1]);
		//if (this->Samples[SamplesIndices[1]].Position.Elements[Axis] > this->Samples[SamplesIndices[2]].Position.Elements[Axis]) std::swap(SamplesIndices[1], SamplesIndices[2]);
		//if (this->Samples[SamplesIndices[0]].Position.Elements[Axis] > this->Samples[SamplesIndices[1]].Position.Elements[Axis]) std::swap(SamplesIndices[0], SamplesIndices[1]);


		this->Nodes[LeftChildNode].SamplesIndices.clear();
		this->Nodes[RightChildNode].SamplesIndices.clear();
		this->Nodes[LeftChildNode].SamplesIndices.reserve(this->BucketSize);
		this->Nodes[RightChildNode].SamplesIndices.reserve(this->BucketSize);
    for (int k = 0; k < SampleIndicesNum; ++k)
    {
      if (k < Median)
      {
        this->Nodes[LeftChildNode].SamplesIndices.push_back(SamplesIndices[k]);
      }
      else
      {
        this->Nodes[RightChildNode].SamplesIndices.push_back(SamplesIndices[k]);
      }
    }

    // this node is now a splitting node
    this->Nodes[NodeIndex].NodeType = NodeTypeSplit;

		// store the threshold
    this->Nodes[NodeIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[Axis] + this->Samples[SamplesIndices[Median]].Position.Elements[Axis]);
    //this->Nodes[NodeIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[Axis] + this->Samples[SamplesIndices[Median]].Position.Elements[Axis]);

		// make child nodes
    this->Nodes[NodeIndex].LeftChild = LeftChildNode;
    this->Nodes[NodeIndex].RightChild = RightChildNode;


		CBBox<DIMENSION> LeftBBox;
		CBBox<DIMENSION> RightBBox;
		LeftBBox = this->Nodes[NodeIndex].BBox;
		RightBBox = this->Nodes[NodeIndex].BBox;
		LeftBBox.Max.Elements[Axis] = this->Nodes[NodeIndex].Threshold;
		RightBBox.Min.Elements[Axis] = this->Nodes[NodeIndex].Threshold;



		//// make a new subtree at the end of the data
    //int NewNodeIndex = this->MakeNode(SamplesIndices, this->Nodes[NodeIndex].Axis);
		this->Nodes[LeftChildNode].NodeType = NodeTypeLeaf;
    this->Nodes[LeftChildNode].Axis = (Axis + 1) % DIMENSION; 
//		this->Nodes[LeftChildNode].Parent = NodeIndex;
		this->Nodes[LeftChildNode].BBox = LeftBBox;

		this->Nodes[RightChildNode].NodeType = NodeTypeLeaf;
    this->Nodes[RightChildNode].Axis = (Axis + 1) % DIMENSION; 
//		this->Nodes[RightChildNode].Parent = NodeIndex;
		this->Nodes[RightChildNode].BBox = RightBBox;

#ifndef CELLSPLIT
		// compute bounding spheres for child nodes
		CVector ZeroVector;
		for (int i = 0; i < DIMENSION; ++i)
		{
			ZeroVector.Elements[i] = 0.0f;
		}

		this->Nodes[LeftChildNode].BSphere.Center = ZeroVector;
		int LeftChildeNodeSamplesNum = (int)(this->Nodes[LeftChildNode].SamplesIndices.size());
		for (int i = 0; i < LeftChildeNodeSamplesNum; i++)
		{
			this->Nodes[LeftChildNode].BSphere.Center = this->Nodes[LeftChildNode].BSphere.Center + this->Samples[this->Nodes[LeftChildNode].SamplesIndices[i]].BSphere.Center;
		}
		this->Nodes[LeftChildNode].BSphere.Center = this->Nodes[LeftChildNode].BSphere.Center * (1.0f / (float)LeftChildeNodeSamplesNum);
		this->Nodes[LeftChildNode].BSphere.SquaredRadius = 0.0f;
		for (int i = 0; i < LeftChildeNodeSamplesNum; i++)
		{
			this->Nodes[LeftChildNode].BSphere.Combine(this->Samples[this->Nodes[LeftChildNode].SamplesIndices[i]].BSphere);
		}
		
		this->Nodes[RightChildNode].BSphere.Center = ZeroVector;
		int RightChildeNodeSamplesNum = (int)(this->Nodes[RightChildNode].SamplesIndices.size());
		for (int i = 0; i < RightChildeNodeSamplesNum; i++)
		{
			this->Nodes[RightChildNode].BSphere.Center = this->Nodes[RightChildNode].BSphere.Center + this->Samples[this->Nodes[RightChildNode].SamplesIndices[i]].BSphere.Center;
		}
		this->Nodes[RightChildNode].BSphere.Center = this->Nodes[RightChildNode].BSphere.Center * (1.0f / (float)RightChildeNodeSamplesNum);
		this->Nodes[RightChildNode].BSphere.SquaredRadius = 0.0f;
		for (int i = 0; i < RightChildeNodeSamplesNum; i++)
		{
			this->Nodes[RightChildNode].BSphere.Combine(this->Samples[this->Nodes[RightChildNode].SamplesIndices[i]].BSphere);
		}

		// compute a current bounding sphere
    this->Nodes[NodeIndex].BSphere.SquaredRadius = 0.0f;
		this->Nodes[NodeIndex].BSphere.Center = 0.5f * (this->Nodes[LeftChildNode].BSphere.Center + this->Nodes[RightChildNode].BSphere.Center);
		this->Nodes[NodeIndex].BSphere.Combine(this->Nodes[LeftChildNode].BSphere);
		this->Nodes[NodeIndex].BSphere.Combine(this->Nodes[RightChildNode].BSphere);

		// propagate the change of bounding sphere toward parent
		while (this->Nodes[NodeIndex].Parent != -1)
		{	
			this->Nodes[this->Nodes[NodeIndex].Parent].BSphere.Combine(this->Nodes[NodeIndex].BSphere);
			NodeIndex = this->Nodes[NodeIndex].Parent;
		}
#else
    this->Heap.Update(NodeIndex, -1.0f);
    //this->Heap.Remove(NodeIndex);

    this->ComputeNodeError(LeftChildNode, Is2D);
    this->Heap.Put(LeftChildNode, this->Nodes[LeftChildNode].Error);

    this->ComputeNodeError(RightChildNode, Is2D);
    this->Heap.Put(RightChildNode, this->Nodes[RightChildNode].Error);
#endif
  }
  else
  {
    // we have enough space, simply add a new sample
    this->Nodes[NodeIndex].SamplesIndices.push_back(NewSampleIndex);

#ifndef CELLSPLIT
		// update bounding sphere
		this->Nodes[NodeIndex].BSphere.Combine(this->Samples[NewSampleIndex].BSphere);

		// propagate the change of bounding sphere toward parent
		while (this->Nodes[NodeIndex].Parent != -1)
		{	
			this->Nodes[this->Nodes[NodeIndex].Parent].BSphere.Combine(this->Nodes[NodeIndex].BSphere);
			NodeIndex = this->Nodes[NodeIndex].Parent;
		}
#else
    this->ComputeNodeError(NodeIndex, Is2D);
    this->Heap.Update(NodeIndex, this->Nodes[NodeIndex].Error);
#endif

  }

#ifndef CELLSPLIT
  // update distances
	this->ComputeSampleDistance(NewSampleIndex);
  this->Heap.Put(NewSampleIndex, this->Samples[NewSampleIndex].Distance); 
  for (int i = 0; i < ReverseKNNsNum; ++i)
	{
		int ReverseKNNIndex = ReverseKNNs.SamplesIndices[i];
    if (this->Samples[ReverseKNNIndex].KNNs.Updated)
    {
		  this->ComputeSampleDistance(ReverseKNNIndex);
      this->Heap.Update(ReverseKNNIndex, this->Samples[ReverseKNNIndex].Distance);
    }
	}
#endif












/*

  // find a leaf node that contains a new sample
  NodeIndex = 0;
  while (this->Nodes2D[NodeIndex].NodeType == NodeTypeSplit)
  {
    if (Sample.Position.Elements[this->Nodes2D[NodeIndex].Axis] < this->Nodes2D[NodeIndex].Threshold)
    {
      NodeIndex = this->Nodes2D[NodeIndex].LeftChild;
    }
    else
    {
      NodeIndex = this->Nodes2D[NodeIndex].RightChild;
    }
  }


  if ((int)this->Nodes2D[NodeIndex].SamplesIndices.size() == this->BucketSize)
  {
    // the leaf node is already full, split 

		// gather indices
    std::vector<int> SamplesIndices;
		SamplesIndices.clear();
    for (int k = 0; k < this->BucketSize; k++)
    {
      SamplesIndices.push_back(this->Nodes2D[NodeIndex].SamplesIndices[k]);
    }
    SamplesIndices.push_back(NewSampleIndex);


		// split into two leaves
		int SampleIndicesNum = this->BucketSize + 1;
		int LeftChildNode = this->CurrentNum2DNodes;
		int RightChildNode = this->CurrentNum2DNodes + 1;
    this->CurrentNum2DNodes += 2;

    // split by the median
		int Axis = this->Nodes2D[NodeIndex].Axis;
		Axis = this->Nodes2D[NodeIndex].BBox.GetLargestAxis2D();
    this->Nodes2D[NodeIndex].Axis = Axis;

    int Median = SampleIndicesNum / 2;
		this->SortSamples(SamplesIndices, Axis, 0, SampleIndicesNum - 1);

		this->Nodes2D[LeftChildNode].SamplesIndices.clear();
		this->Nodes2D[RightChildNode].SamplesIndices.clear();
		this->Nodes2D[LeftChildNode].SamplesIndices.reserve(this->BucketSize);
		this->Nodes2D[RightChildNode].SamplesIndices.reserve(this->BucketSize);
    for (int k = 0; k < SampleIndicesNum; ++k)
    {
      if (k < Median)
      {
        this->Nodes2D[LeftChildNode].SamplesIndices.push_back(SamplesIndices[k]);
      }
      else
      {
        this->Nodes2D[RightChildNode].SamplesIndices.push_back(SamplesIndices[k]);
      }
    }

    // this node is now a splitting node
    this->Nodes2D[NodeIndex].NodeType = NodeTypeSplit;

		// store the threshold
    this->Nodes2D[NodeIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[Axis] + this->Samples[SamplesIndices[Median]].Position.Elements[Axis]);
    //this->Nodes[NodeIndex].Threshold = 0.5f * (this->Samples[SamplesIndices[Median - 1]].Position.Elements[Axis] + this->Samples[SamplesIndices[Median]].Position.Elements[Axis]);

		// make child nodes
    this->Nodes2D[NodeIndex].LeftChild = LeftChildNode;
    this->Nodes2D[NodeIndex].RightChild = RightChildNode;


		CBBox LeftBBox;
		CBBox RightBBox;
		LeftBBox = this->Nodes2D[NodeIndex].BBox;
		RightBBox = this->Nodes2D[NodeIndex].BBox;
		LeftBBox.Max.Elements[Axis] = this->Nodes2D[NodeIndex].Threshold;
		RightBBox.Min.Elements[Axis] = this->Nodes2D[NodeIndex].Threshold;



		//// make a new subtree at the end of the data
		this->Nodes2D[LeftChildNode].NodeType = NodeTypeLeaf;
    this->Nodes2D[LeftChildNode].Axis = (Axis + 1) % DIMENSION; 
		this->Nodes2D[LeftChildNode].BBox = LeftBBox;

		this->Nodes2D[RightChildNode].NodeType = NodeTypeLeaf;
    this->Nodes2D[RightChildNode].Axis = (Axis + 1) % DIMENSION; 
		this->Nodes2D[RightChildNode].BBox = RightBBox;

  }
  else
  {
    // we have enough space, simply add a new sample
    this->Nodes2D[NodeIndex].SamplesIndices.push_back(NewSampleIndex);
  }

*/


}




#ifndef CELLSPLIT
CBSphere CKdTree::ComputeBSpheresIteration(const int NodeIndex)
{
	if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
	{
		for (int i = 0; i < DIMENSION; i++)
		{
			this->Nodes[NodeIndex].BSphere.Center.Elements[i] = 0.0f;
		}

		int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
		for (int i = 0; i < LeafNodeSamplesNum; i++)
		{
			this->Nodes[NodeIndex].BSphere.Center = this->Nodes[NodeIndex].BSphere.Center + this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].BSphere.Center;
		}
		this->Nodes[NodeIndex].BSphere.Center = this->Nodes[NodeIndex].BSphere.Center * (1.0f / (float)LeafNodeSamplesNum);

		this->Nodes[NodeIndex].BSphere.SquaredRadius = 0.0f;
		for (int i = 0; i < LeafNodeSamplesNum; i++)
		{
			this->Nodes[NodeIndex].BSphere.Combine(this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].BSphere);
		}
	}
	else
	{
		CBSphere LeftBSphere;
		CBSphere RightBSphere;

		LeftBSphere = this->ComputeBSpheresIteration(this->Nodes[NodeIndex].LeftChild);
		RightBSphere = this->ComputeBSpheresIteration(this->Nodes[NodeIndex].RightChild);

		// compute a current bounding sphere
		this->Nodes[NodeIndex].BSphere.SquaredRadius = 0.0f;
		this->Nodes[NodeIndex].BSphere.Center = 0.5f * (LeftBSphere.Center + RightBSphere.Center);
		this->Nodes[NodeIndex].BSphere.Combine(LeftBSphere);
		this->Nodes[NodeIndex].BSphere.Combine(RightBSphere);
	}

	return (this->Nodes[NodeIndex].BSphere);
}
#endif



#ifndef CELLSPLIT
void CKdTree::ComputeBSpheres()
{
	int SamplesNum = this->CurrentNumSamples;
	for (int n = 0; n < SamplesNum; n++)
	{
		// get sample's bounding sphere
		this->Samples[n].BSphere.Center = this->Samples[n].Position;
		this->Samples[n].BSphere.SquaredRadius = this->Samples[n].KNNs.MaxSquaredDistance;
	}

	// build BVH
	this->ComputeBSpheresIteration(0);
}
#endif



template <int DIMENSION>
void CKdTree<DIMENSION>::SearchTest(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs)
{
  float SquaredDistance = this->Samples[SampleIndex].Position.SquaredDistance(QueryPosition);
  if (SquaredDistance < KNNs->MaxSquaredDistance)
  {
    KNNs->Updated = true;

		KNNs->MaxSampleIndex = KNNs->Heap.Items[0].Index;
		KNNs->SamplesIndices[KNNs->MaxSampleIndex] = SampleIndex;
		KNNs->SamplesDistances[KNNs->MaxSampleIndex] = SquaredDistance;
		KNNs->Heap.Update(KNNs->MaxSampleIndex, SquaredDistance);

		KNNs->MaxSquaredDistance = KNNs->Heap.Items[0].Value;
	}
}

template <int DIMENSION>
void CKdTree<DIMENSION>::SearchTestAAEllipse(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const CVector<DIMENSION> & scale)
{
  // compute scaled difference vector
  CVector<DIMENSION> diff = this->Samples[SampleIndex].Position - QueryPosition;
  for (int i = 0; i < DIMENSION; ++i)
    diff.Elements[i] *= scale.Elements[i];

  float SquaredDistance = diff.DotProduct(diff);
  
  if (SquaredDistance < KNNs->MaxSquaredDistance)
  {
    KNNs->Updated = true;

    KNNs->MaxSampleIndex = KNNs->Heap.Items[0].Index;
    KNNs->SamplesIndices[KNNs->MaxSampleIndex] = SampleIndex;
    KNNs->SamplesDistances[KNNs->MaxSampleIndex] = SquaredDistance;
    KNNs->Heap.Update(KNNs->MaxSampleIndex, SquaredDistance);

    KNNs->MaxSquaredDistance = KNNs->Heap.Items[0].Value;
  }
}

template <int DIMENSION>
void CKdTree<DIMENSION>::SearchTestEllipse(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv)
{
  // compute warped distance between the query position and the sample position
  CVector<DIMENSION> diff = this->Samples[SampleIndex].Position - QueryPosition;
  CVector<DIMENSION> warpedDiff;
  for (int j = 0; j < DIMENSION; ++j)
  {
    for (int i = 0; i < DIMENSION; ++i)
    {
      warpedDiff.Elements[j] = M.Rows[j].Elements[i] * diff.Elements[i];
    }
  }
  float SquaredDistance = warpedDiff.DotProduct(warpedDiff);
  
  if (SquaredDistance < KNNs->MaxSquaredDistance)
  {
    KNNs->Updated = true;

    KNNs->MaxSampleIndex = KNNs->Heap.Items[0].Index;
    KNNs->SamplesIndices[KNNs->MaxSampleIndex] = SampleIndex;
    KNNs->SamplesDistances[KNNs->MaxSampleIndex] = SquaredDistance;
    KNNs->Heap.Update(KNNs->MaxSampleIndex, SquaredDistance);

    KNNs->MaxSquaredDistance = KNNs->Heap.Items[0].Value;
  }
}


template <int DIMENSION>
void CKdTree<DIMENSION>::SearchIteration(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, float DistanceToCell)
{
  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
  {
		int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
    for (int i = 0; i < LeafNodeSamplesNum; ++i)
    {
			int CurrentKNNsNum = (int)KNNs->SamplesIndices.size();

      if (CurrentKNNsNum == RequiredKNNsNum)
      {
        this->SearchTest(this->Nodes[NodeIndex].SamplesIndices[i], QueryPosition, KNNs);
      }
      else
      {
        int NewSampleIndex = this->Nodes[NodeIndex].SamplesIndices[i];
        KNNs->SamplesIndices.push_back(NewSampleIndex);
        KNNs->SamplesDistances.push_back(this->Samples[NewSampleIndex].Position.SquaredDistance(QueryPosition));
				KNNs->Heap.Put(CurrentKNNsNum, KNNs->SamplesDistances[CurrentKNNsNum]);
      }        
    }

  }

  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
    int Axis = this->Nodes[NodeIndex].Axis;
    float SignedDistance = QueryPosition.Elements[Axis] - this->Nodes[NodeIndex].Threshold;
    float SquaredDistance = SignedDistance * SignedDistance;

		float OldOffSet = this->OffSet.Elements[Axis];
		float NewOffSet = SignedDistance;
    
    if (SignedDistance < 0.0f)
    {  
      this->SearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
			DistanceToCell += - OldOffSet * OldOffSet + NewOffSet * NewOffSet;
//      if (SquaredDistance * 10.5f < KNNs->MaxSquaredDistance)
			if (DistanceToCell * (1.0f + 2.0f) < KNNs->MaxSquaredDistance)
      {
				this->OffSet.Elements[Axis] = NewOffSet;
        this->SearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
				this->OffSet.Elements[Axis] = OldOffSet;
      }
    }
    else
    {
      this->SearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
			DistanceToCell += - OldOffSet * OldOffSet + NewOffSet * NewOffSet;
//      if (SquaredDistance * 10.5f < KNNs->MaxSquaredDistance)
			if (DistanceToCell * (1.0f + 2.0f) < KNNs->MaxSquaredDistance)
      {
				this->OffSet.Elements[Axis] = NewOffSet;
        this->SearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
				this->OffSet.Elements[Axis] = OldOffSet;
      }
    }
  }

}






//void CKdTree::SearchIteration(const int NodeIndex, const CVector & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, float DistanceToCell)
//{
//  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
//  {
//		int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
//    for (int i = 0; i < LeafNodeSamplesNum; ++i)
//    {
//			int CurrentKNNsNum = (int)KNNs->SamplesIndices.size();
//
//      if (CurrentKNNsNum == RequiredKNNsNum)
//      {
//        this->SearchTest(this->Nodes[NodeIndex].SamplesIndices[i], QueryPosition, KNNs);
//      }
//      else
//      {
//        int NewSampleIndex = this->Nodes[NodeIndex].SamplesIndices[i];
//        KNNs->SamplesIndices.push_back(NewSampleIndex);
//        KNNs->SamplesDistances.push_back(this->Samples[NewSampleIndex].Position.SquaredDistance(QueryPosition));
//				KNNs->Heap.Put(CurrentKNNsNum, KNNs->SamplesDistances[CurrentKNNsNum]);
//      }        
//    }
//
//  }
//
//  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
//  {
//    int Axis = this->Nodes[NodeIndex].Axis;
//    float SignedDistance = QueryPosition.Elements[Axis] - this->Nodes[NodeIndex].Threshold;
//    float SquaredDistance = SignedDistance * SignedDistance;
//    
//    if (SignedDistance < 0.0f)
//    { 
//      this->SearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
//      if (SquaredDistance < KNNs->MaxSquaredDistance)
//      {
//        this->SearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
//      }
//    }
//    else
//    {
//      this->SearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
//      if (SquaredDistance < KNNs->MaxSquaredDistance)
//      {
//        this->SearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, DistanceToCell);
//      }
//    }
//  }
//
//}
//





template <int DIMENSION>
void CKdTree<DIMENSION>::SearchIterationAAEllipse(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, const CVector<DIMENSION> & scale)
{
  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
  {
    int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
    for (int i = 0; i < LeafNodeSamplesNum; ++i)
    {
      int CurrentKNNsNum = (int)KNNs->SamplesIndices.size();
      if (CurrentKNNsNum == RequiredKNNsNum)
      {
        this->SearchTestAAEllipse(this->Nodes[NodeIndex].SamplesIndices[i], QueryPosition, KNNs, scale);
      }
      else
      {
        int NewSampleIndex = this->Nodes[NodeIndex].SamplesIndices[i];
        KNNs->SamplesIndices.push_back(NewSampleIndex);

        // compute scaled difference vector
        CVector<DIMENSION> diff = this->Samples[NewSampleIndex].Position - QueryPosition;
        for (int i = 0; i < DIMENSION; ++i)
          diff.Elements[i] *= scale.Elements[i];
        
        KNNs->SamplesDistances.push_back(diff.DotProduct(diff));
        KNNs->Heap.Put(CurrentKNNsNum, KNNs->SamplesDistances[CurrentKNNsNum]);
      }        
    }

  }

  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
    int Axis = this->Nodes[NodeIndex].Axis;
    float SignedDistance = (QueryPosition.Elements[Axis] - this->Nodes[NodeIndex].Threshold) * scale.Elements[Axis];
    float SquaredDistance = SignedDistance * SignedDistance;
    
    if (SignedDistance < 0.0f)
    { 
      this->SearchIterationAAEllipse(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, scale);
      if (SquaredDistance < KNNs->MaxSquaredDistance)
      {
        this->SearchIterationAAEllipse(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, scale);
      }
    }
    else
    {
      this->SearchIterationAAEllipse(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, scale);
      if (SquaredDistance < KNNs->MaxSquaredDistance)
      {
        this->SearchIterationAAEllipse(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, scale);
      }
    }
  }
}


template <int DIMENSION>
void CKdTree<DIMENSION>::SearchIterationEllipse(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv)
{
  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
  {
    int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
    for (int i = 0; i < LeafNodeSamplesNum; ++i)
    {
      int CurrentKNNsNum = (int)KNNs->SamplesIndices.size();

      if (CurrentKNNsNum == RequiredKNNsNum)
      {
        this->SearchTestEllipse(this->Nodes[NodeIndex].SamplesIndices[i], QueryPosition, KNNs, M, MInv);
      }
      else
      {
        int NewSampleIndex = this->Nodes[NodeIndex].SamplesIndices[i];
        KNNs->SamplesIndices.push_back(NewSampleIndex);

        // compute warped distance between the query position and the sample position
        CVector<DIMENSION> diff = this->Samples[NewSampleIndex].Position - QueryPosition;
        CVector<DIMENSION> warpedDiff;
        for (int j = 0; j < DIMENSION; ++j)
        {
          for (int i = 0; i < DIMENSION; ++i)
          {
            warpedDiff.Elements[j] = M.Rows[j].Elements[i] * diff.Elements[i];
          }
        }
        
        KNNs->SamplesDistances.push_back(warpedDiff.DotProduct(warpedDiff));
        KNNs->Heap.Put(CurrentKNNsNum, KNNs->SamplesDistances[CurrentKNNsNum]);
      }
    }
  }
  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
    int Axis = this->Nodes[NodeIndex].Axis;
    float split = this->Nodes[NodeIndex].Threshold;

    //
    // Compute signed distance to warped splitting plane
    //
    // Splitting plane impicit form: point p, and normal n.
    //
    // CVector p;
    // p.Elements[Axis] = split;
    // CVector n;
    // p.Elements[Axis] = 1.0f;

    // warped version of splitting plane
    CVector<DIMENSION> PPrime;
    CVector<DIMENSION> NPrime;
    for (int i = 0; i < DIMENSION; ++i)
    {
      NPrime.Elements[i] = MInv.Rows[Axis].Elements[i];        // n' = (M^-1)^T * n
      PPrime.Elements[i] = M.Rows[i].Elements[Axis] * split;   // p' = M * p;
    }

    float SignedDistance = (QueryPosition-PPrime).DotProduct(NPrime);
    

    //
    // warp the plane using the inverse matrix
    //
    
    
    float SquaredDistance = SignedDistance * SignedDistance;
    
    if (SignedDistance < 0.0f)
    { 
      this->SearchIterationEllipse(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, M, MInv);
      if (SquaredDistance < KNNs->MaxSquaredDistance)
      {
        this->SearchIterationEllipse(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, M, MInv);
      }
    }
    else
    {
      this->SearchIterationEllipse(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs, RequiredKNNsNum, M, MInv);
      if (SquaredDistance < KNNs->MaxSquaredDistance)
      {
        this->SearchIterationEllipse(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs, RequiredKNNsNum, M, MInv);
      }
    }
  }

}


#ifndef CELLSPLIT

void CKdTree::ReverseSearchIteration(const int NodeIndex, const CVector & QueryPosition, CKNNs* KNNs)
{
  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
  {
		int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
    for (int i = 0; i < LeafNodeSamplesNum; ++i)
    {
			int NewSampleIndex = this->Nodes[NodeIndex].SamplesIndices[i];
			if (this->Samples[NewSampleIndex].BSphere.Check(QueryPosition))
			{
				KNNs->SamplesIndices.push_back(NewSampleIndex);
			}
    }
  }

  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
		if (this->Nodes[this->Nodes[NodeIndex].LeftChild].BSphere.Check(QueryPosition))
		{
			this->ReverseSearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, KNNs);
		}
		if (this->Nodes[this->Nodes[NodeIndex].RightChild].BSphere.Check(QueryPosition))
		{
			this->ReverseSearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, KNNs);
		}
  }
}
#endif

template <int DIMENSION>
void CKdTree<DIMENSION>::NearestSearchIteration(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, int& NearestIndex, float& MinSquareDistance)
{
  if (this->Nodes[NodeIndex].NodeType == NodeTypeLeaf)
  {
		int LeafNodeSamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();
    for (int i = 0; i < LeafNodeSamplesNum; ++i)
    {
			float SquaredDistance = QueryPosition.SquaredDistance(this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].Position);
			//float dx = QueryPosition.Elements[0] - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].Position.Elements[0];
			//float dy = QueryPosition.Elements[1] - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].Position.Elements[1];
			//float SquaredDistance = dx * dx + dy * dy;

			if (SquaredDistance < MinSquareDistance)
			{
				MinSquareDistance = SquaredDistance;
				NearestIndex = this->Nodes[NodeIndex].SamplesIndices[i];
			}
		}
  }

  else if (this->Nodes[NodeIndex].NodeType == NodeTypeSplit)
  {
    int Axis = this->Nodes[NodeIndex].Axis;
    float SignedDistance = QueryPosition.Elements[Axis] - this->Nodes[NodeIndex].Threshold;
    float SquaredDistance = SignedDistance * SignedDistance;

			if (SignedDistance < 0.0f)
			{ 
				this->NearestSearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, NearestIndex, MinSquareDistance);
				if (SquaredDistance < MinSquareDistance)
				{
					this->NearestSearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, NearestIndex, MinSquareDistance);
				}
			}
			else
			{
				this->NearestSearchIteration(this->Nodes[NodeIndex].RightChild, QueryPosition, NearestIndex, MinSquareDistance);
				if (SquaredDistance < MinSquareDistance)
				{
					this->NearestSearchIteration(this->Nodes[NodeIndex].LeftChild, QueryPosition, NearestIndex, MinSquareDistance);
				}
			}
  }

}


template <int DIMENSION>
int CKdTree<DIMENSION>::NearestSearch(const CVector<DIMENSION>& QueryPosition, float SearchRadius)
{
  float MinSquareDistance = SearchRadius;
	int NearestIndex = 0;

  this->NearestSearchIteration(0, QueryPosition, NearestIndex, MinSquareDistance);

	return NearestIndex;
}

template <int DIMENSION>
void CKdTree<DIMENSION>::KNNSearch(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, float SearchRadius)
{
  KNNs->SamplesDistances.clear();
  KNNs->SamplesDistances.reserve(RequiredKNNsNum);
  KNNs->MaxSquaredDistance = SearchRadius;

  KNNs->SamplesIndices.clear();
  KNNs->SamplesIndices.reserve(RequiredKNNsNum);
  KNNs->MaxSampleIndex = -1;

	KNNs->Heap.Items.resize(RequiredKNNsNum);
	KNNs->Heap.Indices.resize(RequiredKNNsNum);
  KNNs->Heap.CurrentNumItems = 0;

	for (int i = 0; i < DIMENSION; ++i)
	{
		this->OffSet.Elements[i] = 0.0f;
	}

  this->SearchIteration(0, QueryPosition, KNNs, RequiredKNNsNum);

  std::vector<CHeapItem>().swap(KNNs->Heap.Items);
  std::vector<int>().swap(KNNs->Heap.Indices);
  std::vector<float>().swap(KNNs->SamplesDistances);

}



//! Find K-nearest neighbors from QueryPosition using a distance metric warped by an axis-aligned scaling matrix.
/*!
    The vector scale contains the scaling parameters for each of the axes.
*/
template <int DIMENSION>
void CKdTree<DIMENSION>::KNNSearchAAEllipse(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, const CVector<DIMENSION> & scale, float SearchRadius)
{
  KNNs->SamplesDistances.clear();
  KNNs->SamplesDistances.reserve(RequiredKNNsNum);
  KNNs->MaxSquaredDistance = SearchRadius;

  KNNs->SamplesIndices.clear();
  KNNs->SamplesIndices.reserve(RequiredKNNsNum);
  KNNs->MaxSampleIndex = -1;
  
  KNNs->Heap.Items.resize(RequiredKNNsNum);
  KNNs->Heap.Indices.resize(RequiredKNNsNum);
  KNNs->Heap.CurrentNumItems = 0;

  this->SearchIterationAAEllipse(0, QueryPosition, KNNs, RequiredKNNsNum, scale);
}



//! Find K-nearest neighbors from QueryPosition using a distance metric warped by a matrix.
/*!
    The columns of the warping matrix must be orthogonal to each other, but are
    allowed to have non-unit length. This means the warp is a non-uniform scale
    followed by a rotation.

    M contains the warping matrix: distance(a, b) = (b-a)^T * M * (b-a)
    MInv contains its inverse.
*/
template <int DIMENSION>
void CKdTree<DIMENSION>::KNNSearchEllipse(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv, float SearchRadius)
{
  KNNs->SamplesDistances.clear();
  KNNs->SamplesDistances.reserve(RequiredKNNsNum);
  KNNs->MaxSquaredDistance = SearchRadius;

  KNNs->SamplesIndices.clear();
  KNNs->SamplesIndices.reserve(RequiredKNNsNum);
  KNNs->MaxSampleIndex = -1;

  KNNs->Heap.Items.resize(RequiredKNNsNum);
  KNNs->Heap.Indices.resize(RequiredKNNsNum);
  KNNs->Heap.CurrentNumItems = 0;

  this->SearchIterationEllipse(0, QueryPosition, KNNs, RequiredKNNsNum, M, MInv);
}


#ifndef CELLSPLIT

void CKdTree::ReverseKNNSearch(const CVector& QueryPosition, CKNNs* KNNs)
{
  KNNs->SamplesIndices.clear();

  this->ReverseSearchIteration(0, QueryPosition, KNNs);
}
#endif




template <int DIMENSION>
void CKdTree<DIMENSION>::ComputeNodeError(const int NodeIndex, bool Is2D)
{  
	int SamplesNum = (int)this->Nodes[NodeIndex].SamplesIndices.size();

  float AverageLuminance[3];
  AverageLuminance[0] = 0.0f;
  AverageLuminance[1] = 0.0f;
  AverageLuminance[2] = 0.0f;
  for (int i = 0; i < SamplesNum; ++i)
  {
    for (int c = 0; c < 3; c++)
    {
    AverageLuminance[c] += fabs(this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].XYZ[c]);
    }
  }
  AverageLuminance[0] /= (float)(SamplesNum);
  AverageLuminance[1] /= (float)(SamplesNum);
  AverageLuminance[2] /= (float)(SamplesNum);

  this->Nodes[NodeIndex].Error = 0.0f;
  for (int i = 0; i < SamplesNum; ++i)
  {
    for (int c = 0; c < 3; c++)
    {
      float Difference = fabs((AverageLuminance[c] - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].XYZ[c] / (float)(SamplesNum)) * ((float)(SamplesNum) / (float)(SamplesNum - 1)) - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].XYZ[c]) / (AverageLuminance[c] + 1e-10f);
      //float Difference = AverageLuminance[c] - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].XYZ[c];
      this->Nodes[NodeIndex].Error += Difference;// * Difference;
    }
  }
  this->Nodes[NodeIndex].Error /= (float)(SamplesNum);

  //if (this->Nodes[NodeIndex].Error > 0.5f)
  //{
  //  this->Nodes[NodeIndex].Error = 0.5f;
  //}

  this->Nodes[NodeIndex].Error += this->DistanceEpsilon;

  //this->Nodes[NodeIndex].Error = 1E-5f;
  //for (int i = 0; i < SamplesNum; ++i)
  //{
  //  for (int j = 0; j < SamplesNum; ++j)
  //  {
  //    float Difference = this->Samples[this->Nodes[NodeIndex].SamplesIndices[j]].Luminance - this->Samples[this->Nodes[NodeIndex].SamplesIndices[i]].Luminance;
  //    this->Nodes[NodeIndex].Error += Difference * Difference;
  //  }
  //}
  //this->Nodes[NodeIndex].Error /= (float)(SamplesNum * SamplesNum);
  //if (this->Nodes[NodeIndex].Error > 0.05f) this->Nodes[NodeIndex].Error = 0.05f;

  this->Nodes[NodeIndex].Error *= this->Nodes[NodeIndex].BBox.ComputeVolume();
}




template <int DIMENSION>
void CKdTree<DIMENSION>::ComputeErrors()
{
  int NodesNum = this->CurrentNumNodes;
  for (int i = 0; i < NodesNum; ++i)
  {
    if (this->Nodes[i].NodeType != NodeTypeLeaf)
    {
      continue;
    }

    this->ComputeNodeError(i, false);
  }
}




template <int DIMENSION>
void CKdTree<DIMENSION>::ComputeSampleGradientKNNs(const int SampleIndex, const CKNNs &SKNNs)
{
  CVector<DIMENSION> QueryPosition = this->Samples[SampleIndex].Position;
  float QueryLuminance = this->Samples[SampleIndex].XYZ[1];
  CKNNs KNNs;

  for (int i = 0; i < (this->KNNsNum + 1); i++)
  {
    if (SampleIndex == SKNNs.SamplesIndices[i])
    {
      continue;
    }
    KNNs.SamplesIndices.push_back(SKNNs.SamplesIndices[i]);
  }


  // compute differences
  std::vector<CVector<DIMENSION>> PositionDifferences(this->KNNsNum);
  std::vector<float> LuminaceDifferences(this->KNNsNum);

  for (int i = 0; i < this->KNNsNum; ++i)
  {
		int KNNIndex = KNNs.SamplesIndices[i];
    PositionDifferences[i] = this->Samples[KNNIndex].Position - QueryPosition;
    LuminaceDifferences[i] = this->Samples[KNNIndex].XYZ[1] - QueryLuminance;
  }
 


  // construct a matrix and a vector for the least squares method
  float A[DIMENSION * DIMENSION];
  float b[DIMENSION];

  for (int j = 0; j < DIMENSION; ++j)
  {
		for (int i = 0; i < DIMENSION; ++i)
		{
      A[i + j * DIMENSION] = 0.0f;
      for (int k = 0; k < this->KNNsNum; ++k)
      {
        A[i + j * DIMENSION] = A[i + j * DIMENSION] + PositionDifferences[k].Elements[j] * PositionDifferences[k].Elements[i];
      }
    }
  }

  for (int j = 0; j < DIMENSION; ++j)
  {
    b[j] = 0.0f;
    for (int k = 0; k < this->KNNsNum; ++k)
    {
      b[j] = b[j] + PositionDifferences[k].Elements[j] * LuminaceDifferences[k];
    }
  }



  // solve by LAPACK
  long int N = DIMENSION;
  long int nrhs = 1;
  long int lda = DIMENSION;
  long int ipiv[DIMENSION];
  long int ldb = DIMENSION;
  long int info;
  sgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);



	// limit the length of the gradient
  float blen = 0.0f;
  for (int j = 0; j < DIMENSION; ++j)
  {
    blen = blen + b[j] * b[j];
  }
  blen = sqrtf(blen);
  if (blen > this->MaxGradientLength)
  {
    blen = this->MaxGradientLength / blen;
    for (int j = 0; j < DIMENSION; ++j)
    {
      b[j] = b[j] * blen;
    }
  }



	// store the gradient
  for (int j = 0; j < DIMENSION; ++j)
  {
    this->Samples[SampleIndex].Gradient.Elements[j] = b[j];
  }
}




#ifndef CELLSPLIT
void CKdTree::ComputeSampleGradient(const int SampleIndex)
{
	//#ifdef ISOTROPIC
	//	return;
	//#endif

  CVector QueryPosition = this->Samples[SampleIndex].Position;
  float QueryLuminance = this->Samples[SampleIndex].XYZ[1];
  CKNNs& KNNs = this->Samples[SampleIndex].KNNs;



  // compute differences
  std::vector<CVector> PositionDifferences(this->KNNsNum);
  std::vector<float> LuminaceDifferences(this->KNNsNum);

  for (int i = 0; i < this->KNNsNum; ++i)
  {
		int KNNIndex = KNNs.SamplesIndices[i];
    PositionDifferences[i] = this->Samples[KNNIndex].Position - QueryPosition;
    LuminaceDifferences[i] = this->Samples[KNNIndex].XYZ[1] - QueryLuminance;
  }
 


  // construct a matrix and a vector for the least squares method
  float A[DIMENSION * DIMENSION];
  float b[DIMENSION];

  for (int j = 0; j < DIMENSION; ++j)
  {
		for (int i = 0; i < DIMENSION; ++i)
		{
      A[i + j * DIMENSION] = 0.0f;
      for (int k = 0; k < this->KNNsNum; ++k)
      {
        A[i + j * DIMENSION] = A[i + j * DIMENSION] + PositionDifferences[k].Elements[j] * PositionDifferences[k].Elements[i];
      }
    }
  }

  for (int j = 0; j < DIMENSION; ++j)
  {
    b[j] = 0.0f;
    for (int k = 0; k < this->KNNsNum; ++k)
    {
      b[j] = b[j] + PositionDifferences[k].Elements[j] * LuminaceDifferences[k];
    }
  }



  // solve by LAPACK
  int N = DIMENSION;
  int nrhs = 1;
  int lda = DIMENSION;
  int ipiv[DIMENSION];
  int ldb = DIMENSION;
  int info;
  sgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);



	// limit the length of the gradient
  float blen = 0.0f;
  for (int j = 0; j < DIMENSION; ++j)
  {
    blen = blen + b[j] * b[j];
  }
  blen = sqrtf(blen);
  if (blen > this->MaxGradientLength)
  {
    blen = this->MaxGradientLength / blen;
    for (int j = 0; j < DIMENSION; ++j)
    {
      b[j] = b[j] * blen;
    }
  }



	// store the gradient
  for (int j = 0; j < DIMENSION; ++j)
  {
    this->Samples[SampleIndex].Gradient.Elements[j] = b[j];
  }
}
#endif
