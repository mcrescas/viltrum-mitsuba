#ifndef __KDTREE_H
#define __KDTREE_H

#include "global.h"
#include "vector.h"
#include "bsphere.h"
#include "mt.h"
#include "heap.h"

#include <vector>
#include <algorithm>
#include "vector.h"
#include "bsphere.h"
#include "bbox.h"

// defined: luminance based metric
// undefined: gradient based metric (not used anymore)
#define ISOTROPIC

// defined: cells based strategy
// undefined: kNNs based strategy
#define CELLSPLIT


class CKNNs
{
  public:
    std::vector<int> SamplesIndices; 
    std::vector<float> SamplesDistances;
    float MaxSquaredDistance;
    int MaxSampleIndex;
    CHeap Heap;
    bool Updated;

  private:

};


using Float_ = float;

template <int DIMENSION>
class CSample
{
  public:
    CVector<DIMENSION> Position;
    Float_ XYZ[3];
    CVector<DIMENSION> Gradient;
		CKNNs KNNs;

#ifndef CELLSPLIT
    CVector Gradient;
    float Distance;
		CBSphere BSphere;
		CKNNs KNNs;
#endif


  private:

};





enum ENodeType 
{ 
  NodeTypeSplit, 
  NodeTypeLeaf
};





template <int DIMENSION>
class CNode
{
  public:
    ENodeType NodeType;
    float Threshold;
    int LeftChild, RightChild;
    std::vector<int> SamplesIndices;
    int Axis;
#ifndef CELLSPLIT
		CBSphere BSphere;
		int Parent;
#endif
		CBBox<DIMENSION> BBox;
    float Error;

  private:

};








template <int DIMENSION>
class CKdTree
{
  public:
    std::vector<CNode<DIMENSION>> Nodes;
    std::vector<CSample<DIMENSION>> Samples;
    std::vector<CNode<DIMENSION>> Nodes2D;
    CHeap Heap;

    int BucketSize;
    int KNNsNum;
		int CandidatesNum;
		float MaxGradientLength;
		CVector<DIMENSION> SampleExtent;
    int CurrentNumNodes;
    int CurrentNumSamples;
    int CurrentNum2DNodes;
    float DistanceEpsilon;
    float LuminanceLimit;
		CVector<DIMENSION> OffSet;


    CKdTree() 
    {
      this->BucketSize = 4;
      this->KNNsNum = 4;//DIMENSION + 1;
			this->CandidatesNum = 4;//this->KNNsNum * this->KNNsNum;
			this->MaxGradientLength = 5.0f;
      this->CurrentNumNodes = 0;
      this->CurrentNumSamples = 0;
      this->DistanceEpsilon = 1E-5f;
      this->LuminanceLimit = 100.0f;

			for (int i = 0; i < DIMENSION; i++)
			{
				this->SampleExtent.Elements[i] = 1.0f;
			}
    };


    void Build(const std::vector<CSample<DIMENSION>> SamplesList);
    void Build2D();
    void Insert(CSample<DIMENSION>& Sample, bool Is2D);
    void KNNSearch(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, float SearchRadius = 1E+20f);
    void KNNSearchEllipse(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv, float SearchRadius = 1E+20f);
    void KNNSearchAAEllipse(const int RequiredKNNsNum, const CVector<DIMENSION>& QueryPosition, CKNNs* KNNs, const CVector<DIMENSION> & scale, float SearchRadius = 1E+20f);
    int NearestSearch(const CVector<DIMENSION>& QueryPosition, float SearchRadius = 1E+20f);
    void ComputeSampleGradientKNNs(const int SampleIndex, const CKNNs &SKNNs);

#ifndef CELLSPLIT
    void ReverseKNNSearch(const CVector& QueryPosition, CKNNs* KNNs);
    void ReverseSearchIteration(const int NodeIndex, const CVector & QueryPosition, CKNNs* KNNs);
		void ComputeBSpheres(); 
		CBSphere ComputeBSpheresIteration(const int NodeIndex);
    void ComputeSampleGradient(const int SampleIndex);
		inline void ComputeSampleDistance(const int SampleIndex);

    ///
		inline float ComputeDistance(const CVector& Position, const int SampleIndex, const CKNNs& KNNs, const int ParentIndex);
		void ComputeKNNs();
    void ComputeSampleKNNs(const int SampleIndex);
#endif

    CVector<DIMENSION> GenerateNewSamplePosition(MTRand& mtrand, bool Is2D);


    // cells based adaptive sampling
    void ComputeErrors();
    void ComputeNodeError(const int NodeIndex, bool Is2D);
    int MakeNode(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox);
    int MakeNode2D(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox);
    int MakeNodeRecon(std::vector<int> SamplesIndices, const int Axis, const CBBox<DIMENSION> BBox);

    void SortSamples(std::vector<int>& SamplesIndices, int Axis, int li, int ri);

private:
    void SearchTest(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs);
    void SearchTestEllipse(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv);
    void SearchTestAAEllipse(const int SampleIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const CVector<DIMENSION> & scale);
    void SearchIteration(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, float DistanceToCell = 0.0f);
    void SearchIterationEllipse(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, const CMatrix<DIMENSION> & M, const CMatrix<DIMENSION> & MInv);
    void SearchIterationAAEllipse(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, CKNNs* KNNs, const int RequiredKNNsNum, const CVector<DIMENSION> & scale);
    void NearestSearchIteration(const int NodeIndex, const CVector<DIMENSION> & QueryPosition, int& NearestIndex, float& MinSquareDistance);
};

#include "skdtree.cpp"



//inline float CKdTree::ComputeDistance(const CVector& Position, const int SampleIndex, const CKNNs& KNNs, const int ParentIndex)
//{
//	// compute distance between Position and Samples[SampleIndex] by using KNNs for the metric computation
//
//	#ifdef ISOTROPIC
///*   
//		// the compatible isotropic metric
//		// compute Euclidean distance between Position and Samples[SampleIndex]
//		CVector Difference = this->Samples[SampleIndex].Position - Position;
//		float Distance = Difference.DotProduct(Difference);
//
//		// compute gradient-distorted distance without cosine term of dot product
//		float AverageScaledSquaredDistance = 0.0f;
//		for (int j = 0; j < this->KNNsNum; j++)
//		{
//			const CVector& NeighboringGradient = this->Samples[KNNs.SamplesIndices[j]].Gradient;
//			float ScaledSquaredDistance = Distance * NeighboringGradient.DotProduct(NeighboringGradient);
//
//			AverageScaledSquaredDistance = AverageScaledSquaredDistance + ScaledSquaredDistance;
//		}
//		AverageScaledSquaredDistance = AverageScaledSquaredDistance / (float)(this->KNNsNum);
//		Distance = 1E-8f * Distance + AverageScaledSquaredDistance;
//*/
//		
//
//    // using luminance difference (error based metric)
//		// compute Euclidean distance
//		CVector Difference = this->Samples[SampleIndex].Position - Position;
//		float Distance = Difference.DotProduct(Difference);
///*
//		// compute average luminance
//		float AverageLuminance = 0.0f;
//		for (int j = 0; j < this->KNNsNum; ++j)
//		{
//			AverageLuminance = AverageLuminance + this->Samples[KNNs.SamplesIndices[j]].Luminance;
//		}
//		AverageLuminance = AverageLuminance / (float)(this->KNNsNum);
//
//		float MinLuminance = 1.0E+20f;
//		float MaxLuminance = 0.0f;
//*/
//  float AverageLuminance[3];
//  AverageLuminance[0] = 0.0f;
//  AverageLuminance[1] = 0.0f;
//  AverageLuminance[2] = 0.0f;
//		for (int j = 0; j < this->KNNsNum; ++j)
//  {
//    for (int c = 0; c < 3; c++)
//    {
//    AverageLuminance[c] += fabs(this->Samples[KNNs.SamplesIndices[j]].XYZ[c]);
//    }
//  }
//  AverageLuminance[0] /= (float)(this->KNNsNum);
//  AverageLuminance[1] /= (float)(this->KNNsNum);
//  AverageLuminance[2] /= (float)(this->KNNsNum);
//
//		// compute the average of squared differences of luminance
//		float AverageSquaredDifference = 0.0f;
//		for (int j = 0; j < this->KNNsNum; ++j)
//		{
//      for (int c = 0; c < 3; c++)
//      {
//			float Luminance = this->Samples[KNNs.SamplesIndices[j]].XYZ[c];
//			//if (Luminance < MinLuminance) MinLuminance = Luminance;
//			//if (Luminance > MaxLuminance) MaxLuminance = Luminance;
//
//			//float SquaredDifference = (Luminance - AverageLuminance);
//			float SquaredDifference = fabs(Luminance - this->Samples[ParentIndex].XYZ[c]);
//			AverageSquaredDifference = AverageSquaredDifference + SquaredDifference / (AverageLuminance[c] + 1e-10f);// * SquaredDifference;
//      }
//		}
//		AverageSquaredDifference = AverageSquaredDifference / (float)(this->KNNsNum);
//		//if (AverageSquaredDifference > this->LuminanceLimit) AverageSquaredDifference = this->LuminanceLimit; 
//
//		//AverageSquaredDifference = AverageSquaredDifference / ((MaxLuminance - MinLuminance) * (MaxLuminance - MinLuminance) + 1E-20f);
//		Distance = this->DistanceEpsilon * Distance + Distance * AverageSquaredDifference;
//
//	#else
//		// the anisotropic metric
//		// compute Euclidean distance
//		CVector Difference = this->Samples[SampleIndex].Position - Position;
//		float Distance = Difference.DotProduct(Difference);
//
//		// compute gradient-distorted distance
//		float AverageDistortedSquaredDistance = 0.0f;
//		for (int j = 0; j < this->KNNsNum; j++)
//		{
//			const CVector& NeighboringGradient = this->Samples[KNNs.SamplesIndices[j]].Gradient;
//			float DistortedSquaredDistance = Difference.DotProduct(NeighboringGradient);
//			DistortedSquaredDistance = DistortedSquaredDistance * DistortedSquaredDistance;
//
//			AverageDistortedSquaredDistance = AverageDistortedSquaredDistance + DistortedSquaredDistance;
//		}
//		AverageDistortedSquaredDistance = AverageDistortedSquaredDistance / (float)(this->KNNsNum);
//		Distance = 1E-8f * Distance + AverageDistortedSquaredDistance;
//	#endif
//
//  return Distance;
//}
//
//
//
//inline void CKdTree::ComputeSampleDistance(const int SampleIndex)
//{
//	this->Samples[SampleIndex].Distance = 0.0f;
//	for (int i = 0; i < this->KNNsNum; ++i)
//	{
//		float Distance = this->ComputeDistance(this->Samples[SampleIndex].Position, this->Samples[SampleIndex].KNNs.SamplesIndices[i], this->Samples[SampleIndex].KNNs, SampleIndex);
//		this->Samples[SampleIndex].Distance = this->Samples[SampleIndex].Distance + Distance;
//	}
//}
// 


#endif
