// Samoa2
#include <ctime>
#include <iostream>
#include "dynload.h"

#include "sampling.h"
#include "paramset.h"
#include "film.h"

#include "mt.h"
#include "vector.h"
#include "kdtree.h"
#include <vector>
#include <fstream>
#include <iostream>

 extern "C"
 {
   #include "f2c.h"
   #include "clapack.h"
 }

struct Splat
{
  float params[5];
  Spectrum L0, L1, L2, L3;
};


struct PixelCount
{
  std::vector<int>Px;
  std::vector<int>Py;
};



// #define AREA_LIGHT



#define VORONOI


 


class Samoa2Sampler: public Sampler 
{
  public:
	  Samoa2Sampler(int candidatenum, int knnsnum, int bucketnum, int cellknnsnum, float maxgs, int optstep, int initnum, float imagescale, float luminancelimit, float luminancescale, float distanceepsilon, int xstart, int xend, int ystart, int yend, int ps, int width, int height, string reconfile, string densityfile);
    ~Samoa2Sampler(){}

	  int RoundSize(int size) const;

	  bool GetNextSample(Sample *sample);
	  Spectrum SetPrevSample(Spectrum Ls);

    CKdTree KdTree;
    MTRand Rand;
    std::vector<CSample> InitialSamples;
    std::vector<Spectrum> SamplesLs;
    int InitialSampleNum;
    CSample CurrentSample;

    int NumAdaptiveSamples;
    float LuminanceScale;
    int CellKNNsNum;
    float MaxGaussianScale;
    int NumOptimizationStep;
    float Anisotropy; 


  private:
    int width;
    int height;

	  int totalSamples;
	  int initialSamples;
	  int sampleNum;
	  FILE * sampfile;
  	
	  clock_t startTime,finishTime;
	  float elapsedTimeGet, elapsedTimeSet;
    
    std::string reconfile;
    std::string densityfile;
};





int Samoa2Sampler::RoundSize(int size) const 
{
  return size;
}





Samoa2Sampler::Samoa2Sampler(int candidatenum, int knnsnum, int bucketnum, int cellknnsnum, float maxgs, int optstep, int initnum, float imagescale, float luminancelimit, float luminancescale, float distanceepsilon, int xstart, int xend, int ystart, int yend, int ps, int width, int height, string reconfile, string densityfile)
: Sampler(xstart, xend, ystart, yend, ps), reconfile(reconfile), densityfile(densityfile)
{
  this->sampleNum = 0;
	this->totalSamples = this->TotalSamples();

  this->KdTree.CandidatesNum = candidatenum;
  this->KdTree.KNNsNum = knnsnum;
  this->KdTree.BucketSize = bucketnum;
  this->InitialSampleNum = initnum;
  this->initialSamples = initnum;

  this->width = (this->xPixelEnd - this->xPixelStart - 4);
  this->height = (this->yPixelEnd - this->yPixelStart - 4);

  float ImageWidth = (float)width;
  float ImageHeight = (float)height;


  float MinImageAxis = (float)(this->width);
  if (this->height < MinImageAxis)
  {
    MinImageAxis = (float)(this->height);
  }
  this->KdTree.SampleExtent.Elements[0] = ((float)this->width / MinImageAxis) * imagescale * (float)((float)this->width / ImageWidth);
  this->KdTree.SampleExtent.Elements[1] = ((float)this->height / MinImageAxis) * imagescale * (float)((float)this->height / ImageHeight);
  this->LuminanceScale = luminancescale;
  this->CellKNNsNum = cellknnsnum;
  this->KdTree.DistanceEpsilon = distanceepsilon;
  this->KdTree.LuminanceLimit = luminancelimit;

  // not used anymore
  this->MaxGaussianScale = maxgs;
  this->NumOptimizationStep = optstep;

  std::cout << this->xPixelStart << std::endl;
  std::cout << this->xPixelEnd << std::endl;
  std::cout << this->yPixelStart << std::endl;
  std::cout << this->yPixelEnd << std::endl;

  // 
  this->Anisotropy = maxgs;

  std::cout << "candidate: " << this->KdTree.CandidatesNum << std::endl;
  std::cout << "knns: " << this->KdTree.KNNsNum << std::endl;
  std::cout << "bucket: " << this->KdTree.BucketSize << std::endl;
  std::cout << "initial: " << this->InitialSampleNum << std::endl;
  std::cout << "scalex: " << this->KdTree.SampleExtent.Elements[0] << std::endl;
  std::cout << "scaley: " << this->KdTree.SampleExtent.Elements[1] << std::endl;
  std::cout << "scalelum: " << this->LuminanceScale << std::endl;
  std::cout << "cellknns: " << this->CellKNNsNum << std::endl;
  std::cout << "total: " << this->TotalSamples()  << std::endl;
  std::cout << "w: " << width  << std::endl;
  std::cout << "h: " << height << std::endl;

  Rand.seed(570905);

  // allocate all memory
  this->KdTree.Nodes.resize(this->TotalSamples() * 3);
  this->KdTree.Samples.resize(this->TotalSamples());
  this->KdTree.Heap.Items.resize(this->TotalSamples() * 2 + 1);
  this->KdTree.Heap.Indices.resize(this->TotalSamples() * 2 + 1);
  this->SamplesLs.resize(this->TotalSamples());


  //this->KdTree.Nodes2D.resize(this->TotalSamples() * 2);
  //this->KdTree.CurrentNum2DNodes = 0;


	//this->sampfile = fopen(samplefile.c_str(),"w");


	//this->NumAdaptiveSamples = width * height * 100; //not used anymore
}





bool Samoa2Sampler::GetNextSample(Sample *sample) 
{
	startTime = clock();
	this->sampleNum++;


	if (this->sampleNum == (this->totalSamples)) 
  {
    long recon_start_time = clock();

    std::vector<CHeapItem>().swap(this->KdTree.Heap.Items);
    std::vector<int>().swap(this->KdTree.Heap.Indices);

    //// -------- output sample positions in 3D ---------------
    //std::vector<CVector> points; 
    //for (int i = 0; i < (int)this->KdTree.Samples.size(); i++)
    //{
    //  if ((i % 10) == 0)
    //  {
    //    points.push_back(this->KdTree.Samples[i].Position);
    //  }
    //}

    //std::ofstream pointset("test.mqo", std::ios::out);
    //pointset << "Metasequoia Document\n";
    //pointset << "Format Text Ver 1.0\n";
    //pointset << "\n";
    //pointset << "Scene {\n";
    //pointset << "	pos 0.0000 0.0000 1500.0000\n";
    //pointset << "	lookat 0.0000 0.0000 0.0000\n";
    //pointset << "	head -0.5236\n";
    //pointset << "	pich 0.5236\n";
    //pointset << "	ortho 0\n";
    //pointset << "	zoom2 5.0000\n";
    //pointset << "	amb 0.250 0.250 0.250\n";
    //pointset << "}\n";
    //pointset << "Object ""obj1"" {\n";
    //pointset << "	visible 15\n";
    //pointset << "	locking 0\n";
    //pointset << "	shading 1\n";
    //pointset << "	facet 59.5\n";
    //pointset << "	color 0.898 0.498 0.698\n";
    //pointset << "	color_type 0\n";
    //pointset << "	vertex " << points.size() << " {\n";
    //for (int i = 0; i < (int)points.size(); i++)
    //{
    //  pointset << "		";
    //  for (int j = 0; j < DIMENSION; j++)
    //  {
    //    pointset << " " << points[i].Elements[j] * 400.0f - 200.0f;
    //  }
    //  pointset << "\n";
    //}
    ////pointset << "		-26.8427 124.7665 -67.6795\n";
    ////pointset << "		-19.3180 115.2255 -65.6633\n";
    ////pointset << "		-111.4174 112.6568 -10.7772\n";
    ////pointset << "		-107.1877 100.5471 -5.1461\n";
    //pointset << "	}\n";
    //pointset << "	face " << points.size() << " {\n";
    //for (int i = 0; i < (int)points.size(); i++)
    //{
    //  pointset << "		2 V(" << i << " " << i << ")\n";
    //}
    ////pointset << "		2 V(0 1)\n";
    ////pointset << "		2 V(2 3)\n";
    //pointset << "	}\n";
    //pointset << "}\n";
    //pointset << "Eof\n";




	  char *searchpath = getenv("PBRT_SEARCHPATH");
	  if (!searchpath)
	  {
      std::cerr << "You need to set PBRT_SEARCHPATH." << std::endl;
		  return false;
	  }
	  UpdatePluginPath(searchpath);
  	
	  int w = (int)(width);
	  int h = (int)(height);	
	  ParamSet film_params;

	  film_params.AddInt("xresolution", &w);
	  film_params.AddInt("yresolution", &h);
	  film_params.AddString("filename", &reconfile);


    // -------------------- reconstruction --------------------
	  Film* film = MakeFilm( "samoa_image", film_params,0);
	  float params[5];
	  //for (int i = 0; i < (int)this->KdTree.Samples.size(); i++) 
    //{
      /*
      // compute the maximum gradient length
      float mg = 0.0f;
      for (int e = 0; e < this->KdTree.KNNsNum; e++)
      {
        int KNNIndex = this->KdTree.Samples[i].KNNs.SamplesIndices[e];
        float g = this->KdTree.Samples[KNNIndex].Gradient.DotProduct(this->KdTree.Samples[KNNIndex].Gradient);
        if (g > mg) 
        {
          mg = g;
        }
      }
      mg = mg + 1E-20f;

      double A[DIMENSION * DIMENSION];
      int m=DIMENSION,n=DIMENSION,lda=DIMENSION,info,piv[DIMENSION],lwork=DIMENSION;
      double work[DIMENSION];

      for (int j = 0; j < DIMENSION; ++j)
      {
        for (int k = 0; k < DIMENSION; ++k)
        {
          int index = j + k * DIMENSION;
          A[index] = 0.0;
          for (int e = 0; e < this->KdTree.KNNsNum; e++)
          {
            int KNNIndex = this->KdTree.Samples[i].KNNs.SamplesIndices[e];
            A[index] = A[index] + this->KdTree.Samples[KNNIndex].Gradient.Elements[j] * this->KdTree.Samples[KNNIndex].Gradient.Elements[k] / mg;
          }

          if (j == k)
          {
            A[index] = A[index] + 1.0f;
          }
        }
      }
      dgetrf_( &m, &n, A, &lda, piv, &info);
      dgetri_(&n, A, &lda, piv, work, &lwork, &info);
      */

/*
  		params[0] = this->KdTree.Samples[i].BSphere.SquaredRadius;
  		params[1] = 0.1f;
      params[2] = 0.1f;
      params[3] = 0.1f;
      params[4] = 10.0f;
*/
/*
	float s2_c = 0.3f;
	float q = 0.75f;
	float m = sqrtf(this->KdTree.Samples[i].BSphere.SquaredRadius) * sqrtf(this->width * this->height);

	float gx2, gy2, gb2;
	gx2 = gy2 = gb2 = 0.f;

	float stretch2_inv = 1/(1+gx2+gy2+gb2);
	float s2_x = (m * q)*(m * q);
	float s2_y = s2_x * stretch2_inv * (1+gb2);
	float scx = s2_x+s2_c;
	float scy = s2_y+s2_c;
	float k = m*m* sqrtf(scx*scy*stretch2_inv) / (s2_x*s2_y);
	k*=m*m;

	float scx_inv = 1/scx;
	float scy_inv = 1/scy;
	
	float gpl2 = gx2+gy2;
	
	float ghx2 = 0.;
	float ghy2 = 1.;
	float ghxy = 0.;
	
	//if (gpl2 > 0) {
	//	float gpl2_inv = 1/gpl2;
	//	ghx2 = gx2*gpl2_inv;
	//	ghy2 = gy2*gpl2_inv;
	//	ghxy = gradient[0]*gradient[1]*gpl2_inv;
	//}
	
	float xxF = (ghy2*scx_inv + ghx2*scy_inv);
	float xyF = -2*ghxy*(scx_inv-scy_inv);
	float yyF = (ghx2*scx_inv + ghy2*scy_inv);

	
	params[0]=k;
	params[1]=xxF;
	params[2]=xyF;
	params[3]=yyF;
	params[4]=m*q;
*/

/*
		  Sample sample(this->KdTree.Samples[i].Position.Elements[0] * (float)this->width,
                    this->KdTree.Samples[i].Position.Elements[1] * (float)this->height);
		  Ray unused;

		  film->AddSample( sample, unused, SamplesLs[i], 1, params);
	  }
*/





   // Splat splt;
	  //for (int i = 0; i < (int)this->KdTree.Nodes.size(); ++i)
	  //{
		 // if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) continue;

   //   float x0 = (this->KdTree.Nodes[i].BBox.Min.Elements[0] * (float)this->width);
   //   float y0 = (this->KdTree.Nodes[i].BBox.Min.Elements[1] * (float)this->height);
   //   float x1 = (this->KdTree.Nodes[i].BBox.Max.Elements[0] * (float)this->width);
   //   float y1 = (this->KdTree.Nodes[i].BBox.Max.Elements[1] * (float)this->height);

		 // float v = 1.0f;
		 // for (int j = 2; j < DIMENSION; j++)
		 // {
			//  v = v * (this->KdTree.Nodes[i].BBox.Max.Elements[j] - this->KdTree.Nodes[i].BBox.Min.Elements[j]);
		 // }

	  //  splt.params[0]=x0;
	  //  splt.params[1]=y0;
	  //  splt.params[2]=x1;
	  //  splt.params[3]=y1;
	  //  splt.params[4]=v;

   //   CVector QueryPosition;
   //   CKNNs KNNs;
      


		  //int smp = this->KdTree.Nodes[i].SamplesIndices[0];
	   // Sample sample(this->KdTree.Samples[smp].Position.Elements[0] * (float)this->width,
    //                this->KdTree.Samples[smp].Position.Elements[1] * (float)this->height);
	   // Ray unused;

    //  Spectrum La;
    //  La = SamplesLs[this->KdTree.Nodes[i].SamplesIndices[0]];
    //  for (int j = 1; j < (int)this->KdTree.Nodes[i].SamplesIndices.size(); j++)
    //  {
    //    La = La + SamplesLs[this->KdTree.Nodes[i].SamplesIndices[j]];
    //  }
    //  La = La / (float)(this->KdTree.Nodes[i].SamplesIndices.size());

	   // film->AddSample( sample, unused, La, 1,&splt);
    //}

/*
    // splatting
    for (int y = (int)y0; y < (int)y1; y++)
    {
      if ((y >= this->height) || (y < 0)) continue;
      for (int x = (int)x0; x < (int)x1; x++)
      {
        if ((x >= this->width) || (x < 0)) continue;

				QueryPosition.Elements[0] = x / (float)(this->width);
				QueryPosition.Elements[1] = y / (float)(this->height);
				QueryPosition.Elements[2] = 0.5f * (this->KdTree.Nodes[i].BBox.Min.Elements[2] + this->KdTree.Nodes[i].BBox.Max.Elements[2]);
				this->KdTree.KNNSearch(1, QueryPosition, &KNNs);
        Spectrum La = SamplesLs[KNNs.SamplesIndices[0]];

	      Sample sample(x, y);
	      Ray unused;

        film->AddSample( sample, unused, La, 1,&splt);
      }
    }
*/
    Splat splt;



    // write out sample distribution
	  ParamSet film_params2;

	  film_params2.AddInt("xresolution", &width);
	  film_params2.AddInt("yresolution", &height);
	  film_params2.AddString("filename", &densityfile);
	  Film* film2 = MakeFilm("samoa_image_", film_params2, 0);
    splt.params[0]=0;
    splt.params[1]=0;
    splt.params[2]=0;
    splt.params[3]=0;
    splt.params[4]=1.0f;
    Spectrum WhiteL;
    WhiteL = 1.0f;
    for (int i = 0; i < this->KdTree.CurrentNumSamples; ++i)
    {
		  int x = (int)(this->KdTree.Samples[i].Position.Elements[0] * this->width / this->KdTree.SampleExtent.Elements[0]);
		  int y = (int)(this->KdTree.Samples[i].Position.Elements[1] * this->height / this->KdTree.SampleExtent.Elements[1]);

	    Ray unused;
	    Sample sample(x, y);
      //film2->AddSample( sample, unused, SamplesLs[i], 1,&splt);
      film2->AddSample( sample, unused, WhiteL, 1,&splt);
    }
    film2->WriteImage();




    //this->KdTree.Build2D();

    printf("\n");
    this->KdTree.BucketSize = 1;
    int CurrenNumtKdtreeNodes = (int)this->KdTree.CurrentNumNodes;
	  for (int i = 0; i < CurrenNumtKdtreeNodes; ++i)
	  {
		  if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) continue;
		  if (this->KdTree.Nodes[i].SamplesIndices.size() == 0) continue;

      //if (this->KdTree.Nodes[i].BBox.GetLargestAxis() < 2) continue;

      this->KdTree.MakeNode(this->KdTree.Nodes[i].SamplesIndices, this->KdTree.Nodes[i].Axis, this->KdTree.Nodes[i].BBox);
      if ((i % 1000) == 0) printf("Subdividing Nodes: [%d%%]\r", (int)((i / (float)CurrenNumtKdtreeNodes) * 100.0f));
    }
    printf("Subdividing Nodes: [%d%%]\r", 100);


    printf("\n");
    this->KdTree.KNNsNum = this->CellKNNsNum;
    for (int i = 0; i < this->KdTree.CurrentNumSamples; i++)
    {
      this->KdTree.KNNSearch(this->KdTree.KNNsNum + 1, this->KdTree.Samples[i].Position, &(this->KdTree.Samples[i].KNNs));
      this->KdTree.ComputeSampleGradientKNNs(i, this->KdTree.Samples[i].KNNs);
      if ((i % 1000) == 0) printf("Computing kNNs and gradients: [%d%%]\r", (int)((i / (float)this->KdTree.CurrentNumSamples) * 100.0f));
    }
    printf("Computing kNNs and gradients: [%d%%]\r", 100);




   // std::vector<float> ProjectedVariance(this->width * this->height);
	  //std::vector<float> ProjectedAverage(this->width * this->height);
	  //std::vector<float> ProjectedCount(this->width * this->height);

	  //for (int ix = 0; ix < (this->width * this->height); ix++)
	  //{
		 // ProjectedVariance[ix] = 0.0f;
		 // ProjectedAverage[ix] = 0.0f;
		 // ProjectedCount[ix] = 0.0f;
	  //}

	  //for (int i = 0; i < this->KdTree.CurrentNumSamples; i++)
	  //{
		 // int ix = (int)(this->KdTree.Samples[i].Position.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width);
		 // int iy = (int)(this->KdTree.Samples[i].Position.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height);

   //   if ((ix < 0) || (ix >= this->width) || (iy < 0) || (iy >= this->height))
		 // {
			//  continue;
		 // }
		 // ProjectedAverage[ix + iy * this->width] += this->KdTree.Samples[i].XYZ[1];
		 // ProjectedCount[ix + iy * this->width] += 1.0f;
	  //}
	  //for (int ix = 0; ix < (this->width * this->height); ix++)
	  //{
		 // if (ProjectedCount[ix] > 0.0f)
		 // {
			//  ProjectedAverage[ix] /= ProjectedCount[ix];
		 // }
	  //}

	  //for (int i = 0; i < this->KdTree.CurrentNumSamples; i++)
	  //{
		 // int ix = (int)(this->KdTree.Samples[i].Position.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width);
		 // int iy = (int)(this->KdTree.Samples[i].Position.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height);

   //   if ((ix < 0) || (ix >= this->width) || (iy < 0) || (iy >= this->height))
		 // {
			//  continue;
		 // }
		 // ProjectedVariance[ix + iy * this->width] += fabs(this->KdTree.Samples[i].XYZ[1] - ProjectedAverage[ix + iy * this->width]);
	  //}
	  //for (int ix = 0; ix < (this->width * this->height); ix++)
	  //{
		 // if (ProjectedCount[ix] > 0.0f)
		 // {
			//  ProjectedVariance[ix] /= ProjectedCount[ix];
   //     //ProjectedVariance[ix] = sqrtf(ProjectedVariance[ix]);
		 // }
	  //}

    // write out sample distribution
	  string filename3 = "aniso.exr";
	  ParamSet film_params3;

	  film_params3.AddInt("xresolution", &width);
	  film_params3.AddInt("yresolution", &height);
	  film_params3.AddString("filename", &filename3);
	  Film* film3 = MakeFilm("samoa_image", film_params3, 0);
    splt.params[0]=0;
    splt.params[1]=0;
    splt.params[2]=0;
    splt.params[3]=0;
    splt.params[4]=1.0f;
    for (int i = 0; i < this->KdTree.CurrentNumSamples; ++i)
    {
		  int x = (int)(this->KdTree.Samples[i].Position.Elements[0] * width / this->KdTree.SampleExtent.Elements[0]);
		  int y = (int)(this->KdTree.Samples[i].Position.Elements[1] * height / this->KdTree.SampleExtent.Elements[1]);
      
      float val = 0.0f;
      float MaxScaling = 0.0f;
      float MinScaling = 1e+10f;
      CVector ScalingFactors;
      for (int mm = 0; mm < this->CellKNNsNum; mm++)
	    {
        for (int j = 0; j < DIMENSION; j++)
        {
		      float cval = this->KdTree.Samples[this->KdTree.Samples[i].KNNs.SamplesIndices[mm]].Gradient.Elements[j];
          ScalingFactors.Elements[j] += cval * cval;
          val += cval * cval;
        }
	    }

      for (int j = 0; j < DIMENSION; j++)
      {
        if (ScalingFactors.Elements[j] > MaxScaling)
        {
          MaxScaling = ScalingFactors.Elements[j];
        }
        if (ScalingFactors.Elements[j] < MinScaling)
        {
          MinScaling = ScalingFactors.Elements[j];
        }
      }


      //ScalingFactors = ScalingFactors * (1.0f / (val + 1e-10f));
      //ScalingFactors.Elements[0] = ScalingFactors.Elements[2];
      //ScalingFactors.Elements[1] = ScalingFactors.Elements[2];

      Spectrum Ls(MaxScaling / (MinScaling + 1e-20f));//ScalingFactors.Elements);

	    Ray unused;
	    Sample sample(x, y);
      film3->AddSample(sample, unused, Ls, 1, &splt);
    }
    film3->WriteImage();






   // // write out sample distribution
	  //string filename4 = "spread.exr";
	  //ParamSet film_params4;

	  //film_params4.AddInt("xresolution", &width);
	  //film_params4.AddInt("yresolution", &height);
	  //film_params4.AddString("filename", &filename4);
	  //Film* film4 = MakeFilm("samoa_image", film_params4, 0);
   // splt.params[0]=0;
   // splt.params[1]=0;
   // splt.params[2]=0;
   // splt.params[3]=0;
   // splt.params[4]=1.0f;
   // for (int i = 0; i < this->KdTree.CurrentNumSamples; ++i)
   // {
		 // int x = (int)(this->KdTree.Samples[i].Position.Elements[0] * this->width / this->KdTree.SampleExtent.Elements[0]);
		 // int y = (int)(this->KdTree.Samples[i].Position.Elements[1] * this->height / this->KdTree.SampleExtent.Elements[1]);

   //   float MaxDistance2D = 0.0f;
   //   for (int j = 0; j < this->KdTree.KNNsNum; j++)
   //   {
   //     float dx = this->KdTree.Samples[i].Position.Elements[0] - this->KdTree.Samples[this->KdTree.Samples[i].KNNs.SamplesIndices[j]].Position.Elements[0];
   //     float dy = this->KdTree.Samples[i].Position.Elements[1] - this->KdTree.Samples[this->KdTree.Samples[i].KNNs.SamplesIndices[j]].Position.Elements[1];
   // 
   //     float Distance2D = dx * dx + dy * dy;
   //     //if (MaxDistance2D < Distance2D)
   //     //{
   //     //  MaxDistance2D = Distance2D;
   //     //}
   //     MaxDistance2D += sqrtf(Distance2D) * 10.0f;
   //   }
   //   MaxDistance2D /= (float)(this->KdTree.KNNsNum);

	  //  Ray unused;
	  //  Sample sample(x, y);
   //   Spectrum Ls(MaxDistance2D);
   //   film4->AddSample( sample, unused, Ls, 1,&splt);
   // }
   // film4->WriteImage();


    std::vector<PixelCount> SamplePixelCount(this->KdTree.CurrentNumSamples);
    std::vector<float> SamplePixelDistance(this->KdTree.CurrentNumSamples);
    for (int i = 0; i < this->KdTree.CurrentNumSamples; i++)
    {
      SamplePixelCount[i].Px.clear();
      SamplePixelCount[i].Py.clear();
      SamplePixelDistance[i] = 0.0f;
    }



    printf("\n");
    int KdTreeCellKNNsNum = this->CellKNNsNum;
    // --------------------- splatting --------------------
	  for (int i = 0; i < (int)this->KdTree.CurrentNumNodes; ++i)
	  {
		  if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) continue;
		  if (this->KdTree.Nodes[i].SamplesIndices.size() == 0) continue;

      float x0 = (this->KdTree.Nodes[i].BBox.Min.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width);
      float y0 = (this->KdTree.Nodes[i].BBox.Min.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height);
      float x1 = (this->KdTree.Nodes[i].BBox.Max.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width);
      float y1 = (this->KdTree.Nodes[i].BBox.Max.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height);

      int ix0 = (int)x0;
      int ix1 = (int)x1;
      int iy0 = (int)y0;
      int iy1 = (int)y1;



		  float v = 1.0f;
		  for (int j = 2; j < DIMENSION; j++)
		  {
			  v = v * (this->KdTree.Nodes[i].BBox.Max.Elements[j] - this->KdTree.Nodes[i].BBox.Min.Elements[j]);
		  }

	    splt.params[0] = x0;
	    splt.params[1] = y0;
	    splt.params[2] = x1;
	    splt.params[3] = y1;
	    splt.params[4] = v;

      CVector QueryPosition;
      //CKNNs KNNs;

      // splatting
		  float SearchRadius = 1E+20f;
		  int PrevIndex = -1;
		  float RevWidth = this->KdTree.SampleExtent.Elements[0] / (float)(this->width);
		  float RevHeight = this->KdTree.SampleExtent.Elements[1] / (float)(this->height);

      Spectrum La;
	    Ray unused;

#ifndef VORONOI

      this->KdTree.KNNsNum = this->CellKNNsNum;
      this->KdTree.ComputeSampleKNNs(this->KdTree.Nodes[i].SamplesIndices[0]);
      this->KdTree.ComputeSampleGradient(this->KdTree.Nodes[i].SamplesIndices[0]);	

      CVector Gradient = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].Gradient;


      
      Gradient.Elements[0] = 0.0f;
      Gradient.Elements[1] = 0.0f;
      Gradient.Elements[2] = 1.0f;


      //La = SamplesLs[this->KdTree.Nodes[i].SamplesIndices[0]];
      //for (int j = 1; j < (int)this->KdTree.Nodes[i].SamplesIndices.size(); j++)
      //{
      //  La = La + SamplesLs[this->KdTree.Nodes[i].SamplesIndices[j]];
      //}
      //La = La / (float)(this->KdTree.Nodes[i].SamplesIndices.size());
      Gradient.Normalize();
#endif

		  //for (int j = 0; j < snum; ++j)
		  //{
			  QueryPosition = 0.5f * (this->KdTree.Nodes[i].BBox.Min + this->KdTree.Nodes[i].BBox.Max);
 
			  //QueryPosition.Elements[0] = 0.5f * (this->KdTree.Nodes[i].BBox.Min.Elements[0] + this->KdTree.Nodes[i].BBox.Max.Elements[0]);
			  //QueryPosition.Elements[1] = 0.5f * (this->KdTree.Nodes[i].BBox.Min.Elements[1] + this->KdTree.Nodes[i].BBox.Max.Elements[1]);
			  //QueryPosition.Elements[2] = 0.5f * (this->KdTree.Nodes[i].BBox.Min.Elements[2] + this->KdTree.Nodes[i].BBox.Max.Elements[2]);
			  //QueryPosition.Elements[2] = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[j]].Position.Elements[2];

        //float QueryXYZ[3];
        //for (int c = 0; c < this->KdTree.Nodes[i].SamplesIndices.size(); c++)
        //{
        //  QueryXYZ[c] = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].XYZ[c];
        //}

#ifdef VORONOI
        //CKNNs CellKNNs;
        //this->KdTree.KNNSearch(this->CellKNNsNum, QueryPosition, &CellKNNs);



        CKNNs &CellKNNs = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].KNNs;


      //int KNNs2DNum = this->CellKNNsNum;
      //CKNNs& KNNs2D = CellKNNs;
      ////this->KNNSearch(KNNs2DNum, QueryPosition, &KNNs2D);

      //CVector ScalingFactor;
      //for (int ax = 2; ax < DIMENSION; ax++)
      //{
      //  ScalingFactor.Elements[ax] = 0.0f;
      //  this->KdTree.SortSamples(KNNs2D.SamplesIndices, ax, 0, KNNs2DNum - 1);

      //  for (int c = 0; c < 3; c++)
      //  {
      //    for (int mm = 1; mm < KNNs2DNum; mm++)
      //    {
      //      ScalingFactor.Elements[ax] += fabs(this->KdTree.Samples[KNNs2D.SamplesIndices[mm]].XYZ[c] - this->KdTree.Samples[KNNs2D.SamplesIndices[mm - 1]].XYZ[c]);
      //    }
      //  }
      //  ScalingFactor.Elements[ax] /= (float)(KNNs2DNum) * 3.0f;// * 0.1f;// * 5000.0f;
      //  ScalingFactor.Elements[ax] *= 0.05f;
      //}
      //ScalingFactor.Elements[0] = 1.0f;
      //ScalingFactor.Elements[1] = 1.0f;




      //CVector ScalingFactor;
      //for (int ax = 0; ax < DIMENSION; ax++)
      //{
      //  ScalingFactor.Elements[ax] = 1.0f;
      //}

      //this->KdTree.KNNSearchAAEllipse(this->CellKNNsNum, QueryPosition, &CellKNNs, ScalingFactor);


//        std::vector<CVector> DifferenceVectors(this->CellKNNsNum);
//        for (int m = 0; m < this->CellKNNsNum; m++)
//        {
//          int CellKNNsIndex = CellKNNs.SamplesIndices[m];
//          DifferenceVectors[m] = this->KdTree.Samples[CellKNNsIndex].Position - QueryPosition;
//          DifferenceVectors[m].Normalize();
//
//          float ScalingFactor = 1.0f;
//          for (int c = 0; c < 3; c++)
//          {
//            ScalingFactor += fabs(this->KdTree.Samples[CellKNNsIndex].XYZ[c] - QueryXYZ[c]);
//          }
//          ScalingFactor = pow(ScalingFactor, this->Anisotropy);
//DifferenceVectors[m] = DifferenceVectors[m] * ScalingFactor;
//          //DifferenceVectors[m] = DifferenceVectors[m] * (1.0f / (fabs(this->KdTree.Samples[CellKNNsIndex].XYZ[1] - QueryLuminance) + 1e-10));
//          //DifferenceVectors[m] = DifferenceVectors[m] * (fabs(this->KdTree.Samples[CellKNNsIndex].XYZ[1] - QueryLuminance) + 1e-10);
//        }
//        float DifferenceMatrix[DIMENSION][DIMENSION];
//        for (int mi = 0; mi < DIMENSION; mi++)
//        {
//          for (int mj = 0; mj < DIMENSION; mj++)
//          {
//            DifferenceMatrix[mi][mj] = 0.0f;
//            for (int m = 0; m < this->CellKNNsNum; m++)
//            {
//              DifferenceMatrix[mi][mj] += DifferenceVectors[m].Elements[mi] * DifferenceVectors[m].Elements[mj];
//            }
//          }
//        }


#else
        CKNNs& CellKNNs = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].KNNs;
#endif

#ifndef VORONOI
        float AverageDistance = 0.0f;
        for (int k = 0; k < (int)CellKNNs.SamplesIndices.size(); ++k)
        {
          AverageDistance = AverageDistance + QueryPosition.Distance(this->KdTree.Samples[CellKNNs.SamplesIndices[k]].Position);
        }
        AverageDistance = AverageDistance / (float)(CellKNNs.SamplesIndices.size());

        float RotationMatrix[DIMENSION][DIMENSION];
        CVector RotationBasisVectors[DIMENSION];

        // initialize basis vectors
        RotationBasisVectors[0] = Gradient;
        RotationBasisVectors[0].Normalize();
        for (int ii = 1; ii < DIMENSION; ++ii)
        {
          RotationBasisVectors[ii].UniformRandomHyperSpehere(Rand);
          RotationBasisVectors[ii].Normalize();
        }

        // apply Gram-Schmidt orthonormalization
        for (int ii = 0; ii < DIMENSION; ++ii)
        {
          for (int jj = 0; jj < ii; ++jj)
          {
            RotationBasisVectors[ii] = RotationBasisVectors[ii] - RotationBasisVectors[ii].DotProduct(RotationBasisVectors[jj]) * RotationBasisVectors[jj];
          }
          RotationBasisVectors[ii].Normalize();
        }
 

        // compute optimal scaling
        float GaussianScale = 1.0f;
        float MinimumError = 1.0E20f;
        QueryPosition = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].Position;
        for (int si = 1; si <= this->NumOptimizationStep; si++)
        {
          La = 0.0f;
          float w = 0.0f;
          float ss = 0.0f;
          for (int k = 0; k < (int)CellKNNs.SamplesIndices.size(); ++k)
          {
            CVector ODifference = QueryPosition - this->KdTree.Samples[CellKNNs.SamplesIndices[k]].Position;
            CVector Difference;
            for (int ii = 0; ii < DIMENSION; ++ii)
            {
              Difference.Elements[ii] = ODifference.DotProduct(RotationBasisVectors[ii]);
            }

            if (Gradient.DotProduct(Gradient) != 0.0f)
            {
              ss = (float)si / (float)this->NumOptimizationStep * this->MaxGaussianScale;
              ss = 1.0f / ss;
              
              float sm = pow(1.0f / ss, 1.0f / (float)(DIMENSION));
              Difference.Elements[0] *= ss;
              for (int ii = 1; ii < DIMENSION; ++ii)
              {
                Difference.Elements[ii] *= sm;
              }
            }
            else
            {
              Difference = ODifference;
            }

            for (int ii = 0; ii < DIMENSION; ++ii)
            {
              Difference.Elements[ii] /= AverageDistance;
            }

            float cw = exp(-min(0.5f * Difference.DotProduct(Difference), 30.0f));
            w += cw;
            La = La + SamplesLs[CellKNNs.SamplesIndices[k]] * cw;
          }
          La = La / w;
          La = La - SamplesLs[this->KdTree.Nodes[i].SamplesIndices[0]];



          if (MinimumError > fabs(La.y()))
          {
            MinimumError = fabs(La.y());
            GaussianScale =  ss;
          }
        }
			  QueryPosition = 0.5f * (this->KdTree.Nodes[i].BBox.Min + this->KdTree.Nodes[i].BBox.Max);
#endif


       // std::vector<float> BaseDistances(this->CellKNNsNum * this->CellKNNsNum);
       // for (int k = 0; k < this->CellKNNsNum; ++k)
       // { 
       //   int idx = CellKNNs.SamplesIndices[k];
       //   float val = 0.0f;

			    //CVector DifferenceVector;
			    //DifferenceVector = this->KdTree.Samples[idx].Position - QueryPosition;
       //   DifferenceVector.Elements[0] = 0.0f;
       //   DifferenceVector.Elements[1] = 0.0f;

			    //for (int mm = 0; mm < this->CellKNNsNum; mm++)
			    //{
				   // float cval = DifferenceVector.DotProduct(this->KdTree.Samples[CellKNNs.SamplesIndices[mm]].Gradient);
				   // //val += cval * cval;
       //     BaseDistances[k * this->CellKNNsNum + mm] = cval;
			    //}
       // }




			  for (int x = ix0; x < ix1; ++x)
			  {
          QueryPosition.Elements[0] = (x + 0.5f) * RevWidth;
				  //if ((x >= this->width) || (x < 0)) continue;
				  for (int y = iy0; y < iy1; ++y)
				  {
					  //if ((y >= this->height) || (y < 0)) continue;

#ifndef VORONOI
					  QueryPosition.Elements[0] = (x + 0.5f) * RevWidth * this->KdTree.SampleExtent.Elements[0];
					  QueryPosition.Elements[1] = (y + 0.5f) * RevHeight * this->KdTree.SampleExtent.Elements[1];

            La = 0.0f;
            float w = 0.0f;
            for (int k = 0; k < (int)CellKNNs.SamplesIndices.size(); ++k)
            {
              CVector ODifference = QueryPosition - this->KdTree.Samples[CellKNNs.SamplesIndices[k]].Position;
              CVector Difference;
              for (int ii = 0; ii < DIMENSION; ++ii)
              {
                Difference.Elements[ii] = ODifference.DotProduct(RotationBasisVectors[ii]);
              }

              if (Gradient.DotProduct(Gradient) != 0.0f)
              {
                float ss = GaussianScale;
                float sm = pow(1.0f / ss, 1.0f / (float)(DIMENSION));
                Difference.Elements[0] *= ss;
                for (int ii = 1; ii < DIMENSION; ++ii)
                {
                  Difference.Elements[ii] *= sm;
                }
              }
              else
              {
                Difference = ODifference;
              }
   

              Difference = ODifference;
              for (int ii = 0; ii < DIMENSION; ++ii)
              {
                Difference.Elements[ii] /= (AverageDistance / 5.0f);
              }
 

              //Difference.Elements[0] = ODifference.Elements[0] * (float)this->width / this->KdTree.SampleExtent.Elements[0];
              //Difference.Elements[1] = ODifference.Elements[1] * (float)this->height / this->KdTree.SampleExtent.Elements[1];

              float cw = exp(-min(0.5f * Difference.DotProduct(Difference), 30.0f));
              //if ((Difference.Elements[0] * Difference.Elements[0] + Difference.Elements[1] * Difference.Elements[1]) < 1.0f)
              //{
              //  cw = 1.0f;
              //}
              //else
              //{
              //  cw = 0.0f;
              //}
              w += cw;
              La = La + SamplesLs[CellKNNs.SamplesIndices[k]] * cw;
            }
            La = La / (w + 1e-10f);
#else

            /// ----------- exact closest sample query --------------
					  //if (PrevIndex != -1) 
       //     {
       //       //float dx = QueryPosition.Elements[0] - this->KdTree.Samples[PrevIndex].Position.Elements[0];
       //       //float dy = QueryPosition.Elements[1] - this->KdTree.Samples[PrevIndex].Position.Elements[1];
       //       //SearchRadius = (dx * dx + dy * dy) * 1.01f;

       //       SearchRadius = QueryPosition.SquaredDistance(this->KdTree.Samples[PrevIndex].Position) * 1.1f;
       //     }
					  //PrevIndex = this->KdTree.NearestSearch(QueryPosition, BBox, SearchRadius);


            //------------- approximate closest sample query --------------
			      // generate a new position within the bounding box
			      //for (int ei = 0; ei < DIMENSION; ++ei)
			      //{
				     // QueryPosition.Elements[ei] = (this->Nodes[i].BBox.Max.Elements[ei] - this->Nodes[i].BBox.Min.Elements[ei]) * mtrand.rand() + this->Nodes[i].BBox.Min.Elements[ei];
			      //}
					  QueryPosition.Elements[1] = (y + 0.5f) * RevHeight;


            //this->KdTree.KNNSearch(this->CellKNNsNum, QueryPosition, &CellKNNs);

            float mind = 1E+20f;
            PrevIndex = 0;//this->KdTree.Nodes[i].SamplesIndices[0];
            for (int k = 0; k < KdTreeCellKNNsNum; ++k)
            { 
              int idx = CellKNNs.SamplesIndices[k];
              float val = 0.0f;//
              //val = QueryPosition.SquaredDistance(this->KdTree.Samples[idx].Position);
						  //float dx = QueryPosition.Elements[0] - this->KdTree.Samples[idx].Position.Elements[0];
						  //float dy = QueryPosition.Elements[1] - this->KdTree.Samples[idx].Position.Elements[1];
						  //float val = dx * dx + dy * dy;


 



              bool isNotNearest = false;
					    val = 0.0f;
					    CVector DifferenceVector;
					    DifferenceVector = this->KdTree.Samples[idx].Position - QueryPosition;
              //int DistanceIndex = k * this->CellKNNsNum;
              for (int mm = 0; mm < KdTreeCellKNNsNum; mm++)
					    {
                //float cdx = dx * this->KdTree.Samples[CellKNNs.SamplesIndices[mm]].Gradient.Elements[0];
                //float cdy = dy * this->KdTree.Samples[CellKNNs.SamplesIndices[mm]].Gradient.Elements[1];
                //float cval = BaseDistances[DistanceIndex] + cdx + cdy;
              
                //val += cval * cval;
                //DistanceIndex++;

						    float cval = DifferenceVector.DotProduct(this->KdTree.Samples[CellKNNs.SamplesIndices[mm]].Gradient);
						    val += cval * cval;
                if (val >= mind)
                {
                  isNotNearest = true;
                  break;
                }
              }
              if (isNotNearest)
              {
                continue;
              }






              //float val = 0.0f;
              //for (int mm = 0; mm < 2; mm++)
              //{
              //  float dx = QueryPosition.Elements[mm] - this->KdTree.Samples[idx].Position.Elements[mm];
              //  val = val + dx * dx;
              //}
              //for (int mm = 2; mm < DIMENSION; mm++)
              //{
              //  float dx = QueryPosition.Elements[mm] - this->KdTree.Samples[idx].Position.Elements[mm];
              //  val = val + dx * dx * ProjectedVariance[x + y * this->width];
              //}
              //for (int mm = 0; mm < DIMENSION; mm++)
              //{
              //  float dx = QueryPosition.Elements[mm] - this->KdTree.Samples[idx].Position.Elements[mm];
              //  val = val + dx * dx * ScalingFactor.Elements[mm] * ScalingFactor.Elements[mm];
              //}

              //CVector DifferenceVector;
              //DifferenceVector = this->KdTree.Samples[idx].Position - QueryPosition;
              //CVector TransformedVector;
              //for (int m = 0; m < DIMENSION; m++)
              //{
              //  TransformedVector.Elements[m] = 0.0f;
              //}
              //for (int m = 0; m < DIMENSION; m++)
              //{
              //  for (int mi = 0; mi < DIMENSION; mi++)
              //  {
              //    TransformedVector.Elements[m] += DifferenceMatrix[m][mi] * DifferenceVector.Elements[mi];
              //  }
              //}
              //val = fabs(TransformedVector.DotProduct(DifferenceVector));


 
              if (val < mind)
              {
                mind = val;
                PrevIndex = idx;
              }
            }
            //La = SamplesLs[PrevIndex];
#endif
            //splt.params[4] = 1.0f;
				    float dx = (QueryPosition.Elements[0] - this->KdTree.Samples[PrevIndex].Position.Elements[0]) / this->KdTree.SampleExtent.Elements[0] * (float)(this->width);
				    float dy = (QueryPosition.Elements[1] - this->KdTree.Samples[PrevIndex].Position.Elements[1]) / this->KdTree.SampleExtent.Elements[1] * (float)(this->height);
				    float val = (dx * dx + dy * dy);

            SamplePixelDistance[PrevIndex] += val;
          
	          Sample sample(x, y);
            film->AddSample( sample, unused,  SamplesLs[PrevIndex], 1,&splt);
            
            //bool isSamePixel = false;
            //for (int m = 0; m < (int)SamplePixelCount[PrevIndex].Px.size(); m++)
            //{
            //  if ((SamplePixelCount[PrevIndex].Px[m] == x) && (SamplePixelCount[PrevIndex].Py[m] == y))
            //  {
            //    isSamePixel = true;
            //    break;
            //  }
            //}
            //if (!isSamePixel)
            //{
              SamplePixelCount[PrevIndex].Px.push_back(x);
              SamplePixelCount[PrevIndex].Py.push_back(y);
            //}


            //Spectrum Ls(val);
            //film->AddSample( sample, unused, Ls, 1,&splt);

				  } 
        }


		  //}
       
        if ((i % 1000) == 0) printf("Reconstructing: [%d%%]\r", (int)((i / (float)this->KdTree.CurrentNumNodes) * 100.0f));

    }
    printf("Reconstructing the image: [%d%%]\r", 100);


	  film->WriteImage();
	


    // write out sample distribution
    string filename5 = "2dspread.exr";
	  ParamSet film_params5;

	  film_params5.AddInt("xresolution", &width);
	  film_params5.AddInt("yresolution", &height);
	  film_params5.AddString("filename", &filename5);
	  Film* film5 = MakeFilm("samoa_image", film_params5, 0);
    splt.params[0]=0;
    splt.params[1]=0;
    splt.params[2]=0;
    splt.params[3]=0;
    splt.params[4]=1.0f;
    WhiteL = 1.0f;
    for (int i = 0; i < this->KdTree.CurrentNumSamples; ++i)
    {
		  int x = (int)(this->KdTree.Samples[i].Position.Elements[0] * this->width / this->KdTree.SampleExtent.Elements[0]);
		  int y = (int)(this->KdTree.Samples[i].Position.Elements[1] * this->height / this->KdTree.SampleExtent.Elements[1]);

	    Ray unused;
	    Sample sample(x, y);
      //film2->AddSample( sample, unused, SamplesLs[i], 1,&splt);
      //film5->AddSample( sample, unused, SamplePixelCount[i].Px.size() / (float)(this->width * this->height), 1,&splt);
      film5->AddSample( sample, unused, SamplePixelDistance[i] / ((float)SamplePixelCount[i].Px.size() + 1e-10f), 1,&splt);
   }
    film5->WriteImage();





		//fprintf(sampfile, "\nFinalState\n");
		//fclose (sampfile);
	  //printf("\n");

		//  net->reconstruct();
		float reconTime = ((float)clock()-(float)startTime)/CLOCKS_PER_SEC;
		printf("\nReconstruction: %.3fsec\n",reconTime);

		return false;
	}



  if (this->sampleNum <= this->InitialSampleNum)
  {
    sample->imageX = (float)(Rand.rand()) * (float)(this->width) + this->xPixelStart + 2;
    sample->imageY = (float)(Rand.rand()) * (float)(this->height) + this->yPixelStart + 2;
    if (DIMENSION == 2) 
    {
      sample->time = 0.0f;
      sample->lensU = 0.0f;
      sample->lensV = 0.0f;
      this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
      this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
    }
    else if (DIMENSION == 3) 
    {
      sample->time = (float)(Rand.rand());
      sample->lensU = 0.0f;
      sample->lensV = 0.0f;
      this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
      this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
      this->CurrentSample.Position.Elements[2] = sample->time * this->KdTree.SampleExtent.Elements[2];
    }
    else if (DIMENSION == 5)
    {
      // dof and motion blur
      sample->time = (float)(Rand.rand());
      sample->lensU = (float)(Rand.rand());
      sample->lensV = (float)(Rand.rand());
      this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
      this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
      this->CurrentSample.Position.Elements[2] = sample->time * this->KdTree.SampleExtent.Elements[2];
      this->CurrentSample.Position.Elements[3] = sample->lensU * this->KdTree.SampleExtent.Elements[3];
      this->CurrentSample.Position.Elements[4] = sample->lensV * this->KdTree.SampleExtent.Elements[4];
    }
    else if (DIMENSION == 6)
    {
      // softshadow and dof
      sample->time = 0.0f;
      sample->twoD[0][0] = (float)(Rand.rand());
      sample->twoD[0][1] = (float)(Rand.rand());
      sample->lensU = (float)(Rand.rand());
      sample->lensV = (float)(Rand.rand());
      this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
      this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
      this->CurrentSample.Position.Elements[2] = sample->twoD[0][0] * this->KdTree.SampleExtent.Elements[2];
      this->CurrentSample.Position.Elements[3] = sample->twoD[0][1] * this->KdTree.SampleExtent.Elements[3];
      this->CurrentSample.Position.Elements[4] = sample->lensU * this->KdTree.SampleExtent.Elements[4];
      this->CurrentSample.Position.Elements[5] = sample->lensV * this->KdTree.SampleExtent.Elements[5];
    }

    else
    {
      sample->time = 0.0f;

#ifdef AREA_LIGHT
      sample->lensU = 0.0f;
      sample->lensV = 0.0f;
      sample->twoD[0][0] = (float)(Rand.rand());
      sample->twoD[0][1] = (float)(Rand.rand());
#else
      sample->lensU = (float)(Rand.rand());
      sample->lensV = (float)(Rand.rand());
#endif

      this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
      this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];

#ifdef AREA_LIGHT
      this->CurrentSample.Position.Elements[2] = sample->twoD[0][0] * this->KdTree.SampleExtent.Elements[2];
      this->CurrentSample.Position.Elements[3] = sample->twoD[0][1] * this->KdTree.SampleExtent.Elements[3];
#else
      this->CurrentSample.Position.Elements[2] = sample->lensU * this->KdTree.SampleExtent.Elements[2];
      this->CurrentSample.Position.Elements[3] = sample->lensV * this->KdTree.SampleExtent.Elements[3];
#endif

    }
 
    this->InitialSamples.push_back(this->CurrentSample);
  }



  else
  {
    CVector NewPosition;

    NewPosition = this->KdTree.GenerateNewSamplePosition(this->Rand, false);

		// ******** uniform sampling ********
		for (int i = 0; i < DIMENSION; i++)
		{
			NewPosition.Elements[i] = Rand.rand() * this->KdTree.SampleExtent.Elements[i];
		}
		// ******** uniform sampling ********


    sample->imageX = NewPosition.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width + this->xPixelStart + 2.0f;
    sample->imageY = NewPosition.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height + this->yPixelStart + 2.0f;
    if (DIMENSION == 2) 
    {
      sample->time = 0.0f;
      sample->lensU = 0.0f;
      sample->lensV = 0.0f;
    }
    else if (DIMENSION == 3) 
    {
      sample->time = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
      sample->lensU = 0.0f;
      sample->lensV = 0.0f;
    }
    else if (DIMENSION == 5)
    {
      sample->time = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
      sample->lensU = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
      sample->lensV = NewPosition.Elements[4] / this->KdTree.SampleExtent.Elements[4];
    }
    else if (DIMENSION == 6)
    {
      sample->time = 0.0f;
      sample->twoD[0][0] = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
      sample->twoD[0][1] = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
      sample->lensU = NewPosition.Elements[4] / this->KdTree.SampleExtent.Elements[4];
      sample->lensV = NewPosition.Elements[5] / this->KdTree.SampleExtent.Elements[5];
    }
    else
    {
      sample->time = 0.0f;
#ifdef AREA_LIGHT
      sample->twoD[0][0] = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
      sample->twoD[0][1] = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
#else
      sample->lensU = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
      sample->lensV = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
#endif
    }

    this->CurrentSample.Position = NewPosition;
  }


	finishTime = clock();
	elapsedTimeGet = ((float)finishTime - (float)startTime) / CLOCKS_PER_SEC;
	return true;
}



 

Spectrum Samoa2Sampler::SetPrevSample(Spectrum Ls) 
{
	//net->setPrevSample(Ls);
  if (this->sampleNum <= this->InitialSampleNum)
  {
    //float xyz[3];
    //Ls.XYZ(xyz);
    Ls.XYZ(this->InitialSamples[this->InitialSamples.size() - 1].XYZ);
    //this->InitialSamples[this->InitialSamples.size() - 1].Luminance = this->LuminanceScale * Ls.y();

    if (this->sampleNum == this->InitialSampleNum)
    {

      this->KdTree.Build(this->InitialSamples);
#ifdef CELLSPLIT
      for (int i = 0; i < (int)this->KdTree.CurrentNumNodes; i++)
      {
        if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) 
        {
          continue;
        }

        this->KdTree.ComputeNodeError(i, false);
        this->KdTree.Heap.Put(i, this->KdTree.Nodes[i].Error); 
      }
#else 
	    this->KdTree.ComputeKNNs();
	    this->KdTree.ComputeBSpheres();
      //for (int i = 0; i < this->InitialSampleNum; i++)
      //{
      //  this->KdTree.ComputeSampleGradient(i);
      //} 
      for (int i = 0; i < this->InitialSampleNum; i++)
      {
        this->KdTree.ComputeSampleDistance(i);
      } 
      for (int i = 0; i < this->InitialSampleNum; i++)
      {
        this->KdTree.Heap.Put(i, this->KdTree.Samples[i].Distance); 
      } 
#endif

    }
  }
  else
  {
    if (this->sampleNum <= this->totalSamples)
    {
      //this->CurrentSample.Luminance = this->LuminanceScale * Ls.y();
      Ls.XYZ(this->CurrentSample.XYZ);

      this->KdTree.Insert(this->CurrentSample, false);
    }

    //if (this->sampleNum == this->NumAdaptiveSamples)
    //{
    //  this->KdTree.Heap.Clear();
    //  for (int i = 0; i < (int)this->KdTree.Nodes.size(); i++)
    //  {
    //    if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) 
    //    {
    //      continue;
    //    }
    //    this->KdTree.ComputeNodeError(i, true);
    //    this->KdTree.Heap.Put(i, this->KdTree.Nodes[i].Error); 
    //  }
    //}
  }
  this->SamplesLs[this->sampleNum - 1] = Ls;

	return Ls;
}






extern "C" DLLEXPORT Sampler *CreateSampler(const ParamSet &params, const Film *film) 
{
  std::cout << "DIMENSIONS: " << DIMENSION << std::endl;
  string reconfile = params.FindOneString("reconfile", "recon.exr");
  string densityfile = params.FindOneString("densityfile", "density.exr");
	int width = params.FindOneInt("width", 500);
	int height = params.FindOneInt("height", 500);

	// initialize common sampler parameters
	int xstart, xend, ystart, yend;
	film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
	int ps = params.FindOneInt("pixelsamples", 2); 
	//return new Samoa2Sampler(xstart, xend, ystart, yend, ps, width, height,samplefile,imagefile);

  int candidatenum = params.FindOneInt("numcandidate", 4); 
  int knnsnum = params.FindOneInt("numknns", 4);
  int bucketnum = params.FindOneInt("numcellsmp", 4); 
  int cellknnsnum = params.FindOneInt("numreconknns", 15);
  int initnum = params.FindOneInt("numinit", 1024);
  float imagescale = params.FindOneFloat("scaleimageaxes", 1.0f);
  float luminancelimit = params.FindOneFloat("limitluminance", 100.0f); 
  float luminancescale = params.FindOneFloat("scaleluminance", 40.0f);
  float distanceepsilon = params.FindOneFloat("epsilondistance", 1E-5f);
  float maxgs = params.FindOneFloat("maxgaussianscale", 10.0f);
  int optstep = params.FindOneInt("optimizationstep", 10);

  return new Samoa2Sampler(candidatenum, knnsnum, bucketnum, cellknnsnum, maxgs, optstep, initnum, imagescale, luminancelimit, luminancescale, distanceepsilon, xstart, xend, ystart, yend, ps, width, height, reconfile, densityfile);
}
