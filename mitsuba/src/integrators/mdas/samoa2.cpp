

#include "samoa2.h"


template <int DIMENSION, bool AREA_LIGHT>
int Samoa2Sampler<DIMENSION, AREA_LIGHT>::RoundSize(int size) const 
{
	return size;
}

template <int DIMENSION, bool AREA_LIGHT>
Samoa2Sampler<DIMENSION, AREA_LIGHT>::Samoa2Sampler(int candidatenum, int knnsnum, int bucketnum, int cellknnsnum, float maxgs, int optstep, // @suppress("Class members should be properly initialized")
							 int initnum, float imagescale, float luminancelimit, float luminancescale, float distanceepsilon,
							 int xstart, int xend, int ystart, int yend, int ps, int width, int height, std::string reconfile,
							 std::string densityfile)
	: reconfile(reconfile), densityfile(densityfile)
{


	this->xStart = xstart;
	this->xEnd = xend;
	this->yStart = ystart;
	this->yEnd = yend;

	this->ps = ps;


	this->xPixelStart = xStart;
	this->xPixelEnd = xEnd;
	this->yPixelStart = yStart;
	this->yPixelEnd = yEnd;

	this->sampleNum = 0;
	this->totalSamples = ps * (xPixelEnd - xPixelStart) * (yPixelEnd - yPixelStart);

	this->KdTree.CandidatesNum = candidatenum;
	this->KdTree.KNNsNum = knnsnum;
	this->KdTree.BucketSize = bucketnum;
	this->InitialSampleNum = initnum;
	this->initialSamples = initnum;

	//this->width = (this->xPixelEnd - this->xPixelStart - 4);
	//this->height = (this->yPixelEnd - this->yPixelStart - 4);
	this->width = (this->xPixelEnd - this->xPixelStart);
	this->height = (this->yPixelEnd - this->yPixelStart);

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
	std::cout << "total: " << this->totalSamples  << std::endl;
	std::cout << "w: " << width  << std::endl;
	std::cout << "h: " << height << std::endl;

	Rand.seed(570905);

	// allocate all memory
	this->KdTree.Nodes.resize(this->totalSamples * 3);
	this->KdTree.Samples.resize(this->totalSamples);
	this->KdTree.Heap.Items.resize(this->totalSamples * 2 + 1);
	this->KdTree.Heap.Indices.resize(this->totalSamples * 2 + 1);
	this->SamplesLs.resize(this->totalSamples);
}




template <int DIMENSION, bool AREA_LIGHT>
bool Samoa2Sampler<DIMENSION, AREA_LIGHT>::GetNextSample(Sample *sample) 
{
	startTime = clock();
	this->sampleNum++;


	if (this->sampleNum == (this->totalSamples)) 
	{
		this->reconstructImage();
		return false;
	}



	if (this->sampleNum <= this->InitialSampleNum)
	{
		this->sampleInitial(sample);
	}
	else
	{
		this->sampleKNN(sample);
	}


	finishTime = clock();
	elapsedTimeGet = ((float)finishTime - (float)startTime) / CLOCKS_PER_SEC;
	return true;
}




template <int DIMENSION, bool AREA_LIGHT>
mitsuba::Spectrum Samoa2Sampler<DIMENSION, AREA_LIGHT>::SetPrevSample(mitsuba::Spectrum Ls)
{
	//net->setPrevSample(Ls);
	if (this->sampleNum <= this->InitialSampleNum)
	{
		//Ls.XYZ(this->InitialSamples[this->InitialSamples.size() - 1].XYZ);
		mitsuba::Float xyz_data[3];
		Ls.toXYZ(xyz_data[0], xyz_data[1], xyz_data[2]);
		this->InitialSamples[this->InitialSamples.size() - 1].XYZ[0] = xyz_data[0];
		this->InitialSamples[this->InitialSamples.size() - 1].XYZ[1] = xyz_data[1];
		this->InitialSamples[this->InitialSamples.size() - 1].XYZ[2] = xyz_data[2];

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
			//Ls.XYZ(this->CurrentSample.XYZ);

			mitsuba::Float xyz__[3];
			Ls.toXYZ(xyz__[0], xyz__[1], xyz__[2]);

			this->CurrentSample.XYZ[0] = xyz__[0];
			this->CurrentSample.XYZ[1] = xyz__[1];
			this->CurrentSample.XYZ[2] = xyz__[2];

			this->KdTree.Insert(this->CurrentSample, false);
		}
	}
	this->SamplesLs[this->sampleNum - 1] = Ls;

	return Ls;
}


template <int DIMENSION, bool AREA_LIGHT>
void Samoa2Sampler<DIMENSION, AREA_LIGHT>::sampleInitial(Sample* sample) {
	//sample->imageX = (float)(Rand.rand()) * (float)(this->width) + this->xPixelStart + 2;
	//sample->imageY = (float)(Rand.rand()) * (float)(this->height) + this->yPixelStart + 2;
	sample->imageX = (float)(Rand.rand()) * (float)(this->width);
	sample->imageY = (float)(Rand.rand()) * (float)(this->height);
	if (DIMENSION == 2)
	{
		sample->time = 0.0f;
		sample->lensU = 0.0f;
		sample->lensV = 0.0f;
		//this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart) / (float)this->height * this->KdTree.SampleExtent.Elements[1];


		//this->CurrentSample.Position.Elements[0] = (sample->imageX) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
	}
	else if (DIMENSION == 3)
	{
		sample->time = (float)(Rand.rand());
		sample->lensU = 0.0f;
		sample->lensV = 0.0f;
		//this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[2] = sample->time * this->KdTree.SampleExtent.Elements[2];
	}
	else if (DIMENSION == 5)
	{
		if (AREA_LIGHT) {
			sample->time = (float)(Rand.rand());
			sample->twoD[0][0] = (float)(Rand.rand());
			sample->twoD[0][1] = (float)(Rand.rand());
		} else {
			// dof and motion blur
			sample->time = (float)(Rand.rand());
			sample->lensU = (float)(Rand.rand());
			sample->lensV = (float)(Rand.rand());
		}

		//this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[2] = sample->time * this->KdTree.SampleExtent.Elements[2];

		if (AREA_LIGHT) {
			this->CurrentSample.Position.Elements[3] = sample->twoD[0][0] * this->KdTree.SampleExtent.Elements[3];
			this->CurrentSample.Position.Elements[4] = sample->twoD[0][1] * this->KdTree.SampleExtent.Elements[4];
		} else {
			this->CurrentSample.Position.Elements[3] = sample->lensU * this->KdTree.SampleExtent.Elements[3];
			this->CurrentSample.Position.Elements[4] = sample->lensV * this->KdTree.SampleExtent.Elements[4];
		}
	}
	else if (DIMENSION == 6)
	{
		// softshadow and dof
		sample->time = 0.0f;
		sample->twoD[0][0] = (float)(Rand.rand());
		sample->twoD[0][1] = (float)(Rand.rand());
		sample->lensU = (float)(Rand.rand());
		sample->lensV = (float)(Rand.rand());
		//this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[2] = sample->twoD[0][0] * this->KdTree.SampleExtent.Elements[2];
		this->CurrentSample.Position.Elements[3] = sample->twoD[0][1] * this->KdTree.SampleExtent.Elements[3];
		this->CurrentSample.Position.Elements[4] = sample->lensU * this->KdTree.SampleExtent.Elements[4];
		this->CurrentSample.Position.Elements[5] = sample->lensV * this->KdTree.SampleExtent.Elements[5];
	}
	else
	{
		sample->time = 0.0f;

		if (AREA_LIGHT) {
			sample->lensU = 0.0f;
			sample->lensV = 0.0f;
			sample->twoD[0][0] = (float)(Rand.rand());
			sample->twoD[0][1] = (float)(Rand.rand());
		} else {
			sample->lensU = (float)(Rand.rand());
			sample->lensV = (float)(Rand.rand());
		}

		//this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart - 2) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		//this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart - 2) / (float)this->height * this->KdTree.SampleExtent.Elements[1];
		this->CurrentSample.Position.Elements[0] = (sample->imageX - this->xPixelStart) / (float)this->width * this->KdTree.SampleExtent.Elements[0];
		this->CurrentSample.Position.Elements[1] = (sample->imageY - this->yPixelStart) / (float)this->height * this->KdTree.SampleExtent.Elements[1];

		if (AREA_LIGHT) {
			this->CurrentSample.Position.Elements[2] = sample->twoD[0][0] * this->KdTree.SampleExtent.Elements[2];
			this->CurrentSample.Position.Elements[3] = sample->twoD[0][1] * this->KdTree.SampleExtent.Elements[3];
		} else {
			this->CurrentSample.Position.Elements[2] = sample->lensU * this->KdTree.SampleExtent.Elements[2];
			this->CurrentSample.Position.Elements[3] = sample->lensV * this->KdTree.SampleExtent.Elements[3];
		}

	}

	this->InitialSamples.push_back(this->CurrentSample);
}


/*
 * Get a new sample given KNN
 *
 */
template <int DIMENSION, bool AREA_LIGHT>
void Samoa2Sampler<DIMENSION, AREA_LIGHT>::sampleKNN(Sample* sample) {
	CVector<DIMENSION> NewPosition;

	NewPosition = this->KdTree.GenerateNewSamplePosition(this->Rand, false);

	// ******** uniform sampling ********
	// for (int i = 0; i < DIMENSION; i++)
	// {
	// 	NewPosition.Elements[i] = Rand.rand() * this->KdTree.SampleExtent.Elements[i];
	// }
	// ******** uniform sampling ********


	//sample->imageX = NewPosition.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width + this->xPixelStart + 2.0f;
	//sample->imageY = NewPosition.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height + this->yPixelStart + 2.0f;
	sample->imageX = NewPosition.Elements[0] / this->KdTree.SampleExtent.Elements[0] * (float)this->width + this->xPixelStart;
	sample->imageY = NewPosition.Elements[1] / this->KdTree.SampleExtent.Elements[1] * (float)this->height + this->yPixelStart;
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
		if (AREA_LIGHT) {
			sample->time = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
			sample->twoD[0][0] = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
			sample->twoD[0][1] = NewPosition.Elements[4] / this->KdTree.SampleExtent.Elements[4];
		} else {
			sample->time = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
			sample->lensU = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
			sample->lensV = NewPosition.Elements[4] / this->KdTree.SampleExtent.Elements[4];
		}
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
		if (AREA_LIGHT) {
			sample->twoD[0][0] = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
			sample->twoD[0][1] = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
		} else {
			sample->lensU = NewPosition.Elements[2] / this->KdTree.SampleExtent.Elements[2];
			sample->lensV = NewPosition.Elements[3] / this->KdTree.SampleExtent.Elements[3];
		}
	}

	this->CurrentSample.Position = NewPosition;
}

template <int DIMENSION, bool AREA_LIGHT>
void Samoa2Sampler<DIMENSION, AREA_LIGHT>::reconstructImage() {
	long recon_start_time = clock();

	std::vector<CHeapItem>().swap(this->KdTree.Heap.Items);
	std::vector<int>().swap(this->KdTree.Heap.Indices);

	/*char *searchpath = getenv("PBRT_SEARCHPATH");
	if (!searchpath)
	{
		std::cerr << "You need to set PBRT_SEARCHPATH." << std::endl;
		return false;
	}
	UpdatePluginPath(searchpath);*/

	//int w = (int)(width);
	//int h = (int)(height);
	//samoaFilm = createSamoaFilm(w, h);

	//float params[5];
	Splat splt;



	// write out sample distribution
	/*ParamSet film_params2;

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
		film2->AddSample( sample, unused, WhiteL, 1,&splt);
	}
	film2->WriteImage();*/




	//printf("\n");
	this->KdTree.BucketSize = 1;
	int CurrenNumtKdtreeNodes = (int)this->KdTree.CurrentNumNodes;
	for (int i = 0; i < CurrenNumtKdtreeNodes; ++i)
	{
		if (this->KdTree.Nodes[i].NodeType != NodeTypeLeaf) continue;
		if (this->KdTree.Nodes[i].SamplesIndices.size() == 0) continue;

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


	// write out sample distribution
	/*string filename3 = "aniso.exr";
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


		Spectrum Ls(MaxScaling / (MinScaling + 1e-20f));

		Ray unused;
		Sample sample(x, y);
		film3->AddSample(sample, unused, Ls, 1, &splt);
	}
	film3->WriteImage();*/



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

		CVector<DIMENSION> QueryPosition;
		//CKNNs KNNs;

		// splatting
		float SearchRadius = 1E+20f;
		int PrevIndex = -1;
		float RevWidth = this->KdTree.SampleExtent.Elements[0] / (float)(this->width);
		float RevHeight = this->KdTree.SampleExtent.Elements[1] / (float)(this->height);

		mitsuba::Spectrum La;

		QueryPosition = 0.5f * (this->KdTree.Nodes[i].BBox.Min + this->KdTree.Nodes[i].BBox.Max);

		CKNNs &CellKNNs = this->KdTree.Samples[this->KdTree.Nodes[i].SamplesIndices[0]].KNNs;



		for (int x = ix0; x < ix1; ++x)
		{
			QueryPosition.Elements[0] = (x + 0.5f) * RevWidth;
			//if ((x >= this->width) || (x < 0)) continue;
			for (int y = iy0; y < iy1; ++y)
			{
				//if ((y >= this->height) || (y < 0)) continue;

				QueryPosition.Elements[1] = (y + 0.5f) * RevHeight;


				float mind = 1E+20f;
				PrevIndex = 0;
				for (int k = 0; k < KdTreeCellKNNsNum; ++k)
				{
					int idx = CellKNNs.SamplesIndices[k];
					float val = 0.0f;


					bool isNotNearest = false;
					val = 0.0f;
					CVector<DIMENSION> DifferenceVector;
					DifferenceVector = this->KdTree.Samples[idx].Position - QueryPosition;
					//int DistanceIndex = k * this->CellKNNsNum;
					for (int mm = 0; mm < KdTreeCellKNNsNum; mm++)
					{

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



					if (val < mind)
					{
						mind = val;
						PrevIndex = idx;
					}
				}
				float dx = (QueryPosition.Elements[0] - this->KdTree.Samples[PrevIndex].Position.Elements[0]) / this->KdTree.SampleExtent.Elements[0] * (float)(this->width);
				float dy = (QueryPosition.Elements[1] - this->KdTree.Samples[PrevIndex].Position.Elements[1]) / this->KdTree.SampleExtent.Elements[1] * (float)(this->height);
				float val = (dx * dx + dy * dy);

				SamplePixelDistance[PrevIndex] += val;

				Sample sample(x, y);
				samoaFilm.AddSample( sample,  SamplesLs[PrevIndex], 1.0, &splt);
				//imageblock->put(mitsuba::Point2(x,y), SamplesLs[PrevIndex], 1.0);

				SamplePixelCount[PrevIndex].Px.push_back(x);
				SamplePixelCount[PrevIndex].Py.push_back(y);

			}
		}


		//}

		if ((i % 1000) == 0) printf("Reconstructing: [%d%%]\r", (int)((i / (float)this->KdTree.CurrentNumNodes) * 100.0f));

	}
	printf("Reconstructing the image: [%d%%]\r", 100);



	// write out sample distribution
	/*string filename5 = "2dspread.exr";
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

		film5->AddSample( sample, unused, SamplePixelDistance[i] / ((float)SamplePixelCount[i].Px.size() + 1e-10f), 1,&splt);
	}
	film5->WriteImage();*/

	float reconTime = ((float)clock()-(float)startTime)/CLOCKS_PER_SEC;
	printf("\nReconstruction: %.3fsec\n",reconTime);
}

template <int DIMENSION, bool AREA_LIGHT>
void Samoa2Sampler<DIMENSION, AREA_LIGHT>::generateSampleDistribution(mitsuba::ImageBlock *block) {
	for (int i = 0; i < this->KdTree.CurrentNumSamples; ++i)
	{
		int x = (int)(this->KdTree.Samples[i].Position.Elements[0] * this->width / this->KdTree.SampleExtent.Elements[0]);
		int y = (int)(this->KdTree.Samples[i].Position.Elements[1] * this->height / this->KdTree.SampleExtent.Elements[1]);

		mitsuba::Point2i pos(x,y);
		auto pixelValue = block->getBitmap()->getPixel(pos);
		pixelValue += mitsuba::Spectrum(1);
		block->getBitmap()->setPixel(pos, pixelValue);
	}
}



// template <int DIMENSION, bool AREA_LIGHT>
// Samoa2Sampler<DIMENSION, AREA_LIGHT> createSamoaSampler (const mitsuba::Properties &params, mitsuba::Film* film) {

// 	int width = film->getSize().x;
// 	int height = film->getSize().y;

// 	int xstart = 0, xend = width, ystart = 0, yend = height;
// 	//film->GetSampleExtent(&xstart, &xend, &ystart, &yend);
// 	int ps = params.getInteger("pixelsamples", 2);

// 	std::string reconfile = params.getString("reconfile", "recon.exr");
// 	std::string densityfile = params.getString("densityfile", "density.exr");

// 	int candidatenum = params.getInteger("numcandidate", 4);
// 	int knnsnum = params.getInteger("numknns", 4);
// 	int bucketnum = params.getInteger("numcellsmp", 4);
// 	int cellknnsnum = params.getInteger("numreconknns", 15);
// 	int initnum = params.getInteger("numinit", 1024);
// 	float imagescale = params.getFloat("scaleimageaxes", 1.0f);
// 	float luminancelimit = params.getFloat("limitluminance", 100.0f);
// 	float luminancescale = params.getFloat("scaleluminance", 40.0f);
// 	float distanceepsilon = params.getFloat("epsilondistance", 1E-5f);
// 	float maxgs = params.getFloat("maxgaussianscale", 10.0f);
// 	int optstep = params.getInteger("optimizationstep", 10);

// 	return Samoa2Sampler<DIMENSION, AREA_LIGHT>(candidatenum, knnsnum, bucketnum, cellknnsnum, maxgs, optstep, initnum, imagescale, luminancelimit, luminancescale, distanceepsilon, xstart, xend, ystart, yend, ps, width, height, reconfile, densityfile);
// }

