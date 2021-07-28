/*
 * samoa2.h
 *
 *  Created on: Sep 3, 2019
 *      Author: miguel
 */

#ifndef SRC_INTEGRATORS_MDAS_SAMOA2_H_
#define SRC_INTEGRATORS_MDAS_SAMOA2_H_


// Samoa2
#include <ctime>
#include <iostream>

#include "mt.h"
#include "vector.h"
#include "kdtree.h"
#include <vector>
#include <fstream>
#include <iostream>


#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/render/film.h>
#include "sample.h"
#include "samoa_image.h"

extern "C"
{
	//#include "f2c.h"
	//#include "clapack.h"
}

/*struct Splat
{
	float params[5];
	mitsuba::Spectrum L0, L1, L2, L3;
};*/


struct PixelCount
{
	std::vector<int>Px;
	std::vector<int>Py;
};



// #define AREA_LIGHT



#define VORONOI




template <int DIMENSION, bool AREA_LIGHT>
class Samoa2Sampler
{
public:
	Samoa2Sampler(int candidatenum, int knnsnum, int bucketnum, int cellknnsnum, float maxgs, int optstep, int initnum,
			      float imagescale, float luminancelimit, float luminancescale, float distanceepsilon, int xstart, int xend,
				  int ystart, int yend, int ps, int width, int height, std::string reconfile, std::string densityfile);
	~Samoa2Sampler(){}

	int RoundSize(int size) const;
	bool GetNextSample(Sample *sample);
	mitsuba::Spectrum SetPrevSample(mitsuba::Spectrum Ls);
	void sampleKNN(Sample *sample);
	void sampleInitial(Sample *sample);
	void reconstructImage();

	void generateSampleDistribution(mitsuba::ImageBlock *block);

	CKdTree<DIMENSION> KdTree;
	MTRand Rand;
	std::vector<CSample<DIMENSION>> InitialSamples;
	std::vector<mitsuba::Spectrum> SamplesLs;
	int InitialSampleNum;
	CSample<DIMENSION> CurrentSample;

	int NumAdaptiveSamples;
	float LuminanceScale;
	int CellKNNsNum;
	float MaxGaussianScale;
	int NumOptimizationStep;
	float Anisotropy;

	int totalSamples;

	// Films buffers
	SamoaImageFilm samoaFilm;

private:
	int width;
	int height;

	// Added from Sampler implementation
	int xStart;
	int xEnd;
	int yStart;
	int yEnd;
	int ps;
	int xPixelStart;
	int xPixelEnd;
	int yPixelStart;
	int yPixelEnd;
	std::string reconfile;
	std::string densityfile;

	int initialSamples;
	int sampleNum;
	FILE * sampfile;

	clock_t startTime,finishTime;
	float elapsedTimeGet, elapsedTimeSet;


};

// template <int DIMENSION, bool AREA_LIGHT>
// Samoa2Sampler<DIMENSION, AREA_LIGHT> createSamoaSampler (const mitsuba::Properties &params, mitsuba::Film* film);


#include "samoa2.cpp"


#endif /* SRC_INTEGRATORS_MDAS_SAMOA2_H_ */
