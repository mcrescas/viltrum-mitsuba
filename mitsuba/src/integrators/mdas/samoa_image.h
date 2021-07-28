/*
 * samoa_image.h
 *
 *  Created on: Sep 3, 2019
 *      Author: miguel
 */

#ifndef SRC_INTEGRATORS_MDAS_SAMOA_IMAGE_H_
#define SRC_INTEGRATORS_MDAS_SAMOA_IMAGE_H_

#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/imageblock.h>
#include <mitsuba/core/ref.h>
#include <mitsuba/core/properties.h>

#include "sample.h"

struct Splat
{
  float params[5];
  mitsuba::Spectrum L0, L1, L2, L3;
};


// SamoaImageFilm Declarations
class SamoaImageFilm
{
  public:
	  // SamoaImageFilm Public Methods
	  SamoaImageFilm(int xres, int yres, const float crop[4], bool premult);
	  SamoaImageFilm() {};
	  ~SamoaImageFilm() {
		  //delete pixels;
		  //delete isSampled;
	  }

	  //void SamoaImageFilm::AddSample(float imageX, float imageY, float * params, const mitsuba::Spectrum &L, float alpha);
	  void AddSample(float imageX, float imageY, float * params, const mitsuba::Spectrum &L, float alpha, const mitsuba::Spectrum L0, const mitsuba::Spectrum L1, const mitsuba::Spectrum L2, const mitsuba::Spectrum L3);
	  //void AddSample(const Sample &sample, const Ray &ray, const mitsuba::Spectrum &L, float alpha) {Error("Add Sample requires more parameters.");}
	  void AddSample(const Sample &sample, const mitsuba::Spectrum &L, float alpha, void* data);
	  void GetSampleExtent(int *xstart, int *xend, int *ystart, int *yend) const;
	  void Clear();
	  void WriteImage(mitsuba::ImageBlock *imageBlock);



  private:
	  // SamoaImageFilm Private Data
	  bool premultiplyAlpha;
	  float cropWindow[4];
	  int xPixelStart, yPixelStart, xPixelCount, yPixelCount;

	  int xResolution , yResolution;

	  struct Pixel
	  {
		  Pixel() : L(0.0f)
      {
			  alpha = 0.0f;
			  weightSum = 0.0f;
		  }
		  mitsuba::Spectrum L;
		  float alpha;
		  float weightSum;
	  };

	  //BlockedArray<Pixel> *pixels;
	  //BlockedArray<bool> *isSampled;
	  std::vector<std::vector<Pixel>> pixels;
	  std::vector<std::vector<bool>> isSampled;
};

SamoaImageFilm createSamoaFilm(int w, int h, float crop[4]);


#endif /* SRC_INTEGRATORS_MDAS_SAMOA_IMAGE_H_ */
