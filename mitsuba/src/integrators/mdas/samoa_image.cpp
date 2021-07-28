
#include "samoa_image.h"

// SamoaImageFilm Method Definitions
SamoaImageFilm::SamoaImageFilm(int xres, int yres, const float crop[4], bool premult)
{
  memcpy(cropWindow, crop, 4 * sizeof(float));

  premultiplyAlpha = premult;

  xResolution = xres;
  yResolution = yres;

  // Compute film image extent
  xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
  xPixelCount = std::max(1, Ceil2Int(xResolution * cropWindow[1]) - xPixelStart);
  yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
  yPixelCount = std::max(1, Ceil2Int(yResolution * cropWindow[3]) - yPixelStart);

  // Allocate film image storage
  //pixels = new BlockedArray<Pixel>(xPixelCount, yPixelCount);
  //isSampled = new BlockedArray<bool>(xPixelCount, yPixelCount);

  pixels.resize(xPixelCount);
  for(auto &arr : pixels) arr.resize(yPixelCount);
  isSampled.resize(xPixelCount);
  for(auto &arr : isSampled) arr.resize(yPixelCount);
}





void SamoaImageFilm::AddSample(float imageX, float imageY, float * params, const mitsuba::Spectrum &L, float alpha, const mitsuba::Spectrum L0, const mitsuba::Spectrum L1, const mitsuba::Spectrum L2, const mitsuba::Spectrum L3)
{
  
  int x = (int)imageX;
  int y = (int)imageY;
  x = std::max(x, this->xPixelStart);
  x = std::min(x, this->xPixelStart + this->xPixelCount - 1);
  y = std::max(y, this->yPixelStart);
  y = std::min(y, this->yPixelStart + this->yPixelCount - 1);

  // Update pixel values with filtered sample contribution
  //Pixel &pixel = (*pixels)(x - xPixelStart, y - yPixelStart);
  Pixel &pixel = pixels[x - xPixelStart][y - yPixelStart];

  pixel.alpha += alpha * params[4];
  pixel.weightSum += params[4];   
  //pixel.L.AddWeighted(params[4], L);
  pixel.L.addWeighted(params[4], L);
}





void SamoaImageFilm::AddSample(const Sample &sample, const mitsuba::Spectrum &L, float alpha, void* data)
{
	Splat* splt = (Splat*)data;
	//float* sampleData = (float*)data;
	AddSample(sample.imageX, sample.imageY, &(splt->params[0]), L, alpha, splt->L0, splt->L1, splt->L2, splt->L3);
}





void SamoaImageFilm::Clear() 
{
  for (int y = 0; y < yPixelCount; ++y) 
  {
    for (int x = 0; x < xPixelCount; ++x)
    {
		  //Pixel &pixel = (*pixels)(x, y);
    	  Pixel &pixel = pixels[x][y];
		  pixel.L = mitsuba::Spectrum();
		  pixel.alpha = 1.0f;
		  pixel.weightSum = 0;
      
		  //(*isSampled)(x, y) = false;
		  isSampled[x][y] = false;
    } 
	}
}





void SamoaImageFilm::GetSampleExtent(int *xstart, int *xend, int *ystart, int *yend) const 
{
	*xstart = Floor2Int(xPixelStart + .5f - 3.f);
	*xend   = Floor2Int(xPixelStart + .5f + xPixelCount + 3.f);
	*ystart = Floor2Int(yPixelStart + .5f - 3.f);
	*yend   = Floor2Int(yPixelStart + .5f + yPixelCount + 3.f);
}



void SamoaImageFilm::WriteImage(mitsuba::ImageBlock *imageblock)
{
	// Convert image to RGB and compute final pixel values
	int nPix = xPixelCount * yPixelCount;
	float *rgb = new float[3*nPix];
	float *alpha = new float[nPix];
	int offset = 0;

	float MaxWeight = 0.0f;
	for (int y = 0; y < yPixelCount; ++y)
	{
		for (int x = 0; x < xPixelCount; ++x)
		{
			//float weightSum = (*pixels)(x, y).weightSum;
			float weightSum = pixels[x][y].weightSum;
			if (weightSum > MaxWeight)
			{
				MaxWeight = weightSum;
			}
		}
	}

	printf("%s\n", imageblock->toString().c_str());

	for (int y = 0; y < yPixelCount; ++y)
	{

		for (int x = 0; x < xPixelCount; ++x)
		{
			// Convert pixel spectral radiance to RGB
			mitsuba::Float xyz[3];
			//(*pixels)(x, y).L.XYZ(xyz);

			//(pixels[x][y]).L.toXYZ(xyz[0], xyz[1], xyz[2]);
			mitsuba::Spectrum l = pixels[x][y].L;
			l.toXYZ(xyz[0],xyz[1],xyz[2]);

			const float rWeight[3] = { 3.240479f, -1.537150f, -0.498535f };
			const float gWeight[3] = {-0.969256f,  1.875991f,  0.041556f };
			const float bWeight[3] = { 0.055648f, -0.204043f,  1.057311f };

			//rgb[3*offset  ] = rWeight[0]*xyz[0] + rWeight[1]*xyz[1] + rWeight[2]*xyz[2];
			//rgb[3*offset+1] = gWeight[0]*xyz[0] + gWeight[1]*xyz[1] + gWeight[2]*xyz[2];
			//rgb[3*offset+2] = bWeight[0]*xyz[0] + bWeight[1]*xyz[1] + bWeight[2]*xyz[2];
			alpha[offset] = 1.0f;//(*pixels)(x, y).alpha;

			rgb[3*offset  ] = rWeight[0]*xyz[0] + rWeight[1]*xyz[1] + rWeight[2]*xyz[2];
			rgb[3*offset+1] = gWeight[0]*xyz[0] + gWeight[1]*xyz[1] + gWeight[2]*xyz[2];
			rgb[3*offset+2] = bWeight[0]*xyz[0] + bWeight[1]*xyz[1] + bWeight[2]*xyz[2];

			// Normalize pixel with weight sum
			//float weightSum = (*pixels)(x, y).weightSum;
			float weightSum = pixels[x][y].weightSum;
			if (weightSum != 0.0f) 
			{
				float invWt = 1.0f / weightSum;
				//invWt = 1.0f;// / MaxWeight;
				//rgb[3*offset  ] = Clamp(rgb[3*offset  ] * invWt, 0.0f, INFINITY);
				//rgb[3*offset+1] = Clamp(rgb[3*offset+1] * invWt, 0.0f, INFINITY);
				//rgb[3*offset+2] = Clamp(rgb[3*offset+2] * invWt, 0.0f, INFINITY);

				rgb[3*offset  ] = Clamp(rgb[3*offset  ] * invWt, 0.0f, INFINITY);
				rgb[3*offset+1] = Clamp(rgb[3*offset+1] * invWt, 0.0f, INFINITY);
				rgb[3*offset+2] = Clamp(rgb[3*offset+2] * invWt, 0.0f, INFINITY);

				//rgb[3*offset  ] = Clamp(weightSum, 0.0f, INFINITY);
				//rgb[3*offset+1] = Clamp(weightSum, 0.0f, INFINITY);
				//rgb[3*offset+2] = Clamp(weightSum, 0.0f, INFINITY);

				alpha[offset] = 1.0f;//Clamp(alpha[offset] * invWt, 0.0f, 1.0f);
			}
			// Compute premultiplied alpha color
			if (premultiplyAlpha) 
			{
				rgb[3*offset  ] *= alpha[offset];
				rgb[3*offset+1] *= alpha[offset];
				rgb[3*offset+2] *= alpha[offset];
			}

			mitsuba::Point2i sample(x + xPixelStart, y + yPixelStart);
			mitsuba::Spectrum spec;
			spec.fromLinearRGB(rgb[3*offset], rgb[3*offset+1], rgb[3*offset+2]);

			float alpha2 = alpha[offset];

			//printf("Sample %s  |  %s\n", sample.toString().c_str(), spec.toString().c_str());

			imageblock->getBitmap()->setPixel(sample, spec);
			//imageblock->put(sample, spec, alpha2);

			++offset;
		}
	}

	// Write RGBA image
	//WriteRGBAImage(filename+ ".exr", rgb, alpha, xPixelCount, yPixelCount, xResolution, yResolution, xPixelStart, yPixelStart);

	// Release temporary image memory
	delete[] alpha;
	delete[] rgb;
}

SamoaImageFilm createSamoaFilm(int w, int h, float crop[4]) {
	//float crop[4] = { 0, 1, 0, 1 };

	// Defaults to true
	const bool premultiplyAlpha = true;

	return SamoaImageFilm(w, h, crop, premultiplyAlpha);
}

