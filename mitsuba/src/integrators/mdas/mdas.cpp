/*
 * mdas.cpp
 *
 *  Created on: Sep 2, 2019
 *      Author: miguel
 */


#include <mitsuba/render/film.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include "mdas.h"

#include "../quad/pathTracer.h"
#include "../quad/pathMIS.h"

#include "samoa2.h"
#include "tmp_sampler.h"

MTS_NAMESPACE_BEGIN

class MdasIntegrator : public Integrator {
public:
	MdasIntegrator(const Properties& props): Integrator(props) {
		//  Save properties

		m_hideEmitters = props.getBoolean("hide_emitters", false);

		num_blocks = static_cast<float>(props.getInteger("num_blocks", 1));

		ps = props.getFloat("spp", 2);
		max_depth = props.getInteger("maxDepth", 1);
		area_light = props.getBoolean("area_light", false);
		n_dims = props.getInteger("n_dims", 2);
		m_typeSubIntegrator = props.getString("typeSubIntegrator", "path");

		reconfile = props.getString("reconfile", "recon.exr");
		densityfile = props.getString("densityfile", "density.exr");

		m_sampleDensity = props.getBoolean("sampleDensity", false);
		m_solid_angle = props.getBoolean("solid_angle_sampling", false);

		candidatenum = props.getInteger("numcandidate", 4);
		knnsnum = props.getInteger("numknns", 4);
		bucketnum = props.getInteger("numcellsmp", 4);
		cellknnsnum = props.getInteger("numreconknns", 15);
		initnum = props.getInteger("numinit", 1024);
		imagescale = props.getFloat("scaleimageaxes", 1.0f);
		luminancelimit = props.getFloat("limitluminance", 100.0f);
		luminancescale = props.getFloat("scaleluminance", 40.0f);
		distanceepsilon = props.getFloat("epsilondistance", 1E-5f);
		maxgs = props.getFloat("maxgaussianscale", 10.0f);
		optstep = props.getInteger("optimizationstep", 10);
	}

	MdasIntegrator(Stream* stream, InstanceManager* manager) : Integrator(stream, manager) {

	}

	virtual bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
	                    int sceneResID, int sensorResID, int samplerResID) {

		return true;


		// Set here options -> no need to enter in mitsuba object
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
	                int sceneResID, int sensorResID, int samplerResID) {

		m_scene = scene;
		m_film = scene->getFilm();

		if (m_typeSubIntegrator == "path") {
			// Create integrator
			ref<PathTracerQuadrature> tmp_integrator = static_cast<PathTracerQuadrature *> (PluginManager::getInstance()->
											createObject(MTS_CLASS(MonteCarloIntegrator), Properties("pathQuadrature")));

			if (area_light && m_solid_angle) {
				tmp_integrator->m_solid_angle_sampling = true;
			}
			tmp_integrator->m_directSampling = (max_depth == 1) ? true : false;

			m_integrator = tmp_integrator;
		} else if (m_typeSubIntegrator == "pathMIS") {
			ref<PathTracerMIS> tempPtr = static_cast<PathTracerMIS *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(MonteCarloIntegrator), Properties("pathQuadratureMIS")));
			
			m_integrator = tempPtr;
		} else {
			Log(EError, "Type of subIntegrator undefined");
		}

		// Only direct illumination
		m_integrator->m_maxDepth = max_depth;
		m_integrator->m_hideEmitters = m_hideEmitters;
		

		// m_integrator = static_cast<MonteCarloIntegrator *> (PluginManager::getInstance()->
		// 				                createObject(MTS_CLASS(MonteCarloIntegrator), Properties("path")));
		// // Only direct illumination
		// m_integrator->m_maxDepth = (max_depth == 1) ? max_depth + 1 : max_depth;
		// m_integrator->m_hideEmitters = m_hideEmitters;

		if (n_dims == 2) {
			renderInternal<2>();
		} else if (n_dims == 3) {
			renderInternal<3>();
		} else if (n_dims == 4) {
			renderInternal<4>();
		} else if (n_dims == 5) {
			renderInternal<5>();
		} else if (n_dims == 6) {
			renderInternal<6>();
		} else {
			Log(EError, "Max number of dimensions is 6. Modify code to increase it on demand");
		}

		/*
		 * Reconstruct image
		 */
		//ref<mitsuba::ImageBlock> image = new ImageBlock(Bitmap::EPixelFormat::ESpectrum,film->getSize(), nullptr);

		return true;


		// Initialize "sampler of mdas"
		// Rewrite film to mitsuba style

		// Main loop
		// For every sample
		//		Get data from sampler
		//		Call integrator => direct lightning in first
		// 		Set last spectrum to data structure
		//		Do it meanwhile new samples
		// 		Finally reconstruct the image
		//		Optionally implement density functions -> compare with ours
	}

	template <int DIMENSIONS>
	void renderInternal() {
		if (area_light) {
			Log(EInfo, "Rendering %d dimensions with area light sampling [Depth : %d]", DIMENSIONS, max_depth);
			renderInternal2<DIMENSIONS, true>();
		} else {
			Log(EInfo, "Rendering %d dimensions without area light sampling [Depth : %d]", DIMENSIONS, max_depth);
			renderInternal2<DIMENSIONS, false>();
		}
	}

	template <int DIMENSIONS, bool AREA_LIGHT>
	void renderInternal2() {
		Sensor *sensor = m_scene->getSensor();
		Film* film = sensor->getFilm();

		int width = m_film->getSize().x;
		int height = m_film->getSize().y;

		int xstart = 0;
		int xend;
		int ystart = 0;
		int yend;

		m_film->clear();

		// Need to divide image in blocks because of memory
		// initnum / blocks ?

		int X_INC = width / num_blocks;
		int Y_INC = height / num_blocks;

		float crop[4];

		for (int xB = 1; xB <= num_blocks; xB++){

			if (xB == num_blocks) xend = width;
			else xend = xB * X_INC;

			for(int yB = 1; yB <= num_blocks; yB++){

				if (yB == num_blocks) yend = height;
				else yend = yB * Y_INC;

				// Process area
				Samoa2Sampler<DIMENSIONS, AREA_LIGHT> samoa2sampler (candidatenum, knnsnum, bucketnum, cellknnsnum, maxgs, optstep, initnum, imagescale,
								             luminancelimit, luminancescale, distanceepsilon, xstart, xend, ystart, yend, ps, width,
											 height, reconfile, densityfile);


				crop[0] = (xB-1) / num_blocks;
				crop[1] = (xB) / num_blocks;
				crop[2] = (yB-1) / num_blocks;
				crop[3] = (yB) / num_blocks;

				Log(EInfo, "Range x: %d %d  y: %d %d  crop: %f %f %f %f", xstart, xend, ystart, yend, crop[0], crop[1], crop[2], crop[3]);
				Log(EInfo, "%d %d", xB, yB);

				samoa2sampler.samoaFilm = createSamoaFilm(width, height, crop);
				//samoa2sampler.imageblock = new ImageBlock(Bitmap::EPixelFormat::ESpectrum, Vector2i(xend-xstart, yend-ystart), nullptr);
				//samoa2sampler.imageblock->setOffset(Point2i(xstart, ystart));

				Sample sample;
				bool stillWork;

				//for (int ID_sample = 0; ID_sample < samoa2sampler.totalSamples; ID_sample++) {
				while(true) {
					stillWork = samoa2sampler.GetNextSample(&sample);
					if (!stillWork) break;

					TMP_SAMPLER tmp_sampler(sample.twoD[0]);
					RadianceQueryRecord rRec(m_scene, &tmp_sampler);
					rRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());

					Point2 sample_point(sample.imageX, sample.imageY);
					Point2 apertureSample(0.5f);
					float timeSample = 0.5f;

					if (sensor->needsApertureSample()) {
					   //apertureSample = m_sampler->next2D();
					   apertureSample = Point2(sample.lensU, sample.lensV);
					}

					if (sensor->needsTimeSample()) {
					   //timeSample = m_sampler->next1D();
					   timeSample = sample.time;
					}

					//Log(EDebug, "Obtenidas muestras de tiempo");
					RayDifferential eyeRay;
					eyeRay.scaleDifferential(0);
					Spectrum sampleValue = sensor->sampleRayDifferential(
										eyeRay, sample_point, apertureSample, timeSample);

					sampleValue *= m_integrator->Li(eyeRay, rRec);

					// Save back result
					mitsuba::Spectrum Ls = samoa2sampler.SetPrevSample(sampleValue);
				}
				/*ref<ImageBlock> imageBlock = new ImageBlock(Bitmap::EPixelFormat::ESpectrum, Vector2i(xend-xstart, yend-ystart), nullptr);
				imageBlock->setOffset(Point2i(xstart, ystart));
				imageBlock->setSize(Vector2i(xend-xstart, yend-ystart));

				samoa2sampler.samoaFilm.WriteImage(imageBlock->getBitmap()->getFloatData());
				//film->put(samoa2sampler.imageblock);
				film->addBitmap(imageBlock->getBitmap(), 1.0f);*/

				ref<ImageBlock> imageBlock = new ImageBlock(Bitmap::ESpectrum, m_film->getSize());
				samoa2sampler.samoaFilm.WriteImage(imageBlock.get());

				//film->put(imageBlock);
				m_film->addBitmap(imageBlock->getBitmap());

				if (m_sampleDensity) {
					ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, m_film->getSize(), nullptr);
					result->clear();

					samoa2sampler.generateSampleDistribution(result);
					generateFilmData(result, "density.exr");
				}

				ystart = yend;
			}
			xstart = xend;
			ystart = 0;
		}
		
	}

	void generateFilmData(ImageBlock* imageBlock, std::string filename) {
		Properties props("hdrfilm");
		props.setBoolean("banner", false);
		props.setInteger("width", m_film->getSize().x);
		props.setInteger("height", m_film->getSize().y);
		ref<Film> frame = static_cast<Film *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), props));
		
		fs::path frameDestination = m_scene->getDestinationFile().parent_path();
		frameDestination /= filename;
		frame->setDestinationFile(frameDestination, 0);

		frame->setBitmap(imageBlock->getBitmap());
		frame->develop(m_scene, 0);
	}

	void cancel() {

	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MdasIntegrator[" << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	bool m_hideEmitters;
	float num_blocks;
	int max_depth;
	mitsuba::Float ps;
	bool area_light;
	int n_dims;

	int optstep;
	std::string densityfile;
	std::string reconfile;
	int candidatenum;
	int knnsnum;
	int bucketnum;
	int cellknnsnum;
	int initnum;
	mitsuba::Float imagescale;
	mitsuba::Float luminancelimit;
	mitsuba::Float distanceepsilon;
	mitsuba::Float luminancescale;
	mitsuba::Float maxgs;

	Scene *m_scene;
	Film *m_film;
	ref<MonteCarloIntegrator> m_integrator;
	//ref<MonteCarloIntegrator> m_integrator;

	bool m_sampleDensity;
	bool m_solid_angle;
	std::string m_typeSubIntegrator;

};

MTS_IMPLEMENT_CLASS_S(MdasIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(MdasIntegrator, "Multidimensional Hachisuka");

MTS_NAMESPACE_END


