
#include <string>

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>

#include "tupleMath.h"
#include "quadSampler.h"
#include "pathTracer.h"
#include "singleScattering.h"
#include "singleScatteringEquiangular.h"
#include "singleScatteringVRL.h"


#include "../../../ext/viltrum/quadrature/integrate-bins.h"
#include "../../../ext/viltrum/quadrature/nested.h"
#include "../../../ext/viltrum/quadrature/monte-carlo.h"
#include "../../../ext/viltrum/quadrature/integrate-bins-adaptive.h"
#include "../../../ext/viltrum/quadrature/control-variates.h"
#include "../../../ext/viltrum/quadrature/integrate-adaptive-control-variates.h"
#include "../../../ext/viltrum/quadrature/integrate-bins-adaptive-precalculate.h"
#include "../../../ext/viltrum/quadrature/integrate-optimized-adaptive-stratified-control-variates.h"
#include "../../../ext/viltrum/quadrature/munoz2014.h"

using namespace viltrum;

#include <chrono>
using milliseconds = std::chrono::milliseconds;
using microseconds = std::chrono::microseconds;

MTS_NAMESPACE_BEGIN

// SFINAE vodoo to check C is iterable
template<typename C>
struct is_iterable
{
  typedef long false_type;
  typedef char true_type;

  template<class T> static false_type check(...);
  template<class T> static true_type  check(int,
                    typename T::const_iterator = C().end());

  enum { value = sizeof(check<C>(0)) == sizeof(true_type) };
};

template<class T, class R = void>
struct enable_if_type { typedef R type; };

template<class T, class Enable = void>
struct is_integrator : std::false_type {};

template<class T>
struct is_integrator<T, typename enable_if_type<typename T::is_integrator_tag>::type> : std::true_type
{};

// Compile time check of the type of VILTRUM's machinery
#define ENABLE_IF(...) \
  typename std::enable_if<__VA_ARGS__::value>::type* = nullptr

/*
 * Interface that connects VILTRUM with Mitsuba
 * It generates a function F that calls Mitsuba's scene to generate the value of the render integral
 * in the points needed by our method
 */
template<int NDIM, typename F>
class FunctionWrapper {
	F f;
	mutable std::shared_ptr<unsigned long> nsamples;
	int sppQuadSamples;
public:
	FunctionWrapper(const F& f, int _sppQuadSamples) : f(f), sppQuadSamples(_sppQuadSamples), nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper(F&& f, int _sppQuadSamples) : f(std::forward<F>(f)), sppQuadSamples(_sppQuadSamples), nsamples(std::make_shared<unsigned long>(0)) {}

	template <typename RNG>
	TupleMath operator()(const std::array<double, NDIM>& pos,
						 RNG& rng,
						 bool customRNG) const {
		TupleMath value = f(pos, rng, customRNG);
		(*nsamples) += 1;
		return value;
	}

	TupleMath operator()(const std::array<double, NDIM>& pos) const {
		(*nsamples) += 1;

		std::default_random_engine rng (0);
		if (sppQuadSamples != 0) {
			TupleMath value(0.0);
			for (int i=0; i<sppQuadSamples; i++) {
				value += f(pos, rng, true);
			}
			return value / static_cast<double>(sppQuadSamples);
		} else {
			TupleMath value = f(pos, rng, false);
			return value;
		}
	}
	unsigned long samples() const { return (*nsamples); }
	void setSamples(unsigned long new_samples) { *nsamples = new_samples; }
};


class FunctionResidual {
public:
	mutable ref<MonteCarloIntegrator> ref_integrator;
	int depth;
	int depthQuadrature;

	template<typename F, typename Float, std::size_t DIM, typename RNG>
	auto sample(const F& f, const Range<Float,DIM>& range, RNG& rng) const {
		std::array<Float,DIM> sample;
	    for (std::size_t i=0;i<DIM;++i) {
		    std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
		    sample[i] = dis(rng);
		}

		if (depth != depthQuadrature) {
			ref_integrator->m_maxDepth = depth;
		}

		return std::make_tuple(f(sample, rng, true),sample);
	}

	template<typename F, typename Float, std::size_t DIM, typename RNG>
	auto sampleLow(const F& f, const Range<Float,DIM>& range, RNG& rng) const {
		std::array<Float,DIM> sample;
	    for (std::size_t i=0;i<DIM;++i) {
		    std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
		    sample[i] = dis(rng);
		}

		ref_integrator->m_maxDepth = depthQuadrature;
		return std::make_tuple(f(sample, rng, true),sample);
	}

	template<typename F, typename Float, std::size_t DIM, typename RNG>
	auto sampleHigher(const F& f, const Range<Float,DIM>& range, RNG& rng) const {
		std::array<Float,DIM> sample;
	    for (std::size_t i=0;i<DIM;++i) {
		    std::uniform_real_distribution<Float> dis(range.min(i),range.max(i));
		    sample[i] = dis(rng);
		}

		if (depth != depthQuadrature) {
			// Interface to not loose augmenting depth and using it in the residual
			static_cast<MonteCarloIntegrator *>(ref_integrator.get())->maxDepthLower = depthQuadrature;
			ref_integrator->m_maxDepth = depth;
		}

		auto valueFunction = f(sample, rng, true);
		auto valueHigherFunction = static_cast<MonteCarloIntegrator *>(ref_integrator.get())->LiHigher;
		return std::make_tuple(valueFunction,sample,TupleMath(valueHigherFunction));
	}

	FunctionResidual(ref<MonteCarloIntegrator> _ref, int _depth, int _depthQuadrature): ref_integrator(_ref), depth(_depth), depthQuadrature(_depthQuadrature) {}
};


StatsCounter statistic_totalSamples("Primary-Space Adaptive CV", "Total number of samples", EStatsType::ENumberValue);

/*
 * Implements Primary-Space Adaptive Control Variates as a plugin in Mitsuba
 */
class PSACV : public Integrator {
public:
	enum DATA_IMAGE {
		ERROR,
		DIVISIONS,
		QUAD
	};

	PSACV(const Properties& props): Integrator(props) {
        m_spp_cv = props.getFloat("spp_cv", -1.0);
        m_spp = props.getFloat("spp", 1.0);
		m_numberDimensions = props.getInteger("n_dims", 2);
		m_maxIterations = props.getInteger("maxIterations", 1);
		m_useMonteCarlo = props.getBoolean("monteCarlo", false);
		m_maxDepth = props.getInteger("maxDepth", 1);
		m_typeIntegrator = props.getString("typeIntegrator", "cv_stratified");
		m_typeSubIntegrator = props.getString("typeSubIntegrator", "path");
		m_higherQuadRule = props.getString("higherQuadRule", "simpson");
		m_hideEmitters = props.getBoolean("hideEmitters", false);
		m_error_size_weight = props.getFloat("error_size_weight", 0.0001);
		m_solid_angle_sampling = props.getBoolean("solid_angle_sampling", false);
		m_scaleSpp = props.getFloat("scaleSpp", -1.0);
		m_strictNormals = props.getBoolean("strictNormals", false);
		m_renderQuad = props.getBoolean("renderQuad", false);
		m_maxDepthQuad = props.getInteger("maxDepthQuad", -1);
		m_modeMIS = props.getString("modeMIS", "all");


		m_error_image = props.getBoolean("error_image", false);
		m_divisions_image = props.getBoolean("divisions_image", false);
		m_path_images = props.getString("path_images", "");
		m_isVideo = props.getBoolean("isVideo", false);
		m_FPS = props.getInteger("fps", 30);

		m_sppQuadSamples = props.getInteger("sppQuadSamples", 0);

		//Single scattering props
        m_mapping_0_ts = props.getBoolean("mappingTs", false);

		//VRL props
		m_raySpectrum = props.getSpectrum("vrlSpectrum", Spectrum(1.0));
		m_rayDir = props.getVector("vrlDir", Vector(1.0,0.0,0.0));
		m_rayOrigin = props.getPoint("vrlOrigin", Point(0.0,0.0,0.0));
		m_vrlDistance = props.getFloat("vrlDistance", 10.0);
		m_nb_samples_vrl = props.getInteger("samples_vrl", 1);
	}

    PSACV(Stream* stream, InstanceManager* manager) : Integrator(stream, manager) {}

    void serialize(Stream* stream, InstanceManager* manager) const {
        Integrator::serialize(stream, manager);
    }

    virtual bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
                    int sceneResID, int sensorResID, int samplerResID) {
        Integrator::preprocess(scene, queue, job, sceneResID,
                               sensorResID, samplerResID);

        return true;
    }

    const Integrator *getSubIntegrator(int idx) const {
    	if (idx != 0){
    		return nullptr;
    	}
    	return m_subIntegrator.get();
    }

    void addChild(const std::string &name, ConfigurableObject *child) {}

    void cancel() {
        m_subIntegrator->cancel();
    }

	auto errorEstimator() {
		/* Standard absolute error */
		return error_single_dimension_size(m_error_size_weight);
	}

    /*
     * Configure all
     */
    bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
                int sceneResID, int sensorResID, int samplerResID) {

    	ref<Scheduler> scheduler = Scheduler::getInstance();

		m_scene = scene;
    	m_sensor = scene->getSensor();
    	m_sampler = m_sensor->getSampler();
		m_film = m_sensor->getFilm();
		m_job = job;
		m_useCustomSampler = false;

		/* Compute real size if cropped film */
		realSizeImage = m_film->getCropSize();

		realMinRangeImage.x = static_cast<double>(m_film->getCropOffset().x) / static_cast<double>(m_film->getSize().x);
		realMinRangeImage.y = static_cast<double>(m_film->getCropOffset().y) / static_cast<double>(m_film->getSize().y);

		realMaxRangeImage.x = (static_cast<double>(m_film->getCropOffset().x) + static_cast<double>(m_film->getCropSize().x)) / (static_cast<double>(m_film->getSize().x));
		realMaxRangeImage.y = (static_cast<double>(m_film->getCropOffset().y) + static_cast<double>(m_film->getCropSize().y)) / (static_cast<double>(m_film->getSize().y));

		// Get scale of samples given
		if (m_scaleSpp != -1.0) {
			m_spp_cv = m_spp * (1.0 / m_scaleSpp);
		}

        // Calculate different
        if (m_spp_cv != -1.0) {
			m_spp = m_spp - m_spp_cv;

			if (m_spp <= 0) {
				Log(EInfo, "Initialization spp for CV was negative or zero : %f", m_spp);
			}
        }

		Log(EInfo, "SPP %f   SPP_CV %f  SCALE %f", m_spp, m_spp_cv, 1.0/m_scaleSpp);

		/*
		 *  Create subintegrator
		 */

		if (m_typeSubIntegrator == "path") {
			ref<PathTracerQuadrature> tempPtr = static_cast<PathTracerQuadrature *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(MonteCarloIntegrator), Properties("pathQuadrature")));

			if (m_maxDepthQuad == -1) {
				tempPtr->m_directSampling = (m_maxDepth == 1) ? true : false;
			} else {
				tempPtr->m_directSampling = (m_maxDepthQuad == 1) ? true : false;
			}
			if (m_solid_angle_sampling) {
				tempPtr->m_solid_angle_sampling = m_solid_angle_sampling;
				tempPtr->m_directSampling = true;
			}

			m_subIntegrator = tempPtr;
		} else if (m_typeSubIntegrator == "singleScattering") {
			ref<SingleScatteringQuadrature> tempPtr = static_cast<SingleScatteringQuadrature *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(MonteCarloIntegrator), Properties("singleScatteringQuadrature")));

			tempPtr->m_directSampling = (m_maxDepth == 1) ? true : false;
			tempPtr->m_mapping_0_ts = m_mapping_0_ts;
			if (m_solid_angle_sampling) {
				tempPtr->m_solid_angle_sampling = m_solid_angle_sampling;
				tempPtr->m_directSampling = true;
			}

			m_subIntegrator = tempPtr;
		} else if (m_typeSubIntegrator == "singleScatteringEquiangular") {
			ref<SingleScatteringQuadratureEquiangular> tempPtr = static_cast<SingleScatteringQuadratureEquiangular *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(MonteCarloIntegrator), Properties("singleScatteringQuadratureEquiangular")));

			tempPtr->m_directSampling = (m_maxDepth == 1) ? true : false;
			tempPtr->m_mapping_0_ts = m_mapping_0_ts;
			if (m_solid_angle_sampling) {
				tempPtr->m_solid_angle_sampling = m_solid_angle_sampling;
				tempPtr->m_directSampling = true;
			}

			m_subIntegrator = tempPtr;
		}  else if (m_typeSubIntegrator == "singleScatteringVRL") {
			ref<SingleScatteringQuadratureVRL> tempPtr = static_cast<SingleScatteringQuadratureVRL *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(MonteCarloIntegrator), Properties("singleScatteringQuadratureVRL")));

			tempPtr->m_directSampling = (m_maxDepth == 1) ? true : false;
			tempPtr->m_raySpectrum = m_raySpectrum;
			tempPtr->m_rayDir = m_rayDir;
			tempPtr->m_rayOrigin = m_rayOrigin;
			tempPtr->m_vrlDistance = m_vrlDistance;
			tempPtr->m_nb_samples_vrl = m_nb_samples_vrl;

			m_subIntegrator = tempPtr;
		} else if (m_typeSubIntegrator == "ao") {
			// Create ambient oclusion subintegrator
			Properties props("ao");
			props.setSize("shadingSamples", 1);
			ref<SamplingIntegrator> tempPtr = static_cast<SamplingIntegrator *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(SamplingIntegrator), props));
			tempPtr->preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
			m_subIntegrator = tempPtr;
		} else {
			Log(EError, "Type of subIntegrator undefined. Check documentation!");
		}

		MonteCarloIntegrator* subIntegrator = dynamic_cast<MonteCarloIntegrator*>(m_subIntegrator.get());

		if (subIntegrator) {

			if (m_maxDepthQuad == -1) {
				m_maxDepthQuad = m_maxDepth;
			}

			subIntegrator->m_maxDepth = (m_maxDepthQuad != -1) ? m_maxDepthQuad : m_maxDepth;

			if (m_sppQuadSamples != 0) {
				subIntegrator->m_maxDepth = m_maxDepth;
			}

			subIntegrator->m_rrDepth = 5000;
			subIntegrator->m_hideEmitters = m_hideEmitters;
			subIntegrator->m_strictNormals = m_strictNormals;

			if (m_modeMIS == "all") {
				subIntegrator->m_modeMIS = MonteCarloIntegrator::MIS::ALL;
			} else if (m_modeMIS == "emitter") {
				subIntegrator->m_modeMIS = MonteCarloIntegrator::MIS::EMITTER;
			} else if (m_modeMIS == "brdf") {
				subIntegrator->m_modeMIS = MonteCarloIntegrator::MIS::BRDF;
			} else if (m_modeMIS == "emitter_weight") {
				subIntegrator->m_modeMIS =  MonteCarloIntegrator::EMITTER_WEIGHT;
			} else if (m_modeMIS == "brdf_weight") {
				subIntegrator->m_modeMIS =  MonteCarloIntegrator::MIS::BRDF_WEIGHT;
			} else {
				Log(EError, "[modeMIS] Type of MIS undefined  | Available : all , emitter , brdf , emitter_weight , brdf_weight");
			}
		}

		renderStart();

		timeSampling *= 0.001;
		Log(EInfo, "[TIME *] Sampling routine callback %ld milliseconds", timeSampling);

		// If choosed, render also the only quad result
		if (m_renderQuad) {
			Log(EInfo, "Rendering adicional quadrature ONLY result");
			m_typeIntegrator = "quad";
			// Override spp used because in quad, stepping use it
			m_spp = m_spp_cv;
			m_spp_cv = 0;
			renderStart("quad");
		}

        return true;
    }

	void renderStart(std::string pathResult="") {

		int N;
		if (m_higherQuadRule == "boole") {
			N = Boole::samples;
		} else if (m_higherQuadRule == "simpson") {
			N = Simpson::samples;
		}

		const Vector2i size_image = realSizeImage;
		const int D = m_numberDimensions;
		const int samples_per_iteration = (N-1)*std::pow(N, D-1);
		m_total_bins = size_image.x * size_image.y;
		if (m_isVideo) m_total_bins *= m_FPS * m_sensor->getShutterOpenTime();

		m_maxIterations = (m_spp_cv * m_total_bins) / samples_per_iteration;
		m_maxIterationsQuad = (m_spp * m_total_bins) / samples_per_iteration;

		Log(EInfo, "Using %d iterations from %f spp requested for CV init", m_maxIterations, m_spp_cv);

		if (m_higherQuadRule == "boole") {
			Log(EInfo, "Using boole-simpsons rules");
			m_samples_rule = Boole::samples;
			renderInternal(nested(boole, simpson), pathResult);
		} else if (m_higherQuadRule == "simpson") {
			Log(EInfo, "Using simpsons-trapezoidal rules");
			m_samples_rule = Simpson::samples;
			renderInternal(nested(simpson, trapezoidal), pathResult);
		} else {
			Log(EError, "Higher quadrature rule undefined");
		}
	}

	template <typename NESTEDQ>
	void renderInternal(NESTEDQ&& nestedQ, std::string pathResult="") {

		// Functor for calculating the residual + higher order radiance
		ref<MonteCarloIntegrator> refTmp = static_cast<MonteCarloIntegrator*>(m_subIntegrator.get());
		FunctionResidual customSampler(refTmp, m_maxDepth, m_maxDepthQuad);

		if (m_typeIntegrator == "quad") {
			auto stepper = stepper_bins_adaptive(
				std::forward<NESTEDQ>(nestedQ),
				errorEstimator());
			renderInternalDim(stepper, pathResult);
		} else if (m_typeIntegrator == "mc") {
			Log(EInfo, "Using stepper MC [%f spp]", m_spp);
			size_t seed = std::random_device()();
			auto stepper = stepper_bins_per_bin(stepper_monte_carlo_uniform(seed));
			renderInternalDim(stepper, pathResult);
		} else if (m_typeIntegrator == "cv_pixel_alphaOpt") {
			auto seed = std::random_device()();
			std::mt19937_64 rng (seed);
			auto rg = region_generator(std::forward<NESTEDQ>(nestedQ), errorEstimator(), m_maxIterations);

			auto integrator = integrator_stratified_pixel_control_variates(std::forward<decltype(rg)>(rg),
			AlphaOptimized(),
			std::forward<decltype(customSampler)>(customSampler),
			std::forward<decltype(rng)>(rng), m_spp);
			renderInternalDim(integrator, pathResult);
		} else if (m_typeIntegrator == "munoz2014") {
			const unsigned long munoz_spp_pixel = 8;
			auto integrator = integrator_bins_munoz_2014(m_spp, munoz_spp_pixel, m_error);
			renderInternalDim(integrator, pathResult);
		} else {
			Log(EError, "Type of integrator undefined");
		}
	}

	template <typename STEPPER>
	void renderInternalDim(STEPPER& stepper, std::string pathResult="") {
		if (m_numberDimensions == 2) {
			m_nDim = 2;
			auxRenderInternal<2>(stepper, pathResult);
		} else if (m_numberDimensions == 3) {
			m_nDim = 3;
			auxRenderInternal<3>(stepper, pathResult);
		}
		else if (m_numberDimensions == 4) {
			m_nDim = 4;
			auxRenderInternal<4>(stepper, pathResult);
		}
		else if (m_numberDimensions == 5) {
			m_nDim = 5;
			auxRenderInternal<5>(stepper, pathResult);
		}
		else if (m_numberDimensions == 6) {
			m_nDim = 6;
			auxRenderInternal<6>(stepper, pathResult);
		}
		else if (m_numberDimensions == 7) {
			m_nDim = 7;
			auxRenderInternal<7>(stepper, pathResult);
		}
		else if (m_numberDimensions == 8) {
			m_nDim = 8;
			auxRenderInternal<6>(stepper, pathResult);
		}
		else {
			Log(EError, "Max number of dimensions is 8. Modify code to increase it on demand!");
		}
	}

	template<int NDIMS, typename STEPPER>
	void auxRenderInternal(STEPPER& stepper, std::string pathResult="") {
		const Vector2i size_image = realSizeImage;

		if constexpr(NDIMS >= 3) {
			if (m_isVideo) {
				float shutterInterval = m_sensor->getShutterOpenTime();
				size_t totalFrames = m_FPS * shutterInterval;

				std::array<std::size_t, 3> resolutionBins = {std::size_t(size_image.x), std::size_t(size_image.y), totalFrames};
				auxRenderInternalFinal<NDIMS, true>(stepper, resolutionBins, pathResult);
			} else {
				std::array<std::size_t, 2> resolutionBins = {std::size_t(size_image.x), std::size_t(size_image.y)};
				auxRenderInternalFinal<NDIMS, false>(stepper, resolutionBins, pathResult);
			}
		} else {
			std::array<std::size_t, 2> resolutionBins = {std::size_t(size_image.x), std::size_t(size_image.y)};
			auxRenderInternalFinal<NDIMS, false>(stepper, resolutionBins, pathResult);
		}
	}

	template <int NDIMS, bool ISVIDEO, typename STEPPER, typename RESOLUTION,
	ENABLE_IF(is_integrator<STEPPER>)>
	void auxRenderInternalFinal(STEPPER& integrator, RESOLUTION& resolutionBins, std::string pathResult="") {

		auto f = [this] (const std::array<double, NDIMS>& positions, auto rng, bool customRNG) -> TupleMath { return this->callback_samplePoint(positions, rng, customRNG); };

		FunctionWrapper<NDIMS, decltype(f)> wrapper(f, m_sppQuadSamples);

		std::array<double, NDIMS> range_min, range_max;
		range_min.fill(0.0);
		range_max.fill(1.0);

		auto Frange = range(range_min, range_max);

		if constexpr(ISVIDEO) {
			float shutterInterval = m_sensor->getShutterOpenTime();
			size_t totalFrames = m_FPS * shutterInterval;

			Vector2i size_image = realSizeImage;
			std::vector<std::vector<std::vector<TupleMath> > > data_frames (totalFrames,std::vector<std::vector<TupleMath> >(size_image.x,std::vector <TupleMath>(size_image.y,0.0)));

			auto vb = [&data_frames] (const std::array<std::size_t,3>& i) -> TupleMath& { return data_frames[i[2]][i[0]][i[1]]; };

			Log(EInfo, "Starting integrating video [Frames %i] [Duration %f]", totalFrames, m_sensor->getShutterOpenTime());
			integrator.integrate(vb, resolutionBins, wrapper, Frange);
			Log(EInfo, "Done integrating video");

			int ID_frame = 0;
			for (std::size_t ff=0; ff<totalFrames; ff++) {
				Log(EInfo, "Saving frame %d", ID_frame);

				Properties props("ldrfilm");
				props.setBoolean("banner", false);
				props.setInteger("width", realSizeImage.x);
				props.setInteger("height", realSizeImage.y);
				ref<Film> frame = static_cast<Film *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(Film), props));

				fs::path frameDestination = m_scene->getDestinationFile().parent_path();
				frameDestination /= "frame_" + std::to_string(ID_frame) + ".png";
				frame->setDestinationFile(frameDestination, 0);

				ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

				for(int px=0; px<size_image.x; px++) {
					for(int py=0; py<size_image.y; py++) {
						result->getBitmap()->setPixel(Point2i(px,py), data_frames[ID_frame][px][py].spectrum());
					}
				}
				frame->setBitmap(result->getBitmap());
				frame->develop(m_scene, 0);

				ID_frame++;
			}
		} else {
			Vector2i size_image = realSizeImage;
			std::vector<std::vector<TupleMath>> data_pixels(size_image.x);
			for(auto& pc : data_pixels) {
				pc.resize(size_image.y, 0.0);
			}

			// Wrapper to final bins structure
			auto vb = [&data_pixels] (const std::array<std::size_t,2>& i) -> TupleMath& { return data_pixels[i[0]][i[1]]; };

			unsigned long samples_backup = wrapper.samples();
			wrapper.setSamples(0);

			Log(EInfo, "Launching");
			integrator.integrate(vb, resolutionBins, wrapper, Frange);
			Log(EInfo, "Done");

			Log(EInfo, "Total samples used :  Init(%d)  Iterations(%d)  Total(%f)", samples_backup, wrapper.samples(), ((m_spp_cv < 0) ? m_spp : m_spp_cv + m_spp) * m_total_bins);
			statistic_totalSamples += samples_backup + wrapper.samples();

			if (pathResult == "") {
				ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

				for(int px=0; px<size_image.x; px++) {
					for(int py=0; py<size_image.y; py++) {
						result->getBitmap()->setPixel(Point2i(px,py), data_pixels[px][py].spectrum());
					}
				}
				m_film->setBitmap(result->getBitmap());
			} else {
				renderImage(data_pixels, pathResult);
			}
		}
	}

	template <int NDIMS, bool ISVIDEO, typename STEPPER, typename RESOLUTION,
	ENABLE_IF(!is_integrator<STEPPER>)>
	void auxRenderInternalFinal(STEPPER& stepper, RESOLUTION& resolutionBins, std::string pathResult="") {
		auto f = [this] (const std::array<double, NDIMS>& positions, auto& rng, bool customRNG) -> TupleMath { return this->callback_samplePoint(positions, rng, customRNG); };

		FunctionWrapper<NDIMS, decltype(f)> wrapper(f, m_sppQuadSamples);

		std::array<double, NDIMS> range_min, range_max;
		range_min.fill(0.0);
		range_max.fill(1.0);

		auto Frange = range(range_min, range_max);

		Log(EInfo, "Init Stepper %f spp_cv [%d dimensions]", m_spp_cv, m_nDim);
		auto data_regions = stepper.init(resolutionBins, wrapper, Frange);

		unsigned long samples_backup = wrapper.samples();
		wrapper.setSamples(0);

		if (m_maxDepth != m_maxDepthQuad) {
			MonteCarloIntegrator* tmpRef = dynamic_cast<MonteCarloIntegrator*>(m_subIntegrator.get());
			if (tmpRef != nullptr) {
				tmpRef->m_maxDepth = m_maxDepth;
			}
			m_useCustomSampler = true;
		}

		unsigned int iterations = 0;
		unsigned int LIMIT_ITERATIONS = (m_typeIntegrator == "quad") ? m_maxIterationsQuad : m_spp;

		ref<ProgressReporter> progressReporter = new ProgressReporter("Processing QUAD stepper", LIMIT_ITERATIONS, m_job);
		Log(EInfo, "Iterating Stepper %d spp [%d dimensions]", LIMIT_ITERATIONS, m_nDim);
		while(iterations < LIMIT_ITERATIONS) {
			stepper.step(resolutionBins,wrapper,Frange,data_regions);
			iterations++;
			progressReporter->update(iterations);
		}
		progressReporter->finish();

		Log(EInfo, "Total samples used :  Init(%d)  Iterations(%d)  Total(%f)", samples_backup, wrapper.samples(), ((m_spp_cv < 0) ? m_spp : m_spp_cv + m_spp) * m_total_bins);
		statistic_totalSamples += samples_backup + wrapper.samples();
		Log(EInfo, "Total iterations used : %d", iterations);

		Log(EInfo, "Reconstructing Quad %f spp [%d dimensions]", m_spp, NDIMS);
		//progressReporter = new ProgressReporter("Reconstructing QUAD", data_regions.size(), job);

		if constexpr(ISVIDEO) {
			reconstructVideo(stepper, data_regions, wrapper, resolutionBins, Frange, pathResult);
		} else {
			reconstructImage(stepper, data_regions, wrapper, resolutionBins, Frange, pathResult);
		}
	}

	template <typename STEPPER, typename REGIONS, typename WRAPPER, typename RESOLUTION, typename FRANGE>
	void reconstructImage(STEPPER& stepper, REGIONS& data_regions, WRAPPER& wrapper, RESOLUTION& resolutionBins, FRANGE& Frange, std::string pathResult="") {
		const Vector2i size_image = realSizeImage;

		std::vector<std::vector<TupleMath>> data_pixels(size_image.x);
		for(auto& pc : data_pixels) {
			pc.resize(size_image.y, 0.0);
		}

		auto vb = [&data_pixels] (const std::array<std::size_t,2>& i) -> TupleMath& { return data_pixels[i[0]][i[1]]; };
		stepper.integral(vb, resolutionBins, wrapper, Frange, data_regions);

		if (pathResult == "") {
			ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

			for(int px=0; px<size_image.x; px++) {
				for(int py=0; py<size_image.y; py++) {
					result->getBitmap()->setPixel(Point2i(px,py), data_pixels[px][py].spectrum());
				}
			}
			m_film->setBitmap(result->getBitmap());
		} else {
			renderImage(data_pixels, pathResult);
		}

		// If requested generate images with error and divisions
		if (m_error_image) {
			generateImage(data_regions, wrapper, DATA_IMAGE::ERROR, "error.exr");
		}
	}

	template <typename STEPPER, typename REGIONS, typename WRAPPER, typename RESOLUTION, typename FRANGE>
	void reconstructVideo(STEPPER& stepper, REGIONS& data_regions, WRAPPER& wrapper, RESOLUTION& resolutionBins, FRANGE& Frange, std::string pathResult="frame") {
		float shutterInterval = m_sensor->getShutterOpenTime();
		size_t totalFrames = m_FPS * shutterInterval;

		ref<ProgressReporter> progressReporter = new ProgressReporter("Frame reconstruction", totalFrames, m_job);

		const Vector2i size_image = realSizeImage;
		std::vector<std::vector<std::vector<TupleMath> > > data_frames (totalFrames,std::vector<std::vector<TupleMath> >(size_image.x,std::vector <TupleMath>(size_image.y,0.0)));

		auto vb_temp = [&data_frames] (const std::array<std::size_t,3>& i) -> TupleMath& { return data_frames[i[2]][i[0]][i[1]]; };
		stepper.integral(vb_temp, resolutionBins, wrapper, Frange, data_regions);

		int ID_frame = 0;
		for (std::size_t ff=0; ff<totalFrames; ff++) {
			Log(EInfo, "Saving frame %d", ID_frame);

			Properties props("ldrfilm");
			props.setBoolean("banner", false);
			props.setInteger("width", realSizeImage.x);
			props.setInteger("height", realSizeImage.y);
			ref<Film> frame = static_cast<Film *> (PluginManager::getInstance()->
				createObject(MTS_CLASS(Film), props));

			fs::path frameDestination = m_scene->getDestinationFile().parent_path();
			frameDestination /= pathResult + "_" + std::to_string(ID_frame) + ".png";
			frame->setDestinationFile(frameDestination, 0);

			ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

			for(int px=0; px<size_image.x; px++) {
				for(int py=0; py<size_image.y; py++) {
					result->getBitmap()->setPixel(Point2i(px,py), data_frames[ID_frame][px][py].spectrum());
				}
			}
			frame->setBitmap(result->getBitmap());
			frame->develop(m_scene, 0);

			ID_frame++;
			progressReporter->update(ff);
		}
		progressReporter->finish();
	}

	template <typename REGIONS, typename WRAPPER,
			  typename std::enable_if<is_iterable<REGIONS>::value,
                                 REGIONS>::type* = nullptr>
	void generateImage(REGIONS& regions, WRAPPER& wrapper, PSACV::DATA_IMAGE type, std::string name) {
		Properties props("hdrfilm");
		props.setBoolean("banner", false);
		props.setInteger("width", realSizeImage.x);
		props.setInteger("height", realSizeImage.y);
		ref<Film> frame = static_cast<Film *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), props));

		fs::path frameDestination = m_scene->getDestinationFile().parent_path();
		std::string data_name = m_scene->getDestinationFile().filename().string();
		if (m_path_images != "") {
			frameDestination /= fs::path(m_path_images);
		} else {
			frameDestination /= data_name + "_data.exr";
		}
		frame->setDestinationFile(frameDestination, 0);
		Log(EInfo, "Saving image %s", frameDestination.c_str());

		std::vector<std::vector<TupleMath>> bins;
		bins.resize(realSizeImage.x);
		for (auto& v : bins) {
			v.resize(realSizeImage.y, TupleMath(0));
		}

		for (auto r : regions) {
			integrate_2D(bins, r, 0, 1, 0, 1, realSizeImage.x, realSizeImage.y, type);
		}

		ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

		for(int px=0; px<realSizeImage.x; px++) {
			for(int py=0; py<realSizeImage.y; py++) {
				result->getBitmap()->setPixel(Point2i(px,py), bins[px][py].spectrum());
			}
		}
		frame->setBitmap(result->getBitmap());
		frame->develop(m_scene, 0);

	}

	template <typename REGIONS, typename WRAPPER,
			  typename std::enable_if<!is_iterable<REGIONS>::value,
                                 REGIONS>::type* = nullptr>
	void generateImage(REGIONS& regions, WRAPPER& wrapper, PSACV::DATA_IMAGE type, std::string name) {
		Log(EWarn, "Nothing to save in image");
	}

	template <typename DATA_PIXELS>
	void renderImage(DATA_PIXELS& data, std::string path) {
		Properties props("hdrfilm");
		props.setBoolean("banner", false);
		props.setInteger("width", realSizeImage.x);
		props.setInteger("height", realSizeImage.y);
		ref<Film> frame = static_cast<Film *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Film), props));

		fs::path frameDestination = m_scene->getDestinationFile().parent_path();
		std::string data_name = m_scene->getDestinationFile().filename().string();

		frameDestination /= data_name + "_" + path + ".exr";
		frame->setDestinationFile(frameDestination, 0);
		Log(EInfo, "Saving image %s", frameDestination.c_str());

		ref<ImageBlock> result = new ImageBlock(Bitmap::ESpectrum, realSizeImage, nullptr);

		for(int px=0; px<realSizeImage.x; px++) {
			for(int py=0; py<realSizeImage.y; py++) {
				result->getBitmap()->setPixel(Point2i(px,py), data[px][py].spectrum());
			}
		}
		frame->setBitmap(result->getBitmap());
		frame->develop(m_scene, 0);
	}

	template <typename R>
	void integrate_2D(std::vector<std::vector<TupleMath>>& bins, const R& r, double xmin, double xmax, double ymin, double ymax, double sizeX, double sizeY, PSACV::DATA_IMAGE type) {
		double dx = (xmax - xmin)/sizeX;
		double local_xmin = std::get<0>(r.range())[0];
		double local_xmax = std::get<1>(r.range())[0];

		double dy = (ymax - ymin)/sizeY;
		double local_ymin = std::get<0>(r.range())[1];
		double local_ymax = std::get<1>(r.range())[1];

		std::size_t start_binX = std::max(std::size_t(0),std::size_t((local_xmin - xmin)/dx));
		std::size_t end_binX = std::min(std::size_t(sizeX)-1,std::size_t((local_xmax - xmin)/dx));

		std::size_t start_binY = std::max(std::size_t(0),std::size_t((local_ymin - ymin)/dy));
		std::size_t end_binY = std::min(std::size_t(sizeY)-1,std::size_t((local_ymax - ymin)/dy));

		auto errorFunctor = errorEstimator();

		for (std::size_t i = start_binX; i <= end_binX; ++i) {
			for (std::size_t j = start_binY; j <= end_binY; ++j) {
				if (i == start_binX || i == end_binX || j == start_binY || j == end_binY) {
					bins[i][j][0] = 1.0;		// Save division border
				}
				bins[i][j][2] += std::get<0>(errorFunctor(r));
			}
		}
	}

    void wakeup(ConfigurableObject *parent,
            std::map<std::string, SerializableObject *> &params) {
    }

    std::string toString() const {
    	std::ostringstream oss;
        oss << "PSACV[" << endl
            << "]";
        return oss.str();
    }

	template <typename Array, typename RNG>
	TupleMath callback_samplePoint (const Array& values, RNG& rng, bool customRNG) {
		auto t_sampling1 = std::chrono::high_resolution_clock::now();

		QuadratureSampler<values.size(), decltype(rng)> quadratureSampler(values, rng, customRNG);

		RadianceQueryRecord rRec(m_scene, &quadratureSampler);
		rRec.newQuery(RadianceQueryRecord::ESensorRay, m_sensor->getMedium());

		Point2 sample_pos(rRec.nextSample2D());
		Point2 sample_point;
		sample_point.x = sample_pos.x * realSizeImage.x;
		sample_point.y = sample_pos.y * realSizeImage.y;

		Point2 apertureSample(0.5f);
		float timeSample = 0.5f;

		if (m_sensor->needsTimeSample())
		   timeSample = rRec.nextSample1D();

		if (m_sensor->needsApertureSample())
		   apertureSample = rRec.nextSample2D();

		RayDifferential eyeRay;
		eyeRay.hasDifferentials = false;
		eyeRay.scaleDifferential(0);
		Spectrum sampleValue = m_sensor->sampleRay(
							eyeRay, sample_point, apertureSample, timeSample);

		sampleValue *= m_subIntegrator->Li(eyeRay, rRec);

		auto t_sampling2 = std::chrono::high_resolution_clock::now();
		auto d_sampling = std::chrono::duration_cast<microseconds>( t_sampling2 - t_sampling1 ).count();
		timeSampling += d_sampling;

		return TupleMath(sampleValue);
	}

    MTS_DECLARE_CLASS()
private:
    ref<Film> m_film;
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Sampler> m_sampler;
	ref<SamplingIntegrator> m_subIntegrator;

	int64_t timeSampling = 0;

	const RenderJob *m_job;

    // Error thresold
    float m_error;


	double m_spp_cv;
	unsigned int m_total_bins;
	unsigned long m_maxIterationsQuad;

	Vector2i realSizeImage;
	Vector2 realMinRangeImage;
	Vector2 realMaxRangeImage;

	double m_spp;
	double m_sppQuadSamples;
    int m_numberDimensions;
    unsigned long m_maxIterations;
    int m_maxDepth;
    bool m_useMonteCarlo;
	std::string m_typeIntegrator;
	std::string m_typeSubIntegrator;
	std::string m_higherQuadRule;
	bool m_hideEmitters;
	double m_error_size_weight;
	bool m_error_image;
	bool m_divisions_image;
	std::string m_path_images;
	bool m_solid_angle_sampling;

	int m_samples_rule;
	int m_nDim;
	bool m_isVideo;
	int m_FPS;

	double m_scaleSpp;
	bool m_strictNormals;
	bool m_renderQuad;
	int m_maxDepthQuad;
	bool m_useCustomSampler;
	std::string m_modeMIS;

	//Single scattering
	bool m_mapping_0_ts;

	//VRL
	Spectrum m_raySpectrum;
    Vector m_rayDir;
    Point m_rayOrigin;
    Float m_vrlDistance;
    int m_nb_samples_vrl;

};

MTS_IMPLEMENT_CLASS_S(PSACV, false, Integrator)
MTS_EXPORT_PLUGIN(PSACV, "Primary-Space Adaptive Control Variates using Piecewise Polynomials Approximations");

MTS_NAMESPACE_END
