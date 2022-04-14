
#ifndef SRC_INTEGRATORS_SINGLESCATTERINGVRL_H_
#define SRC_INTEGRATORS_SINGLESCATTERINGVRL_H_


/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength_singleScattering_vrl("Single scattering quadrature", "Average path length", EAverage);

/*! \plugin{path}{Path tracer}
 * \order{2}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *        which the implementation will start to use the ``russian roulette''
 *        path termination criterion. \default{\code{5}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See the description below
 *        for details.\default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 *
 * This integrator implements a basic path tracer and is a \emph{good default choice}
 * when there is no strong reason to prefer another method.
 *
 * To use the path tracer appropriately, it is instructive to know roughly how
 * it works: its main operation is to trace many light paths using \emph{random walks}
 * starting from the sensor. A single random walk is shown below, which entails
 * casting a ray associated with a pixel in the output image and searching for
 * the first visible intersection. A new direction is then chosen at the intersection,
 * and the ray-casting step repeats over and over again (until one of several
 * stopping criteria applies).
 * \begin{center}
 * \includegraphics[width=.7\textwidth]{images/integrator_path_figure.pdf}
 * \end{center}
 * At every intersection, the path tracer tries to create a connection to
 * the light source in an attempt to find a \emph{complete} path along which
 * light can flow from the emitter to the sensor. This of course only works
 * when there is no occluding object between the intersection and the emitter.
 *
 * This directly translates into a category of scenes where
 * a path tracer can be expected to produce reasonable results: this is the case
 * when the emitters are easily ``accessible'' by the contents of the scene. For instance,
 * an interior scene that is lit by an area light will be considerably harder
 * to render when this area light is inside a glass enclosure (which
 * effectively counts as an occluder).
 *
 * Like the \pluginref{direct} plugin, the path tracer internally relies on multiple importance
 * sampling to combine BSDF and emitter samples. The main difference in comparison
 * to the former plugin is that it considers light paths of arbitrary length to compute
 * both direct and indirect illumination.
 *
 * For good results, combine the path tracer with one of the
 * low-discrepancy sample generators (i.e. \pluginref{ldsampler},
 * \pluginref{halton}, or \pluginref{sobol}).
 *
 * \paragraph{Strict normals:}\label{sec:strictnormals}
 * Triangle meshes often rely on interpolated shading normals
 * to suppress the inherently faceted appearance of the underlying geometry. These
 * ``fake'' normals are not without problems, however. They can lead to paradoxical
 * situations where a light ray impinges on an object from a direction that is classified as ``outside''
 * according to the shading normal, and ``inside'' according to the true geometric normal.
 *
 * The \code{strictNormals}
 * parameter specifies the intended behavior when such cases arise. The default (\code{false}, i.e. ``carry on'')
 * gives precedence to information given by the shading normal and considers such light paths to be valid.
 * This can theoretically cause light ``leaks'' through boundaries, but it is not much of a problem in practice.
 *
 * When set to \code{true}, the path tracer detects inconsistencies and ignores these paths. When objects
 * are poorly tesselated, this latter option may cause them to lose a significant amount of the incident
 * radiation (or, in other words, they will look dark).
 *
 * The bidirectional integrators in Mitsuba (\pluginref{bdpt}, \pluginref{pssmlt}, \pluginref{mlt} ...)
 * implicitly have \code{strictNormals} set to \code{true}. Hence, another use of this parameter
 * is to match renderings created by these methods.
 *
 * \remarks{
 *    \item This integrator does not handle participating media
 *    \item This integrator has poor convergence properties when rendering
 *    caustics and similar effects. In this case, \pluginref{bdpt} or
 *    one of the photon mappers may be preferable.
 * }
 */
class SingleScatteringQuadratureVRL : public MonteCarloIntegrator {
public:
	SingleScatteringQuadratureVRL(const Properties &props)
        : MonteCarloIntegrator(props) {
            m_directSampling = props.getBoolean("directSampling", false);
            m_raySpectrum = props.getSpectrum("vrlSpectrum", Spectrum(1.0));
            m_rayDir = props.getVector("vrlDir", Vector(1.0,0.0,0.0));
            m_rayOrigin = props.getPoint("vrlOrigin", Point(0.0,0.0,0.0));
            m_vrlDistance = props.getFloat("vrlDistance", 10.0);
            m_nb_samples_vrl = props.getInteger("samples_vrl", 1);
            //m_solid_angle_sampling = props.getBoolean("solid_angle_sampling", false);

            if (m_maxDepth == 1) {
                Log(EInfo, "Enabled direct sampling of lights");
                m_directSampling = true;
            } else {
                Log(EInfo, "Disabled direct sampling of lights");
                m_directSampling = false;
            }

            /*if (m_solid_angle_sampling) {
                Log(EInfo, "Enabled solid angle sampling + direct sampling");
                m_directSampling = true;
            }*/
        }

    //First value : s, second value factor (1/probability)
	std::tuple<Float,Float> sample_distance(const Medium& medium, Float sample) const {
		Float sigma_ref = medium.getSigmaT().max();
        Float d = - std::log(std::max(1.e-20,sample))/sigma_ref;
        return std::tuple<Float,Float>(d, 1.0f/(sample*sigma_ref));
    }


    Float asinh(Float value) const
	{
		if(value>0)
			return log(value + sqrt(value * value + 1));
		else
			return -log(-value + sqrt(value * value + 1));
	}

    // Sample the Virtual Ray Light assuming an isotropic medium. Note that
    // this function samples according to the attenuation IN 3D.
    Float sample_vrl( const Float A0, const Float A1, const Float h,
        const Float sinthita, Float &pdf, Float sample) const
    {
        // Compute vt according to Eq. 13 in [1].
        Float epsilon = std::max(1.e-20,sample);
        Float vt = h*std::sinh(A1*epsilon+A0*(1.f-epsilon))/sinthita;

        // Compute the pdf according to Eq.9 in [1].
        // This is equivalent to (Eq.10)/(Eq.11) in [1]. Note that
        // the terms have been already simplified.
        pdf = sinthita / (sqrtf(h*h+vt*vt*sinthita*sinthita)*(A1-A0));

        return vt;
    }

    // Sample the camera ray assuming an isotropic medium, using Kulla and Fajardo's
    // equiangular sampling [2]. Note that this function samples according to the
    // attenuation IN 3D.
    Float sample_camera_ray( const Float B0, const Float B1, const Float h, Float &pdf, Float sample) const
    {
        // Compute vt according to Eq. 14 in [1], which bases
        // on Eq. 6 in [2].
        Float epsilon = std::max(1.e-20,sample);
        Float ut = h*std::tan(B1*epsilon+B0*(1.f-epsilon));

        // Compute the pdf according to a modified version of
        // Eq.5 in [2].
        pdf = h / ((B1-B0)*(h*h + ut*ut));

        return ut;
    }

    //First value : s, second value factor (1/probability)
	std::tuple<Float,Float> sample_distance_equiangular(const Point lightPos, const Ray &ray, const Medium& medium, Float sample) const {
		// distance from ray to light closest point in ray
        Float t_light = dot(lightPos - ray.o, ray.d);

        // distance of point from light
        float distance = std::max((lightPos - (ray.o + t_light*ray.d)).length(),Epsilon);

        Float thetaA = std::atan(-t_light/distance);
        Float thetaB = std::atan(M_PI/2);

        Float tan_r = std::tan((1-sample)*thetaA + sample*thetaB);
        //Float r = std::atan(r);

        Float t = tan_r* distance;
        return std::tuple<Float,Float>(t_light+t, distance*(t_light/distance + M_PI/2));
    }

    // Auxiliary functions needed to sample the VRL
	Float rays_distance( const Ray &r1, const Ray &r2, Float &t1, Float &t2 )const
	{
		Float a = dot(r1.d,r1.d);
		Float b = dot(r1.d,r2.d);
		Float c = dot(r2.d,r2.d);
		Float d = dot(r1.d,r1.o-Point(0.0));
		Float e = dot(r2.d,r2.o-Point(0.0));

		t1 = (b*e-c*d)/(a*c-b*b);
		t2 = (a*e-b*d)/(a*c-b*b);

		Vector h =	r1.o + t1*r1.d -
						r2.o + t2*r2.d;

		return h.length();
	}

    Spectrum surface_hit_factor(const Medium& medium, Float dist) const {
		return (-dist*medium.getSigmaT()).exp()/(1.0f - std::exp(-dist*medium.getSigmaT().max()));
    }


    /// Unserialize from a binary data stream
	SingleScatteringQuadratureVRL(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;
        Float m_max_dist = m_vrlDistance;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (rRec.medium) {
                // Compute single scattering with vrl radiance
                const Medium* medium = rRec.medium;

                //TODO Make parameter
                RayDifferential vritualRay(m_rayOrigin, m_rayDir, ray.time);
                //rRec.rayIntersect(vritualRay);

                // Assuming a homogeneous medium
	            Spectrum u_t= medium->getSigmaT(),
			    u_s= medium->getSigmaS();

                // Get u_h, v_h, and h (Fig. 3 in [1])
                Float uh, vh;
                Float h = rays_distance(ray, vritualRay, uh, vh);

                // Compute u0t, u1t, v0t and v1t as the change of
                // variables explained in Sec. 4.1 in [1]).
                Float u0t = -uh, u1t = (its.t>m_max_dist?m_max_dist:its.t)-uh,
				v0t = -vh;

                RadianceQueryRecord vRec(rRec);
                vRec.rayIntersect(vritualRay);
				Float v1t = (vRec.its.t>m_max_dist?m_max_dist:vRec.its.t)-vh;
                //v1t = m_max_dist-vh;

                // Compute values needed to apply VRL and Camera Ray sampling
                // according to [1] and [2].
                Float costhita = dot(ray.d, vritualRay.d);
                Float sinthita = sqrt(1-costhita*costhita);
                Float Av0 = asinh(v0t/h*sinthita),
                    Av1 = asinh(v1t/h*sinthita);
                Float Bu0 = std::atan2(u0t,h),
                    Bu1 = std::atan2(u1t,h);

                //TODO Integrator parameter
                Spectrum flux = m_raySpectrum; // VRL Spectrum
			    //Float t_vrl = vritualRay.get_time_of_flight();

                // And then compute the contribution of the VRL by solving
                // Eq. 7 [1] using Eq. 8 [1].
                for( int s=0; s< m_nb_samples_vrl; ++s)
                {
                    Float pdf_v, pdf_u;

                    // Sample v according to Eq.13 in [1]. Note that we are
                    // assuming an isotropic medium.
                    Point2 samples = rRec.nextSample2D();
                    Float vt = sample_vrl(Av0, Av1, h, sinthita, pdf_v, samples[0]);

                    // Sample v according to Eq.14 in [1], formulated in [2].
                    Float ut = sample_camera_ray( Bu0, Bu1, h, pdf_u, samples[1]);

                    if(its.isValid() && ut+uh > its.t){
                        throughput *= ((-(ut+uh)*medium->getSigmaT()).exp()* (pdf_u*pdf_v));
                        continue;
                    }
                    // Evaluate function applying the integrated term in Eq. 7 in [1].
                    // Note that we need to unmake the variable changes performed before
                    // (Sec. 4.1 in [1]).
                    Point	u = r.o+r.d*(ut+uh),
                                v = vritualRay.o +
                                    vritualRay.d*(vt+vh);

                    // Compute distance
                    Float W = (u-v).length();
                    if(W < Epsilon){
                        W = 1.e20;
                    }
                    Vector v2u = normalize(u-v);

                    // Evaluate the visibility term
                    RayDifferential vray(v, v2u, ray.time);
                    vRec.rayIntersect(vray);
                    if (vRec.its.t< W+1.e-8){
                        continue;
                    }

                    //Spectrum Luv(1.f/(D==3?(W*W):W));
                    Spectrum Luv(1.f/(W*W));

                    // Compute phase functions
                    const PhaseFunction *phase = medium->getPhaseFunction();

                    // For Ray phase function
                    MediumSamplingRecord mRec; mRec.medium = medium; mRec.time = ray.time;
                    mRec.t = ut+uh; mRec.p = u;

                    // For VRL phase function
                    MediumSamplingRecord vRec; vRec.medium = medium; vRec.time = ray.time;
                    vRec.t = vt+vh; vRec.p = v;

                    Luv *= phase->eval(PhaseFunctionSamplingRecord(mRec, v2u, -ray.d))*
                        phase->eval(PhaseFunctionSamplingRecord(vRec, vritualRay.d, -v2u));

                    // Compute scattering term (assuming homogeneous medium)
                    Luv *= u_s*u_s;

                    // Finally, compute extinction (assuming homogeneous medium)
                    Float c2u2v2l = (ut+uh)+(vt+vh)+W;
                    Luv *= (-u_t*c2u2v2l).exp();

                    // And store the sample
                    Li += Luv*flux/(pdf_u*pdf_v*(Float)m_nb_samples_vrl);
                    return Li; //We never hit anything else.
                }

            }

            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    //Li += throughput * scene->evalEnvironment(ray);
                	Li = throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered)) {
                //Li += throughput * its.Le(-ray.d);
            	Li = throughput * its.Le(-ray.d);

            	break;
            }

#if 0
            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);
#endif

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);
            if (m_directSampling) {
                /*
                * Implement direct sampling of all emitters so they can contribute when quadrature without illum dim
                */
                if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                    (bsdf->getType() & BSDF::ESmooth)) {

                    const ref_vector<Emitter>& emitters = scene->getEmitters();
                    Spectrum value;
                    for (const ref<Emitter> & emitter : emitters) {

                        if( emitter->isEnvironmentEmitter() ||
                                (!emitter->isDegenerate() && emitter->getShape() == nullptr)) {
                            continue;
                        }

                        DirectSamplingRecord dRec2(its);
                        value = scene->sampleDirectPoint(dRec2, emitter.get());

                        if (!value.isZero()) {
                            /* Allocate a record for querying the BSDF */
                            BSDFSamplingRecord bRec(its, its.toLocal(dRec2.d), ERadiance);

                            /* Evaluate BSDF * cos(theta) */
                            const Spectrum bsdfVal = bsdf->eval(bRec);

                            /* Prevent light leaks due to the use of shading normals */
                            if (!bsdfVal.isZero() && (!m_strictNormals
                                    || dot(its.geoFrame.n, dRec2.d) * Frame::cosTheta(bRec.wo) > 0)) {

                                Spectrum medium_weight(1.0); //We account for the light source extintion
                                if (rRec.medium) medium_weight *= (-dRec.dist*rRec.medium->getSigmaT()).exp();
                                Li += throughput * value * bsdfVal * medium_weight;
                            }
                        }
                    }
                }
            }

            if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                    * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            //Spectrum bsdfEval = bsdf->eval(bRec);
            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */

            throughput *= bsdfWeight;
            //throughput *= bsdfEval;
            eta *= bRec.eta;
#if 1
            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (m_directSampling) {
                if (hitEmitter &&
                    (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {

                    /* Compute the prob. of generating that direction using the
                    implemented direct illumination sampling technique */
                    /*const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                        scene->pdfEmitterDirect(dRec) : 0;
                    Li += throughput * value * miWeight(bsdfPdf, lumPdf);*/
                    //Li = throughput * value;
                    break;
                }
            } else {
                if (hitEmitter &&
                                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                    Li += throughput * value;
                    break;
                }
            }

#endif

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            if (m_directSampling) {
                rRec.type = RadianceQueryRecord::ERadianceNoEmission;
            }

            //if (rRec.depth++ >= m_rrDepth) {
            if (rRec.depth++ >= 5000) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }

            // Update path length
            //rRec.depth++;
        }

        /* Store statistics */
        avgPathLength_singleScattering_vrl.incrementBase();
        avgPathLength_singleScattering_vrl += rRec.depth;

        // Try to use it to discard paths that doesn't contribute to the final result

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "PathTracer Quadrature[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    bool m_directSampling;
    Spectrum m_raySpectrum;
    Vector m_rayDir;
    Point m_rayOrigin;
    Float m_vrlDistance;
    int m_nb_samples_vrl;
    //bool m_solid_angle_sampling;

    MTS_DECLARE_CLASS()
};

MTS_NAMESPACE_END

#endif /* SRC_INTEGRATORS_PTRACER_PATHTRACER_H_ */
