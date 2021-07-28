
#ifndef SRC_INTEGRATORS_SINGLESCATTERING_H_
#define SRC_INTEGRATORS_SINGLESCATTERING_H_


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

#include "singleScatteringSampling.h"

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength_singleScattering("Single scattering quadrature", "Average path length", EAverage);

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
class SingleScatteringQuadrature : public MonteCarloIntegrator {
public:
	SingleScatteringQuadrature(const Properties &props)
        : MonteCarloIntegrator(props) {
            m_directSampling = props.getBoolean("directSampling", false);
            m_mapping_0_ts = props.getBoolean("mappingTs", false);
            m_solid_angle_sampling = props.getBoolean("solid_angle_sampling", false);

            if (m_maxDepth == 1) {
                Log(EInfo, "Enabled direct sampling of lights");
                m_directSampling = true;
            } else {
                Log(EInfo, "Disabled direct sampling of lights");
                m_directSampling = false;
            }

            if (m_solid_angle_sampling) {
                Log(EInfo, "Enabled solid angle sampling + direct sampling");
                m_directSampling = true;
            }
        }

    /// Unserialize from a binary data stream
	SingleScatteringQuadrature(Stream *stream, InstanceManager *manager)
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


        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (rRec.medium) {
                Float medium_sample = rRec.nextSample1D();

                // [0, Ts] Mapping
                Float mfp_ts_sample = SingleScatteringSampling::cdf_mfp(*rRec.medium, its.t);
                if(m_mapping_0_ts) medium_sample*=(mfp_ts_sample);

                auto [medium_interaction_distance, single_scattering_probability_factor] = SingleScatteringSampling::sample_distance_mfp(*rRec.medium, medium_sample);
                if ((!its.isValid()) || m_mapping_0_ts || (medium_interaction_distance < its.t)) {  //We hit with the medium so we calculate is contribution and NEVER HIT ANYTHING ELSE

                    throughput *= (-(medium_interaction_distance)*rRec.medium->getSigmaT()).exp();

                    // [0, Ts] Mapping
                    if(m_mapping_0_ts)throughput *= mfp_ts_sample;

                    if (m_solid_angle_sampling) {
                        const PhaseFunction *phase = rRec.medium->getPhaseFunction();
                        MediumSamplingRecord mRec; mRec.t = medium_interaction_distance; mRec.p = ray(mRec.t); mRec.medium = rRec.medium; mRec.time = ray.time;


                        DirectSamplingRecord dRec(mRec);
                        int interactions = m_maxDepth - rRec.depth;
                        Spectrum value = scene->sampleAttenuatedEmitterDirect(
                            dRec, its, rRec.medium, interactions,
                            rRec.nextSample2D(), rRec.sampler);

                        //Spectrum value = scene->sampleAttenuatedEmitterDirect(dRec, rRec.nextSample2D());
                        if (dRec.dist<=Epsilon) value = Spectrum(1.e20);
                        if (value.isZero()) return value;
                        else return single_scattering_probability_factor * value * throughput * rRec.medium->getSigmaS() * phase->eval(PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));
                    } else {
                       const ref_vector<Emitter>& emitters = scene->getEmitters();
                       const PhaseFunction *phase = rRec.medium->getPhaseFunction();
                        MediumSamplingRecord mRec; mRec.t = medium_interaction_distance; mRec.p = ray(mRec.t); mRec.medium = rRec.medium; mRec.time = ray.time;
                        for (const ref<Emitter> & emitter : emitters) {
                            if( emitter->isEnvironmentEmitter() || (!emitter->isDegenerate() && emitter->getShape() == nullptr)) {
                                continue;
                            }

                            DirectSamplingRecord dRec(mRec);
                            int interactions = m_maxDepth - rRec.depth - 1;

                            // Direct light sampling
                            Spectrum value = scene->sampleDirectPoint(dRec, emitter.get());

                            // Transmittance
                            value *= (-(dRec.dist)*rRec.medium->getSigmaT()).exp();
                            /*value *= scene->evalTransmittance(dRec.ref, false,
                            dRec.p, emitter->isOnSurface(), dRec.time, rRec.medium,
                            interactions, rRec.sampler);*/

                            if (dRec.dist<=Epsilon) value = Spectrum(1.e20);
                            if (!value.isZero()) {
                                Li += single_scattering_probability_factor * value * throughput * rRec.medium->getSigmaS() * phase->eval(PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));
                            }
                        }
                        /*if(!m_mapping_0_ts)
                            return Li;*/
                    }
                } else {
                    //throughput *= single_scattering_probability_factor;
                    //throughput *= 1.0;
                    //throughput *= surface_hit_factor(*rRec.medium, medium_interaction_distance);
                    //throughput *= (-its.t*rRec.medium->getSigmaT()).exp() / (1.0 - SingleScatteringSampling::cdf_mfp(*rRec.medium, its.t));
                    //throughput *= (-(its.t)*rRec.medium->getSigmaT()).exp()*single_scattering_probability_factor;
                    //throughput *= (-(its.t)*rRec.medium->getSigmaT()).exp()*single_scattering_probability_factor*std::exp(-(medium_interaction_distance - its.t));

                    // LATEST ONE
                    throughput *= (-its.t*rRec.medium->getSigmaT()).exp() / (1.0 - SingleScatteringSampling::cdf_mfp(*rRec.medium, its.t));
                }


                if(its.isValid()){
                    throughput = (-its.t*rRec.medium->getSigmaT()).exp();
                }else{
                    throughput *= 0;
                }
            }

            if (!its.isValid()) {
                /* If no intersection could be found, possibly return
                    attenuated radiance from a background luminaire */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered)) {
                    Spectrum value = throughput * scene->evalEnvironment(ray);
                    if (rRec.medium)
                        value *= rRec.medium->evalTransmittance(ray, rRec.sampler);
                    Li += value;
                }

                break;
            }


            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered)) {
                //Li += throughput * its.Le(-ray.d);
            	Li += throughput * its.Le(-ray.d);

            	break;
            }

#if 0
            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

            if (rRec.depth >= m_maxDepth && m_maxDepth != -1)
                break;

            /* Prevent light leaks due to the use of shading normals */
            Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
                    wiDotShN  = Frame::cosTheta(its.wi);
            if (wiDotGeoN * wiDotShN < 0 && m_strictNormals)
                break;
#endif

            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            const BSDF *bsdf = its.getBSDF(ray);
            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);
            if (m_directSampling) {

                if (m_solid_angle_sampling) {
                    if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                    (bsdf->getType() & BSDF::ESmooth)) {
                        int interactions = m_maxDepth - rRec.depth - 1;
                        Spectrum value = scene->sampleAttenuatedEmitterDirect(
                            dRec, its, rRec.medium, interactions,
                            rRec.nextSample2D(), rRec.sampler);
                        //Spectrum value = scene->sampleAttenuatedEmitterDirect(dRec, rRec.nextSample2D());
                        if (!value.isZero()) {
                            const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                            /* Allocate a record for querying the BSDF */
                            BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                            /* Evaluate BSDF * cos(theta) */
                            const Spectrum bsdfVal = bsdf->eval(bRec);

                            /* Prevent light leaks due to the use of shading normals */
                            if (!bsdfVal.isZero() && (!m_strictNormals
                                    || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                                /* Calculate prob. of having generated that direction
                                using BSDF sampling */
                                Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                    ? bsdf->pdf(bRec) : 0;

                                /* Weight using the power heuristic */
                                Float weight = miWeight(dRec.pdf, bsdfPdf);
                                //Li += throughput * value * bsdfVal * weight;
                                Spectrum medium_weight(1.0); //We account for the light source extintion
                                if (rRec.medium) medium_weight *= (-dRec.dist*rRec.medium->getSigmaT()).exp();
                                Li += throughput * value * bsdfVal * medium_weight;
                            }
                        }
                    }
                } else {
                    /*
                    * Implement direct sampling of all emitters so they can contribute when quadrature without illum dim
                    */

                    /* Estimate the direct illumination if this is requested */
                    if ((rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
                        (bsdf->getType() & BSDF::ESmooth)) {
                        const ref_vector<Emitter>& emitters = scene->getEmitters();
                        for (const ref<Emitter> & emitter : emitters) {
                            if( emitter->isEnvironmentEmitter() ||
                                    (!emitter->isDegenerate() && emitter->getShape() == nullptr)) {
                                continue;
                            }

                            int interactions = m_maxDepth - rRec.depth - 1;

                            // Direct light sampling
                            DirectSamplingRecord dRec2(its);
                            Spectrum value = scene->sampleDirectPoint(dRec2, emitter.get());

                            // Transmittance
                            value *= (-(dRec.dist)*rRec.medium->getSigmaT()).exp();
                            /*const Medium *medium = rRec.medium;
                            if (its.shape && its.isMediumTransition())
                                medium = its.getTargetMedium(dRec2.d);

                            value *= scene->evalTransmittance(its.p, true, dRec2.p, emitter->isOnSurface(),
                            dRec2.time, medium, interactions, rRec.sampler);*/


                            if (!value.isZero()) {
                                /* Evaluate BSDF * cos(theta) */
                                BSDFSamplingRecord bRec(its, its.toLocal(dRec2.d), ERadiance);

                                /* Evaluate BSDF * cos(theta) */
                                const Spectrum bsdfVal = bsdf->eval(bRec);

                                /* Prevent light leaks due to the use of shading normals */
                                if (!bsdfVal.isZero() && (!m_strictNormals ||
                                    dot(its.geoFrame.n, dRec2.d) * Frame::cosTheta(bRec.wo) > 0)) {

                                    /* Weight using the power heuristic */
                                    Spectrum medium_weight(1.0); //We account for the light source extintion
                                    Li += throughput * value * bsdfVal * medium_weight;
                                }
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
        avgPathLength_singleScattering.incrementBase();
        avgPathLength_singleScattering += rRec.depth;

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
    bool m_mapping_0_ts;
    bool m_solid_angle_sampling;

    MTS_DECLARE_CLASS()
};

MTS_NAMESPACE_END

#endif /* SRC_INTEGRATORS_PTRACER_PATHTRACER_H_ */
