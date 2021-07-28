#pragma once

#include "../../../../ext/viltrum/quadrature/integrate-bins-adaptive.h"

using namespace viltrum;

template<typename Nested, typename Error, typename RNG>
class AdaptiveResidualRatioTracking {
    mutable RNG rng;
    double maj_sigma_t;
    FILE *m_f;


    StepperAdaptive<Nested,Error> cv_stepper;
    unsigned long adaptive_iterations;

    template<typename R, typename Result>
    struct Data {
        std::vector<R> regions;

        Result Tc;
        Result sumatory;
        unsigned long counter;
        unsigned long counter_queries;

        Data( std::vector<R>&& rs ) :
            regions(std::forward<std::vector<R>>(rs)),
            sumatory(0),counter(0),counter_queries(0) { }
    };
public:
    void set_file(const char* name_file)
    {
        if(m_f)
            fclose(m_f);

        m_f = fopen(name_file, "w");
    }

    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
        auto regions = cv_stepper.init(f,range);
        using VECTOR_TYPE = typename decltype(regions)::value_type;

        return Data<VECTOR_TYPE, decltype(f(std::array<Float,1>{{range.min(0)}}))> (std::move(regions));
    }

    template<typename F, typename Float, typename R>
    void step(const F& f, const Range<Float,1>& range, Data<R,Float>& data, bool verbose = false) const
    {
        // Compute the polynomial approximation
        if (!data.counter) {
            for(unsigned int i=0; i<adaptive_iterations; ++i)
                cv_stepper.step(f,range,data.regions);

            data.Tc = exp(-cv_stepper.integral(f,range,data.regions));

            data.counter_queries = f.samples();

            if( verbose )
                printf("Computing Polynomial Approximation: %d iterations -- %d samples\n", adaptive_iterations, data.counter_queries);

            if(m_f)
            {
                Float dt = (range.max(0)-range.min(0))/(999.);
                for(Float t=range.min(0); t<= range.max(0); t+= dt)
                    fprintf(m_f,"%f %f\n", t, f(std::array<Float,1>{{t}}));

                fprintf(m_f,"n\n");
                for(Float t=range.min(0); t<= range.max(0); t+= dt)
                {
                    int index = -1;
                    for(unsigned int i=0; i<data.regions.size(); ++i)
                        if(data.regions[i].range().min()[0] <= t && data.regions[i].range().max()[0] >= t)
                        {
                            index = i; break;
                        }
                    if( index < 0)
                    {
                        printf("Outside of range! t=%f", t);
                        continue;
                    }
                    fprintf(m_f,"%f %f\n", t, data.regions[index].approximation_at(t));
                }
                fprintf(m_f,"n\n");

                fprintf(m_f,"%f %f\n%f %f\n", range.min(0), maj_sigma_t, range.max(0), maj_sigma_t);
                fprintf(m_f,"n\n");


                // Transmittance computed with ray marching (biased)
                Float delta_t = (range.max(0)-range.min(0))/(Float)10000;

                Float Tr = 1, Tc = 1, T = 1;
                auto t = range.min(0)+delta_t/2.;

                while(true)
                {

                    if( t>range.max(0) )
                        break;

                    T *= exp(-f(std::array<Float,1>{{t}})*delta_t);

                    int index = -1;
                    for(unsigned int i=0; i<data.regions.size(); ++i)
                        if(data.regions[i].range().min()[0] <= t && data.regions[i].range().max()[0] >= t)
                        {
                            index = i; break;
                        }
                    if( index < 0)
                    {
                        printf("Outside of range! t=%f", t);
                        continue;
                    }

                    Tc *= exp(-data.regions[index].approximation_at(t)*delta_t);
                    Tr *= exp(-(f(std::array<Float,1>{{t}})-data.regions[index].approximation_at(t))*delta_t);

                    fprintf(m_f, "%f %f %f %f\n", t, T, Tc, Tr);

                    t += delta_t;
                }
                fprintf(m_f,"s\n");


            }
        }

        std::exponential_distribution<Float> dis(maj_sigma_t);

        Float Tr = 1;
        Float t = range.min(0);


        if( m_f )
        {
            fprintf(m_f,"%f %f\n", t, Tr);
        }

        while(true)
        {
            Float dt = dis(rng);
            t += dt;

            if( t>range.max(0) )
                break;

            unsigned int index=-1;
            for(unsigned int i=0; i<data.regions.size(); ++i)
            {
                if(data.regions[i].range().min()[0] < t && data.regions[i].range().max()[0] >= t)
                {
                    index = i; break;
                }
            }

            if(index < 0)
            {
                printf("Error - Outside of the quadrature reconstruction!"); return;
            }

            Float c_sigma_t = data.regions[index].approximation_at(t);
            Float ft = f(std::array<Float,1>{{t}});


            if( m_f )
            {
                fprintf(m_f,"%f %f\n", t, Tr);
            }

            Tr *= (1 - (ft - c_sigma_t)/maj_sigma_t);

            if( m_f )
            {
                fprintf(m_f,"%f %f\n", t, Tr);


            }
            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f; Tc(t) = %f; T_c = %f\n", t, range.min(0), range.max(0), (ft-c_sigma_t)/maj_sigma_t, Tr, c_sigma_t, data.Tc);

            ++ data.counter_queries;
        }

        if( m_f )
        {
            fprintf(m_f,"%f %f\n", range.max(0), Tr);
            fprintf(m_f,"s\n");

        }

        data.sumatory += Tr*data.Tc;
        ++data.counter;

    }

    template<typename F, typename Float, typename R>
    Float integral(const F& f, const Data<R,Float>& data) const {
        return (data.counter==0)?decltype(data.sumatory)(0):(data.sumatory/data.counter);
    }


    template<typename F, typename Float, typename R>
    unsigned long queries(const F& f, const Data<R,Float>& data) const {
        return (data.counter==0)?decltype(data.counter_queries)(0):(data.counter_queries);
    }









    AdaptiveResidualRatioTracking(
        Nested&& nested, Error&& error, RNG&& r, unsigned long ai, double majorant) :
            cv_stepper(std::forward<Nested>(nested), std::forward<Error>(error)), adaptive_iterations(ai),
            rng(std::forward<RNG>(r)), maj_sigma_t(majorant), m_f(0) { }


    ~AdaptiveResidualRatioTracking()
    {
        if(m_f)
            fclose(m_f);

        m_f=0;
    }
};

template<typename Nested, typename Error, typename RNG>
auto adaptive_residual_ratio_tracking(Nested&& nested, Error&& error, unsigned long adaptive_iterations,
                                      RNG&& rng, double majorant ) {
    return AdaptiveResidualRatioTracking<std::decay_t<Nested>,std::decay_t<Error>, std::decay_t<RNG>>(
        std::forward<Nested>(nested),std::forward<Error>(error), std::forward<RNG>(rng), adaptive_iterations, majorant);
}

template<typename Nested, typename RNG>
auto adaptive_residual_ratio_tracking(Nested&& nested, unsigned long adaptive_iterations, RNG&& rng, double majorant)
{
    return adaptive_residual_ratio_tracking(std::forward<Nested>(nested), error_single_dimension_standard(), adaptive_iterations, std::forward<RNG>(rng), majorant);
}






// #pragma once

// #include "integrate-bins-adaptive.h"
// #include "monte-carlo.h"
// #include "sample-vector.h"

// template<typename Nested, typename Error, typename ResidualStepper, typename VectorSampler>
// class AdaptiveResidualRatioTracking {
//     StepperAdaptive<Nested,Error> cv_stepper;
//     ResidualStepper residual_stepper;
//     VectorSampler vector_sampler;
//     unsigned long adaptive_iterations;

//     template<typename R,typename ResData,typename Sampler>
//     struct Data {
//         std::vector<R> regions;
//         unsigned long cv_iterations;

//         Data(std::vector<R>&& rs, ResData&& rd, Sampler&& vs) :
//             regions(std::forward<std::vector<R>>(rs)),
//             residual_data(std::forward<ResData>(rd)),
//             vector_sampler(std::forward<Sampler>(vs)),
//             cv_iterations(0) { }
//     };
// public:
//     template<typename F, typename Float, std::size_t DIM>
//     auto init(const F& f, const Range<Float,DIM>& range) const {
//         auto regions = cv_stepper.init(f,range);
//         //residual_stepper.init should not do any calculation at all (should be MC)
//         // return Data(std::move(regions),
//         //          residual_stepper.init(f, range),
//         //          vector_sampler(regions));

//         auto init = residual_stepper.init(f, range);

//         using VECTOR_TYPE = typename decltype(regions)::value_type;

//         return Data<VECTOR_TYPE,
//                     decltype(init),
//                     decltype(vector_sampler(regions))>
//                     (std::move(regions),
//                     std::move(init),
//                     vector_sampler(regions));
//     }

//     template<typename F, typename Float, std::size_t DIM, typename R, typename Sampler>
//     void step(const F& f, const Range<Float,DIM>& range, Data<R,ResData,Sampler>& data) const
//     {

//         // Compute the polynomial approximation
//         if (!data.cv_iterations) {

//             while (data.cv_iterations < adaptive_iterations) {
//                 cv_stepper.step(f,range,data.regions);
//                 ++data.cv_iterations;
//             }
//             Tc = exp(-cv_stepper.integral(f,range,data.regions));
//         }

//         ++data.cv_iterations;
//         std::exponential_distribution<Float> dis(maj_sigma_t);

//         Float Tr = 1;
//         Float t = range.min(0);

//         while(true)
//         {
//             Float dt = dis(rng);
//             t += dt;

//             if( t>range.max(0) )
//                 break;

//             unsigned int index=-1;
//             for(unsigned int i=0; i<data.regions.size(); ++i)
//             {
//                 if(data.regions[i].range.min() < t && data.regions[i].range.max() => t)
//                 {
//                     index = i; break;
//                 }
//             }

//             if(index < 0)
//             {
//                 printf("Error - Outside of the quadrature reconstruction!"); return 0;
//             }

//             Float c_sigma_t = data.regions[index].approximation_at(t);

//             Tr *= (1 - (f(t) - c_sigma_t)/maj_sigma_t);

//             if( verbose )
//                 printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f; Tc(t) = %f; T_c = %f\n", t, range.min(0), range.max(0), (f(t)-c_sigma_t)/maj_sigma_t, Tr, exp(-c_sigma_t*t), Tc);
// //                printf("Step: %f [%f, %f] - mu_t(t)= %f : T(t) = %f \n", t, range.min(0), range.max(0), f(t)-c_sigma_t, Tr*exp(-t*c_sigma_t));

//             ++ samples.counter_queries;
//         }


//         samples.sumatory += Tr*Tc;
//         ++samples.counter;

//     }

//     template<typename F, typename Float, std::size_t DIM, typename R,typename ResData,typename Sampler>
//     auto integral(const F& f, const Range<Float,DIM>& range, const Data<R,ResData,Sampler>& data) const {
//         return cv_stepper.integral(f,range,data.regions) +
//                 residual_stepper.integral(f,range,data.residual_data);
//     }

//     AdaptiveResidualRatioTracking(
//         Nested&& nested, Error&& error, ResidualStepper&& rs, VectorSampler&& vs, unsigned long ai) :
//             cv_stepper(std::forward<Nested>(nested), std::forward<Error>(error)),
//             residual_stepper(std::forward<ResidualStepper>(rs)),
//             vector_sampler(std::forward<VectorSampler>(vs)),
//             adaptive_iterations(ai) { }
// };

// template<typename Nested, typename Error, typename ResidualStepper, typename VectorSampler>
// auto stepper_adaptive_control_variates(Nested&& nested, Error&& error, ResidualStepper&& residual_stepper, VectorSampler&& vector_sampler, unsigned long adaptive_iterations) {
//     return AdaptiveResidualRatioTracking<std::decay_t<Nested>,std::decay_t<Error>,std::decay_t<ResidualStepper>,std::decay_t<VectorSampler>>(
//         std::forward<Nested>(nested),std::forward<Error>(error),std::forward<ResidualStepper>(residual_stepper),std::forward<VectorSampler>(vector_sampler),adaptive_iterations);
// }

// template<typename Nested>
// auto stepper_adaptive_control_variates(Nested&& nested, unsigned long adaptive_iterations, std::size_t seed_mc = std::random_device()(), std::size_t seed_vs = std::random_device()()) {
//     return stepper_adaptive_control_variates(std::forward<Nested>(nested),
//     error_single_dimension_standard(),stepper_monte_carlo_uniform(seed_mc),vector_sampler_uniform(seed_vs),adaptive_iterations);
// }

