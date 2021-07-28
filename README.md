
<!-- PROJECT LOGO -->
<br />
<p align="center">

  <h1 align="center"><a href="https://mcrescas.github.io/publications/primary-space-cv/">Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations</a></h1>

  <a href="https://mcrescas.github.io/publications/primary-space-cv/">
    <img src="https://mcrescas.github.io/publications/primary-space-cv/figures/socialMedia.png" alt="Logo" width="100%">
  </a>

  <p align="center">
    ACM Transactions on Graphics - 2021
    <br />
    <a href="https://mcrescas.github.io"><strong>Miguel Crespo</strong></a>
    ·
    <a href="http://giga.cps.unizar.es/~ajarabo/"><strong>Adrian Jarabo</strong></a>
    ·
    <a href="http://adolfo-munoz.com/"><strong>Adolfo Muñoz</strong></a>
  </p>

  <p align="center">
    <a href='https://mcrescas.github.io/publications/primary-space-cv/data/crespo2021primary.pdf'>
      <img src='https://img.shields.io/badge/Paper-PDF-red?style=flat-square' alt='Paper PDF'>
    </a>
    <a href='https://mcrescas.github.io/publications/primary-space-cv' style='padding-left: 0.5rem;'>
      <img src='https://img.shields.io/badge/Project-Page-blue?style=flat-square' alt='Project Page'>
    </a>
    <a href='https://github.com/adolfomunoz/viltrum' style='padding-left: 0.5rem;'>
      <img src='https://img.shields.io/badge/VILTRUM-Lib-green?style=flat-square' alt='VILTRUM Lib'>
    </a>
  </p>
</p>

<br />
<br />

<!-- TABLE OF CONTENTS -->
<details open="open" style='padding: 10px; border-radius:5px 30px 30px 5px; border-style: solid; border-width: 1px;'>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
      <ul>
        <li><a href="#citation">Citation</a></li>
        <li><a href="#project-structure">Project structure</a></li>
        <li><a href="#parameters-explanation">Parameters explanation</a></li>
      </ul>
    </li>
    <li>
      <a href="#compilation">Compilation</a>
      <ul>
        <li><a href="#dependencies">Dependencies</a></li>
        <li><a href="#how-to">How to</a></li>
      </ul>
    </li>
    <li>
      <a href="#results">Results</a>
      <ul>
        <li><a href="#transmittance-estimation">Transmittance estimation</a></li>
        <li><a href="#low-order-single-scattering">Low-order single scattering</a></li>
        <li><a href="#direct-illumination">Direct illumination</a></li>
        <li><a href="#distribution-effects">Distribution effects</a></li>
        <li><a href="#higher-order-integrals">Higher-order integrals</a></li>
      </ul>
    </li>
  </ol>
</details>
<br />
<br />

> **⚠ WARNING: Cloning the repository** <br>
> Our library [**VILTRUM**: **V**aried **I**ntegration **L**ayouts for arbi**TR**ary integrals in a **U**nified **M**anner](https://github.com/adolfomunoz/viltrum) is defined as a *submodule*. You need to use the following command to fetch it correctly: <br>
&nbsp; &nbsp; &nbsp; &nbsp; `git clone --recursive https://github.com/mcrescas/viltrum-mitsuba.git` <br>
> If you already cloned the repository and forgot to specify this flag, it’s possible to fix the repository in retrospect using the following command: <br>
&nbsp; &nbsp; &nbsp; &nbsp; `git submodule update --init --recursive`

<br>

## Overview

This repository contains the source code of the paper **Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations** by *Miguel Crespo*, *Adrian Jarabo*, and *Adolfo Muñoz* from ACM Transactions on Graphics.

The implementation is based on [Mitsuba 0.6](https://github.com/mitsuba-renderer/mitsuba), see the README in its repo for more information.

### Citation

```bibtex
@article{crespo21primary,
  title = {Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations.},
  year = {2021},
  journal = {ACM Transactions on Graphics},
  author = {Crespo, Miguel and Jarabo, Adrian and Mu\~{n}oz, Adolfo},
  volume = {40},
  number = {3},
  issn = {0730-0301},
  url = {https://doi.org/10.1145/3450627},
  doi = {10.1145/3450627},
  issue_date = {July 2021},
  month = jul,
  articleno = {25},
  numpages = {15},
}
```

### Project structure

Our technique is mainly contained as custom plugins inside Mitsuba, which are built on top of our custom library [**VILTRUM**: **V**aried **I**ntegration **L**ayouts for arbi**TR**ary integrals in a **U**nified **M**anner](https://github.com/adolfomunoz/viltrum). Note that a few modifications to the source code of Mitsuba were required, so our plugins are not completely independent.

* Our agnostic custom library deals with all the arithmetic required by our system, and can be found in `mitsuba/ext/viltrum` folder.
* Our integration with Mitsuba can be found inside `mitsuba/src/integrators/quad` folder.
    *  We have develop a main entry point for all of our specific integrators inside `quadrature.cpp` file, which defines an integrator called `quad` that launch our technique.
    *  Specific code for each type of integral can be found inside the corresponding file in that folder. This includes:
        1. Our modified **Path Tracer** inside `pathTracer.(h|cpp)`.
        2. Our modified **Volumetric Path Tracer** inside `singleScattering(...).(h|cpp)`.
        3. Our modified **Heterogeneous media** transmittance estimation inside `heterogeneous.cpp`.
        4. Our different transmittance estimation techniques inside `tracking` folder.
* Several extra utilities can be found inside `mitsuba/src/integrators/quad` folder.
    1. Sampler used to pass to the different integrators the evaluation points in the hyper cube inside `quadSampler.(h|cpp)`.
    2. Modified **area light** that supports defining radiance using a texture inside `areaColor.cpp`.
    3. Modified **independent sampler** that uses a random sampler, generating a different noise pattern while rendering using only one thread inside `independent.cpp`.
    4. Integration of [**OpenVDB**](https://www.openvdb.org/) in Mitsuba 0.6 inside `vdbvolume.(h|cpp)` (disabled by default).
* Other modifications of the original code of Mitsuba are related with changes in the interface of its components (e.g: lights interface or monte carlo integrator interface)

The idea of our integration with Mitsuba 0.6 is the following: because our technique cannot be use with the *sampling interface*, we need to take care of all the life cycle of an standalone integrator. In a basic way, our library in `mitsuba/ext/viltrum` works by using a generic function `F`. Our integration with Mitsuba is designed to provide that function (`callback_samplePoint (const Array& values, RNG& rng, bool customRNG)`) and save the result into disk. *We have two interfaces in our custom library (integrator and stepper). This is the main reason of our double implementation of `auxRenderInternalFinal(STEPPER& integrator, RESOLUTION& resolutionBins, std::string pathResult="") `*

Indeed, `quadrature` plugin launch a custom "integrator" of our library `mitsuba/ext/viltrum` while instantiating in addition an integrator of Mitsuba. While most of our parameters can be left as default for all scenes, there are a few of ones that need to be taken into account depending on the problem.

### Parameters explanation

#### Quad plugin

| Name  | Explanation  |
|---|---|
| spp  | Total budget of samples used by the algorithm |
| n_dims  | Number of dimensions for which the control variate will be built  |
| scaleSpp  | *(1/scaleSpp)\*spp* samples of the total will be used to build the control variate. The rest will be used to compute the residual |
| maxDepth | Maximum depth of the *subintegrator*. If different from maxDepthQuad, it means the number of higher dimensions computed with Monte Carlo.
| maxDepthQuad | Depth at which the control is to be constructed. (Optional: fallback maxDepth value)
| typeIntegrator | Name of the integrator to be used
| typeSubIntegrator | Name of the Mitsuba integrator to be used with our technique
|higherQuadRule| What quadrature rule to use as the higher one : `simpson` or `boole`. Defaults to `simpson`
|error_size_weight| Factor of our heuristic that depends on the size of the hyper region. Defaults to `0.00001`

#### Type integrators
| Name | Explanation |
|---|---|
| quad | Only adaptive quadrature without residual integration
|cv_pixel_alphaOpt| Our full technique featuring adaptive quadrature, control variate and optimizing alpha
|munoz2014| [[Muñoz  2014]](http://webdiis.unizar.es/~amunoz/projects/CGF2014_higherorder/) implementation
|mc| Monte Carlo integrator


#### Type Sub Integrator
| Name | Explanation |
|---|---|
|path| Path tracer
|singleScattering| Volumetric Path tracer
|singleScatteringEquiangular| Volumetric Path tracer using equiangular sampling
|singleScatteringVRL| Volumetric Path tracer using Virtual Ray Lights [[Novák 2012]](https://cs.dartmouth.edu/wjarosz/publications/novak12vrls.html)


#### Heterogeneous media transmittance estimator
| Name | Explanation |
|---|---|
|woodcock| Woodcock tracking
|simpson| Simpson quadrature
|ratiotracking| Ratio tracking
|residualratiotracking| Residual Ratio tracking
|adaptivesimpson| Adaptive Simpson quadrature
|adaptiveboole| Adaptive Boole quadrature
|cvadaptive| Our Adaptive Residual Ratio tracking

----

## Compilation

For specific details about the compilation process, we refer to the [documentation](http://mitsuba-renderer.org/docs.html) of Mitsuba 0.6.

Our modifications are integrated with the build system of Mitsuba, so no specific steps are necessary. We have tested our implementation in `Ubuntu Linux` using `GCC-9`.

### Dependencies

We have integrated [OpenVDB](https://www.openvdb.org/) inside this version of Mitsuba (disabled by default):
* If you want to use it, see `config.py` file in the root of the project and fill the path to each component required of OpenVDB and Intel TBB.
* Additionally, go to `src/integrators/Sconscript` file and uncomment line `31` refering to our `vdbvolume` plugin.

Furthermore, dependencies of Mitsuba 0.6 are required.

### How to

Go to `mitsuba` folder and type `scons --parallelize`. This will launch the building process using all the cores available in the computer.


----

## Results

We have a few of our scenes ready for testing our technique. The following sections have examples of how to use our system in several different types of problems.

To compute the results, launch in each folder the following command `mitsuba -p 1 scene_paper.xml`. Notice that our implementation currently **only works in a single thread**.

### Transmittance estimation

Our technique is integrated inside the transmittance estimation of the heterogenous plugin in Mitsuba. There is no limitation in the scene, as long as our plugin `heterogeneousQuad` is used.An example of scene can be found in `scenes/hetvol`.

### Low-order single scattering

Our single scattering implementation is integrated inside our mega plugin `quad`. An example of scene can be found inside `scenes/pumpkin`.

### Direct illumination

Our direct illumination implementation is integrated inside our mega plugin `quad`. An example of scene can be found inside `scenes/dragon`.

### Distribution effects

Our distribution effects implementation is integrated inside our mega plugin `quad`. An example of scene can be found inside `scenes/chess` which features depth of field.
Note that in this case the number of dimensions should be `3` when dealing with motion blur, `4` while dealing with depth of field or `5` while dealing with both effects.


### Higher-order integrals

Our framework can be used to evaluate higher-dimensional integrals while keeping our control variate working in the lower ones. An example can be found in `scenes/chess`. Note that the rendering command is `mitsuba -p 1 scene_paper_higher.xml`.
