# Transient Wave Dynamic Acousto Elastic Testing (TW-DAET)
### Dependence of Elastic Nonlinearity on Aligned Inhomogeneities
The following GitHub repository collects data, documentation, programs, and figures put together during Jacob Newman's M.Sc. (geophysics) research project. The thesis, in its entirety, can be found in the "Thesis" directory; we provide an abstract below.

We have shared this repository and all its contents to encourage Open Science to drive further the continuation of research on elastic nonlinearity in geomaterials and incentive Open Science in general. We provide all the data, programs, and figures included in the thesis here.

## Abstract
The nonlinear elasticity of geomaterials is a prominent indicator of the materials' internal structure (micro-fractures, grain-to-grain boundaries, etc.), both microscopically and macroscopically. Dynamic acousto-elastic testing (DAET) has become a standard tool to assess the nonlinear elasticity of materials. Here, we use a transient wave variation of DAET (TW-DAET) to evaluate whether this type of testing can (1) distinguish the presence of varying types of inhomogeneities and (2) determine the direction in which said inhomogeneities align. We explore two main types of inhomogeneities with TW-DAET; micro-fracturing in sandstone samples (which we attempt to induce via freeze-thaw cycles) and embedded objects within an elastically linear background cement. The sandstone results remain largely inconclusive though evidence for a potential two-mechanism (one increasing and one decreasing nonlinearity) system presents itself in a subset of the data. We can observe the presence and alignment of inhomogeneities for a part of the cement samples we examine with TW-DAET. We also find evidence for a potential nonlinearity dependence on the relationship between the wavelength of the probe and the spatial characteristics of the inhomogeneities.

## Experimental Design
The nonlinear testing method we employ for this work follows a pump-probe dynamic acousto-elastic testing design in which the perturbation field (i.e. the pump) is a transient wave. Any experiment to test the classical elastic nonlinearity of a medium should aim to characterize the elastic nonlinearity coefficient(s) which scale the nonlinear term(s) of the wave equation. Here we simplify this process by using an indirect measure of a material's classical elastic nonlinearity, namely, a change in wave speed and thus time-of-flight across one of our samples induced by an external strain perturbation (pump wave). We can derive the material, classical elastic nonlinearity by comparing the probe signal with and without the pump wave's presence.

The experimental design we use is shown below. We use an arbitrary waveform generator (Keysight 33500B Series) to excite both the pump and the probe signals. We choose the probe wave to be a low-amplitude (5 V), high-frequency (500 kHz) signal relative to the pump wave's high-amplitude (10 V amplified 50x by a TEGAM HIGH Voltage Amplifier), low-frequency (50 kHz) character. We use a relatively low amplitude for the probe wave to preserve the medium's elastic parameters during the probe's transmission (i.e., the source of perturbation is due exclusively to the pump wave and not the probe). We use a high amplitude for the pump to induce a large enough perturbation to be sensed by the probe wave. We choose the probe and pump frequencies such that the pump wavelength is significantly larger (~ x10) than that of the probe wavelength. This choice results in an approximately steady-state pump wave at the timescale of the probe, thus granting us the ability to sense different phases of the pump waveform with the high-frequency probe. 
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Experimental_Setup.jpg" width = 600>
</p>

<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Orientation_diagrams_Rev3.PNG" width = 600>
</p>
          
<p align="center">          
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/All_cement_samples.png" width = 600>
</p>
          
          
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Sand_layers_and_pump.png" width = 600>
</p>
          
          
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/wires_and_pump.png" width = 600>
</p>
          
          
