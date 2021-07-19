# Transient Wave Dynamic Acousto Elastic Testing (TW-DAET)
### Dependence of Elastic Nonlinearity on Aligned Inhomogeneities
The following GitHub repository collects data, documentation, programs, and figures put together during Jacob Newman's M.Sc. (geophysics) research project. The thesis, in its entirety, can be found in the "Thesis" directory; we provide an abstract below.

We have shared this repository and all its contents to encourage Open Science to drive further the continuation of research on elastic nonlinearity in geomaterials and incentive Open Science in general. We provide all the data, programs, and figures included in the thesis here.

## Abstract
The nonlinear elasticity of geomaterials is a prominent indicator of the materials' internal structure (micro-fractures, grain-to-grain boundaries, etc.), both microscopically and macroscopically. Dynamic acousto-elastic testing (DAET) has become a standard tool to assess the nonlinear elasticity of materials. Here, we use a transient wave variation of DAET (TW-DAET) to evaluate whether this type of testing can (1) distinguish the presence of varying types of inhomogeneities and (2) determine the direction in which said inhomogeneities align. We explore two main types of inhomogeneities with TW-DAET; micro-fracturing in sandstone samples (which we attempt to induce via freeze-thaw cycles) and embedded objects within an elastically linear background cement. The sandstone results remain largely inconclusive though evidence for a potential two-mechanism (one increasing and one decreasing nonlinearity) system is presented in a subset of the data. We observe the presence and alignment of inhomogeneities for a sub-set of the cement samples we examine with TW-DAET. We also find evidence for a potential dependence of the nonlinearity on the relative sizes of the probe's wavelength and the spatial characteristics of the inhomogeneities.

## Experimental Design
The nonlinear testing method we employ for this work follows a pump-probe dynamic acousto-elastic testing design in which the perturbation field (i.e. the pump) is a transient wave. Any experiment to test the classical elastic nonlinearity of a medium should aim to characterize the elastic nonlinearity coefficient(s) which scale the nonlinear term(s) of the wave equation. Here we simplify this process by using an indirect measure of a material's classical elastic nonlinearity, namely, a change in wave speed and thus time-of-flight across one of our samples induced by an external strain perturbation (pump wave). We can derive the material, classical elastic nonlinearity by comparing the probe signal with and without the pump wave's presence.

The experimental design we use is shown below. We use an arbitrary waveform generator (Keysight 33500B Series) to excite both the pump and the probe signals. We choose the probe wave to be a low-amplitude (5 V), high-frequency (500 kHz) signal relative to the pump wave's high-amplitude (10 V amplified 50x by a TEGAM HIGH Voltage Amplifier), low-frequency (50 kHz) character. We use a relatively low amplitude for the probe wave to preserve the medium's elastic parameters during the probe's transmission (i.e., the source of perturbation is due exclusively to the pump wave and not the probe). We use a high amplitude for the pump to induce a large enough perturbation to be sensed by the probe wave. We choose the probe and pump frequencies such that the pump wavelength is significantly larger (~ x10) than that of the probe wavelength. This choice results in an approximately steady-state pump wave at the timescale of the probe, thus granting us the ability to sense different phases of the pump waveform with the high-frequency probe.
### Experimental Setup
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Experimental_Setup.jpg" width = 600>
          
          The transient wave dynamic acousto-elastic testing (TW-DAET) experimental setup we employ.
</p>

### Data Collection Protocol
The velocity perturbation and thus the probe travel time delay imposed by the pump that we aim to measure is small, on the order of nanoseconds. To measure such a small signal, we adopt the recording workflow developed again by Gallot et al. (2014). The procedure is as follows:
1. Send and record a probe signal, S1,
2. Send and record a pump signal, S2,
3. Send and record a pump and probe signal simultaneously as S3,
4. Calculate the perturbed probe: S4 = S3 − S2,
5. Calculate the cross-correlation between the probe signal, S1, and the perturbed probe, S4,
6. Calculate a parabolic fit using the four points nearest the maximum of the cross-correlation function.
7. Take the peak of parabolic fit as the time delay between the probe and perturbed probe.

We repeat this process multiple times at different pump transmission delays to evaluate how the elastic nonlinearity changes with varying phases of the pump. In figure
(2.2), we show a schematical view of one of our samples as we progress through three transmission delay values. At 0 μs, the probe interacts with the pump at a minimum, while at 10 μs, the pump wave has progressed such that the probe senses a maximum. Therefore, as we delay transmission, the probe senses multiple transitions between the pump wave’s positive and negative polarity, allowing for a complete characterization of elastic nonlinearity. We layout the steps of the data collection process in figure (2.3).
### Orientation Diagrams
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Orientation_diagrams_Rev3.PNG" width = 600>
          
          Here we show the pump-probe orientation diagrams employed during TW-DAET. Orientation one (left) defines the pump and probe propagation directions parallel to x and y, respectively, with each wave, polarized vice-versa. Orientation two (right) defines propagation directions of the pump and probe in x and z, respectively, again with each wave polarized vice-versa.
</p>

### Cement Samples
<p align="center">          
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/All_cement_samples.png" width = 600>
          
          Schematic depicting the six prismatic inhomogeneous cement samples on which we complete TW-DAET. Samples A, B, and C consist of embedded unconsolidated sand layers of varying thickness and separation. The inhomogeneities in samples D, E, and F consist of two types of metallic cylinders (metal rods; D and copper wires; E and F) again with varying diameter and separation.
</p>
          

### Pump Wave - Sand Layer Interaction
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/Sand_layers_and_pump.png" width = 600>
          
          Schematic depicting how the pump wave (red) interacts with horizontal unconsolidated sand layers in orientation one (left) and orientation two (right).
</p>
          
### Pump Wave - Metal Cylinder Interaction          
<p align="center">
          <img src="https://github.com/jen053/Nonlinearity/blob/master/Images/Set-up/wires_and_pump.png" width = 600>
          
          Schematic depicting how the pump wave (red) interacts with a network of Cu wires in orientation one (left) and orientation two (right).
</p>

### Supporting Documentation
We give an extensive bibliography within the actual thesis itself. Here are a few fundamental articles which layout the project's foundation.

1. Gallot, T., Malcolm, A., Szabo, T., Brown, S., Burns, D., & Fehler, M., 2014. Characterizing the nonlinear interaction of s- and p-waves in a rock sample, Journal of Applied Physics, 117.

2. Guyer, R. & Johnson, P., 1999. Nonlinear mesoscopic elasticity: Evidence for a new class of materials, Physics Today - PHYS TODAY, 52, 30–36.

3. Haupert, S., Riviére, J., Anderson, B., Ohara, Y., Ulrich, T., & Johnson, P., 2014. Optimized dynamic acousto-elasticity applied to fatigue damage and stress corrosion cracking, Journal of Nondestructive Evaluation, 33, 1–13.

4. Hayes, L. O., Malcolm, A., Moravej, K., & Butt, S. D., 2018. Nonlinear interactions of p and s waves under uniaxial stress, Proceedings of Meetings on Acoustics, 34(1), 045012.

5. Khajehpour Tadavani, S., Poduska, K., Malcolm, A., & Melnikov, A., 2020. A nonlinear elastic approach to study the effect of ambient humidity on sandstone, Journal of Applied Physics, 128, 244902.

6. Landau, L. & Lifshitz, E., 1986. Theory of elasticity, Course Theor. Phys., 7.

7. Riviére, J., Renaud, G., Guyer, R., & Johnson, P., 2013. Pump and probe waves in dynamic acousto-elasticity: Comprehensive description and comparison with nonlinear elastic theories, Journal of Applied Physics, 114.

8. Rusmanugroho, H., Malcolm, A., & Darijani, M., 2019. A numerical model for the nonlinear interaction of elastic waves with cracks, Wave Motion, 92, 102444.

9. TenCate, J., Malcolm, A., Feng, X., & Fehler, M., 2016. The effect of crack orientation on the nonlinear interaction of a p wave with an s wave: Modulation of p waves by s waves, Geophysical Research Letters, 43.
