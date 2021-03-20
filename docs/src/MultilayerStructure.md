# Multilayer Structure
A Multilayer structure or stratified media is a stack of parallel films.

<details>
  <summary><b>Theory</b></summary>
  The solution is based on calculating the reflected and transmitted plane waves at each interface and the phase change due to light propagation inside each layer. The solution is a geometric sum which considers all the internal round trips generated inside the multilayer structure as illustrated in the figure above. <br> <br>
  <b> Approximations </b><br>
  None - The implementation solves Maxwell equations for any type of multilayer structure and any incident beam.
  <br><br>
</details>
<details>
  <summary><b>Constructor</b></summary>
  <p>
  <b>MultilayerStructure(n_A, h_A, ref)</b><br>
  <b>n_A</b> - vector specifying the refractive indexes of all layers; <br>
  <b>h_A</b> - vector specifying the thickness of all layers. The layer with refractive index n_A[i] has a thickness of h_A[i-1];<br>
  <b>ref</b> - ReferenceFrame specifying the position and orientation of the first interface;<br>
</details>
<details>
  <summary><b>Implementation</b></summary>
  <p>
  <b> lightinteraction - FieldAngularSpectrumScalar </b> <br>
  <b>Time:</b> very short; scales with N<sub>x</sub> N<sub>y</sub> <br>
  <b>RAM:</b> None <br>
  <b>Convergence</b> sampling of ns<sub>x</sub> and ns<sub>y</sub><br>
  <br>
  <b>Scattering matrix - FieldAngularSpectrumScalar </b><br>
  <b>Type:</b> Transmission and reflection matrices are diagonal <br>
  <b>Time:</b> very short; scales with N<sub>x</sub> N<sub>y</sub> <br>
  <b>RAM:</b> very small; scales with N<sub>x</sub> N<sub>y</sub> <br>
  <b>Convergence</b> sampling of ns<sub>x</sub> and ns<sub>y</sub>
  <br>
  </p>
</details>

<details>
  <summary><b>Applications</b></summary>
  Example of systems studied with similar models:
  <ul>
      <li><a href = https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957>Fabry-Pérot etalons</a></li>
      <li><a href = https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-23/issue-09/090801/Tool-for-simulating-the-focusing-of-arbitrary-vector-beams-in/10.1117/1.JBO.23.9.090801.full?SSO=1>Aberations induced by a coverslip</a></li>
  <ul>
</details>

<details>
  <summary><b>Reference</b></summary>
  If you use this model, please cite:<br>
  <a href = https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957> Marques, Dylan M., et al. "Modelling fabry-pérot etalons illuminated by focussed beams." Optics express 28.5 (2020): 7691-7706. </a><br><br>
  <b>Additional reading:</b><br>
  Chapter 2 of Fujiwara, Hiroyuki. Spectroscopic ellipsometry: principles and applications. John Wiley & Sons, 2007. <br>
  <a href = https://www.osapublishing.org/oe/abstract.cfm?uri=OE-12-7-1281>Van de Nes, A. S., et al. "Calculation of the vectorial field distribution in a stratified focal region of a high numerical aperture imaging system." Optics Express 12.7 (2004): 1281-1293.</a><br>
</details>
