# Multilayer Structure
A Multilayer structure or stratified media is a stack of parallel films.

### Theory
The solution is based on calculating the reflected and transmitted plane waves at each interface and the phase change due to light propagation inside each layer. The solution is a geometric sum which considers all the internal round trips generated inside the multilayer structure as illustrated in the figure above.

---

### Approximations
None. The implementation solves Maxwell equations for any type of multilayer structure and any incident beam.

---

### Constructor
MultilayerStructure(_n__A_, _h__A_, _ref_)
- **n_A** - vector specifying the refractive indexes of all layers;
- **h_A** - vector specifying the thickness of all layers. The layer with refractive index _n__A[i]_ has a thickness of _h__A[i-1]_;
- **ref** - ReferenceFrame specifying the position and orientation of the first interface;

---

### Implementation

**lightinteraction - FieldAngularSpectrumScalar**
- **Time:** very short; scales with Nx Ny
- **RAM:** None
- **Convergence** sampling of nsx and nsy

**Scattering matrix - FieldAngularSpectrumScalar**
- **Type:** Transmission and reflection matrices are diagonal
- **Time:** very short; scales with Nx Ny
- **RAM:** very small; scales with Nx Ny
- **Convergence** sampling of nsx and nsy

---

### Applications
Example of systems studied with similar models:
- [Fabry-Pérot etalon](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957)
- [Aberations induced by a coverslip](https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-23/issue-09/090801/Tool-for-simulating-the-focusing-of-arbitrary-vector-beams-in/10.1117/1.JBO.23.9.090801.full?SSO=1)

---

### Reference
If you use this model, please cite:
[Marques, Dylan M., et al. "Modelling fabry-pérot etalons illuminated by focussed beams." Optics express 28.5 (2020): 7691-7706.](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-5-7691&id=427957)

---

### Additional reading:
- Chapter 2 of Fujiwara, Hiroyuki. Spectroscopic ellipsometry: principles and applications. John Wiley & Sons, 2007.
- [Van de Nes, A. S., et al. "Calculation of the vectorial field distribution in a stratified focal region of a high numerical aperture imaging system." Optics Express 12.7 (2004): 1281-1293.](https://www.osapublishing.org/oe/abstract.cfm?uri=OE-12-7-1281)
