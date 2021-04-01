# Lens


### Theory
The solution is

---

### Approximations
Thin lens approximation

---

### Constructor
```@docs
Lens(::Type, ::Any, ::Any, ::Any)
```

---

### Implementation

```@docs
lightinteraction(::Lens, ::FieldAngularSpectrumScalar)
```
```@docs
Jolab.coefficient_general(::Lens, ::FieldAngularSpectrumScalar)
```
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
- Goodman, Joseph W. Introduction to Fourier optics. Roberts and Company Publishers, 2005.
- [Wolf, Emil. "Electromagnetic diffraction in optical systems-I. An integral representation of the image field." Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences 253.1274 (1959): 349-357.](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.1959.0199)
- [Richards, Bernard, and Emil Wolf. "Electromagnetic diffraction in optical systems, II. Structure of the image field in an aplanatic system." Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences 253.1274 (1959): 358-379.](https://royalsocietypublishing.org/doi/abs/10.1098/rspa.1959.0200)
- Novotny, Lukas, and Bert Hecht. Principles of nano-optics. Cambridge university press, 2012.
