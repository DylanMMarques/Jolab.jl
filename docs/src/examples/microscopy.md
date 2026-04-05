# Microscopy

A example application of optics is microscopy. The aim of microscopy is to image samples with high resolution. Different types of microscopy exist, each with its own advantages and disavantages. This section exemplifies how to use Jolab to simulate different types of microscopy. 

:::tabs

== tab Widefield microscopy

Widefield microscopy is a common type of microscopy where the sample is illuminated with a wide beam. The image is then formed by the light transmitted or reflected by the sample. The simplest widefield microscopy is composed of a sample illuminated by a collimated beam and two lenses that form the image of the sample of a detector. 

```@example widefield
using Jolab, CairoMakie


```

== Dark-field microscopy

The principle of dark-field microscopy is to illuminate the sample with a beam that does not directly reach the detector. This enables to only collect light that is scattered by the sample, which can be low intensity. This type of microscopy is useful to image samples that are transparent or have low contrast in bright-field microscopy. The simplest dark-field microscopy is composed of a sample illuminated by an annular beam and two lenses that form the image of the sample of a detector.

== Phase-contrast microscopy

Phase-contrast microscopy is a type of microscopy that enables imaging of transparent samples. The principle is to convert phase variations of the light transmitted through the sample into intensity variations. This is achieved by placing a phase plate in the back focal plane of the objective lens. The phase plate shifts the phase of the light that is not scattered by the sample, while leaving the phase of the scattered light unchanged. This results in an interference pattern that can be detected by a camera.

== Differential interference contrast microscopy

Hello something

== Confocal microscopy

In confocal microscopy, a sample is illuminated with a focused laser beam. To recover the image, the focussed spot is scanned across the sample. Typically, light emmited by the sample (e.g. fluorescence) is typically collected by the same objective and focused on a pinhole. Light passing through the pinhole is detected by a photodetector. In these examples, light detected is the light reflected by the sample.


:::
