# Interpolation Filters for Antiderivative Antialiasing

This repository includes accompanying code for the DAFx24 paper "Interpolation Filters for Antiderivative Antialiasing". You can access the full text [here](https://media.githubusercontent.com/media/IoSR-Surrey/DAFx24-Proceedings/main/papers/DAFx24_paper_33.pdf).

The AA-IIR method with linear and cubic interpolation is applied for simulation of the hard clipping nonlinearity (`hard_clipper.m`) and the diode clipper circuit (`diode_clipper.m`). 

Included libraries are `interpolation` for cubic interpolation methods, `metrics` for signal-to-noise ratio (SNR) and noise-to-mask ratio (NMR) calculation, `utils` for miscellaneous functions (plotting, signal generation, etc.).
