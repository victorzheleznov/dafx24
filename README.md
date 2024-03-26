# Interpolation Filters for Antiderivative Antialiasing

This repository includes accompanying code for the DAFx24 submission "Interpolation Filters for Antiderivative Antialiasing".

The AA-IIR method with linear and cubic interpolation is applied for simulation of the hard clipping nonlinearity (`hard_clipper.m`) and the diode clipper circuit (`diode_clipper.m`). 

Included libraries are `interpolation` for cubic interpolation methods, `metrics` for signal-to-noise ratio (SNR) and noise-to-mask ratio (NMR) calculation, `utils` for miscellaneous functions (plotting, signal generation, etc.).
