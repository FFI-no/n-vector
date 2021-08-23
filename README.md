# *n*-vector library in Matlab

Author: Kenneth Gade (PhD, principal scientist at FFI)

This library is for performing a range of different geographical position calculations, e.g. adding, subtracting, interpolating, and averaging positions. Exact results are returned for either ellipsoidal or spherical Earth model.

All calculations are based on [*n*-vector](https://www.navlab.net/nvector) (which can replace e.g. latitude and longitude), giving several advantages:

* The calculations are typically **simple and intuitive**
  * Reason: *n*-vector is a 3D vector and hence the powerful vector algebra can be used to solve many position calculations intuitively and with few code lines.
* The calculations are **non-singular** (i.e. work equally well at or near the poles as any other global position)
  * Reason: The *n*-vector representation is inherently non-singular for all Earth positions.
* The calculations have **no discontinuities** (i.e. work equally well across the dateline (±180° longitude meridian) as any other global positions)
  * Reason: The *n*-vector representation has no discontinuities.
   

For more details and 10 examples of usage, see https://www.navlab.net/nvector 

**Reference:** 
> Kenneth Gade (2010): A Non-singular Horizontal Position Representation, *The Journal of Navigation*, Volume 63, Issue 03, pp 395-417, July 2010, [DOI: 10.1017/S0373463309990415](https://doi.org/10.1017/S0373463309990415). <br/> 
*<https://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>*

### Also available in several other programming languages:

There are several *n*-vector libraries (from other authors) available in other programming languages that are either based on this library, or based directly on the reference article ([Gade, 2010](https://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)). Some of these are listed below. Note: the listed code is not verified by the author of this library.

* **C# (C Sharp)**: https://www.navlab.net/nvector/#download
* **C++**: https://www.navlab.net/nvector/#download
* **Python**: https://github.com/pbrod/Nvector
  * Alternative 2: https://github.com/mrJean1/PyGeodesy 
  * Alternative 3: https://github.com/lxnt/ccnvector
* **JavaScript**: https://github.com/chrisveness/geodesy
* **Haskell**: https://github.com/ofmooseandmen/jord
* **Go / Google's Go / Golang**: https://github.com/fortyninemaps/nvector
*	**R**: https://github.com/euctrl-pru/nvctr 
*	**Java**: https://github.com/omeruluoglu/JGeodesy 
*	**Rust**: https://github.com/ofmooseandmen/jord-rs

## Contributing
If you have suggestions or other comments to the code, please contact the author at kenneth.gade@ffi.no. 

## Using this code
The code works with all Matlab versions (and no toolboxes are needed). A simple example of usage is Example 7 from https://www.navlab.net/nvector, where a horizontal midpoint is calculated:

```matlab
% Three positions A, B and C are given as lat/long in degrees.
% Convert all three to radians and then to n-vectors:
n_EA_E = lat_long2n_E(rad(90),rad(0));
n_EB_E = lat_long2n_E(rad(60),rad(10));
n_EC_E = lat_long2n_E(rad(50),rad(-20));
 
% Find the horizontal mean position, M:
n_EM_E = unit(n_EA_E+n_EB_E+n_EC_E);
```

As described below, the file `examples.m` includes the above solution, and it also gives the following plot:

![Earth figure from Example 7](/Example7.png?raw=true)


## Matlab files included

Below is a list of the files that are included in this library (19 files in total). The file name syntax, mathematical symbols and coordinate frames are defined at https://www.navlab.net/nvector.

Convert between lat/long and *n*-vector:

* `lat_long2n_E.m` Converts latitude and longitude to *n*-vector
* `n_E2lat_long.m` Converts *n*-vector to latitude and longitude

Convert between delta (i.e. local position vector) and  *n*-vectors:

* `n_EA_E_and_n_EB_E2p_AB_E.m` From two positions *A* and *B*, finds the delta position
* `n_EA_E_and_p_AB_E2n_EB_E.m` From position *A* and delta, finds position *B*

Convert between *n*-vector and ECEF-vector (i.e. position vector from Earth center, in meters):

* `n_EB_E2p_EB_E.m` Converts *n*-vector to ECEF-vector
* `p_EB_E2n_EB_E.m` Converts ECEF-vector to *n*-vector

Convert between *n*-vector and rotation matrix (i.e. ***R***<sub>*EN*</sub> or ***R***<sub>*EL*</sub>):

* `R_EN2n_E.m` Finds *n*-vector from  ***R***<sub>*EN*</sub>
* `n_E2R_EN.m` Finds ***R***<sub>*EN*</sub> from *n*-vector
* `R_EL2n_E.m` Finds *n*-vector from ***R***<sub>*EL*</sub>
* `n_E_and_wa2R_EL.m` Finds ***R***<sub>*EL*</sub> from *n*-vector and wander azimuth angle

Convert between Euler angles and rotation matrix:

* `xyz2R.m` Creates a rotation matrix from 3 angles about new axes in the xyz order
* `R2xyz.m` Three angles about new axes in the xyz order are found from a rotation matrix
* `zyx2R.m` Creates a rotation matrix from 3 angles about new axes in the zyx order (e.g. yaw-pitch-roll)
* `R2zyx.m` Three angles about new axes in the zyx order (e.g. yaw-pitch-roll) are found from a rotation matrix

Miscellaneous simple utilities:

* `unit.m` Makes input vector unit length (i.e. norm = 1)
* `rad.m` Converts angle from degrees to radians
* `deg.m` Converts angle from radians to degrees
* `R_Ee.m` Selects axes of the coordinate frame *E*

Examples of how to use the *n*-vector library for position calculations:

* `examples.m` Contains solutions to the 10 examples given at https://www.navlab.net/nvector

A note about the axes of coordinate frame *E*: This library uses the most common *E*-frame by default (i.e. the *z*-axis towards the North Pole). However, if a less common (but in many cases better) choice is preferred (with the *x*-axis towards the North Pole), this can be switched in the `R_Ee.m`-file. More details are available in the help text of this file and in Table 2 of [Gade (2010)](https://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf).

## License
This repository is available under the MIT License. See the license file for details.

## Citing this code
The *n*-vector library is based on the following article, and hence it should be cited in publications using the library:

> Kenneth Gade (2010): A Non-singular Horizontal Position Representation, *The Journal of Navigation*, Volume 63, Issue 03, pp 395-417, July 2010, [DOI: 10.1017/S0373463309990415](https://doi.org/10.1017/S0373463309990415).<br/> 
*<https://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>*

Bibtex entry as follows:
```
@article{gade2010non,
  title={A Non-singular Horizontal Position Representation},
  author={Gade, Kenneth},
  journal={The Journal of Navigation},
  volume={63},
  number={3},
  pages={395--417},
  year={2010},
  url={https://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf}, 
  doi={10.1017/S0373463309990415},
  publisher={Cambridge University Press}
}
```
