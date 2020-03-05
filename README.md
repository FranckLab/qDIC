## The Franck Lab q-factor-based Digital Image Correlation Algorithm: qDIC.

## Purpose
The repository contains the MATLAB m-files for qDIC along with synthetic example images. The qDIC algorithm determines 2D displacement fields between consecutive images or from a static reference image to a current image, or using a "hybrid" reference updating scheme.

## Running qDIC

### Software Requirement
MATLAB 2011b (for "griddedInterpolant") and the associated Image Processing Toolbox (for other miscellaneous function calls) are the minimum supported requirements to run this code.  Under some circumstances older versions may function using "interpn", but performance and/or accuracy may suffer (and you may have to implement the change to "interpn").  Development is currently under Matlab 2017a (newer versions should be fine) on CentOS 7 and Windows 7/10 x64.

A "basic" version is available (in beta) that supports base Matlab (i.e., with no Toolboxs) with similar performance from https://github.com/ALandauer/qDICb.  More up-to-date versions may be found at https://github.com/ALandauer/qDIC, although we occasionally reconcile the two versions 

### Input Image Requirements
* To check if the images have the required speckle pattern and intensity values for correlation please use our [DIC simulator](https://github.com/FranckLab/DIC-Simulator).
* We recommend that the input images should have at least 3 times the subset size as the number of pixels in each dimension. The default subset size is 64x64, meaning the minimum input image size should be 192x192.
* Non-square images are acceptable
* The fundamental image type used for input is .tif (or .mat)
* Out-of-the-box qDIC supports most common Matlab-readable images with our `img2mat.m` that calls `imread`, other file formats require simple modification

### Running including example case
1. Make sure that the main files are in the current (working) directory for MATLAB.
2. Copy the desired test images `test_images` directory as needed.
3. Run the `exampleRunFile.m` file to get 2D displacement fields between the two images. Note that the displacement output is in an nx3 cell array.

## Files
* Function files
 - addDisplacements_2D.m
 - checkConvergenceSSD_2D.m
 - DIC.m
 - filterDisplacements_2D.m
 - flagOutliers_2D.m
 - funIDIC.m
 - IDIC.m
 - inc2cum.m
 - FIDICinc2cum.m
 - removeOutliers_2D.m
 - replaceOutliers_2D.m
 - areaMapping_2D.m

* Supplemental .m files from the MATLAB file exchange:
 - inpaint_nans.m
 - mirt2D_mexinterp.m  (Optional, not currently in use)

* Example files to run basic qDIC
 - exampleRunFile.m
 - img2mat.m
 - imageCropping.m
 - image_eval.m
 - Example test images

## FAQ

**What are the requirements for the input images?**

Please refer to [input image requirement](https://github.com/FranckLab/FIDIC#input-image-requirements).

**Can I use qDIC for finding displacement fields in 3D images?**

No. But you can use [FIDVC](https://github.com/FranckLab/FIDVC), this finds 3D displacements in 3D image stack (i.e. a volumetric image). We do not support any 3D-DIC (stereo) functionality.

**Why does the example fail to run?**

In many cases where the example images fail to run, the minimum specifications for MATLAB have not been met.

## Cite
If used please cite:
[Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018). https://doi.org/10.1007/s11340-018-0377-4)

```bibtex
@Article{Landauer2018,
author="Landauer, AK
and Patel, M
and Henann, DL
and Franck, C",
title="A q-Factor-Based Digital Image Correlation Algorithm (qDIC) for Resolving Finite Deformations with Degenerate Speckle Patterns",
journal="Experimental Mechanics",
year="2018",
issn="1741-2765",
doi="10.1007/s11340-018-0377-4",
url="https://doi.org/10.1007/s11340-018-0377-4"
}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/FIDIC#faq) and [Questions/Issues](https://github.com/FranckLab/FIDIC/issues). Add a new question if similar issue hasn't been reported as a GitHub "Issue". The author's contact information can be found at [Franck Lab](francklabbackup.me.wisc.edu) or via GitHub.
