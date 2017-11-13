The q-factor-based Digital Image Correlation Algorithm.

## Purpose
The repository contains the MATLAB m-files for qDIC along with synthetic example images. The qDIC algorithm determines 2D displacement fields between consecutive images or from a static reference image to a current image. 

## Running qDIC

### Software Requirement
MATLAB 2011b (for "griddedInterpolant") and the associated Image Processing Toolbox (for other miscellaneous function calls) are the minimum supported requirements to run this code.  Under some circimstances older versions may function using "interpn", but performance and/or accuracy may suffer (and you'll have to implement the change to "interpn").  Development is currently under Matlab 2017a on CentOS 7 and Window 7 x64.

### Input Image Requirements
* To check if the images have the required speckle pattern and intensity values for correlation please use our [DIC simulator](https://github.com/FranckLab/DIC-Simulator).
* We recommend that the input image stack  should have at least 3 times the subset size as the number of pixels in each dimension. The default subset size is 64x64, meaning the the minimum input image size should be 192x192.
* Non-square images are acceptable
* The fundamental image type used for input is .mat
* Out-of-the-box qDIC supports most common Matlab-readable images with `img2mat.m`, other file formats require simple modification

### Running including example case
1. Make sure that the main files are in the current (working) directory for MATLAB. 
2. Copy the desired test images `test_images` directory as needed.
3. Run the `exampleRunFile.m` file to get 2D displacement fields between the two images. Note that the displacement output is in an nx1x3 cell array.

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

In many cases where the example images fail to run, the minium specifications for MATLAB have not been met.

## Cite
If used please cite:
[Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast iterative digital volume correlation algorithm for large deformations. Experimental Mechanics. doi: 10.1007/s11340-014-9874-2](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst)

```bibtex
@article{bar2014fast,
  title={A fast iterative digital volume correlation algorithm for large deformations},
  author={Bar-Kochba, E and Toyjanova, J and Andrews, E and Kim, K-S and Franck, C},
  journal={Experimental Mechanics},
  pages={1--14},
  year={2014},
  publisher={Springer}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/FIDIC#faq) and [Questions/Issues](https://github.com/FranckLab/FIDIC/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).
