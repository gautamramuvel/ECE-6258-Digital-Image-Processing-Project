This folder contains the source code of adaptive curvelets as described in the following paper:

	Al-Marzouqi, Hasan, and Ghassan AlRegib. "Curvelet transform with learning-based tiling," 
        Signal Processing: Image Communication 53 (2017): 24-39.
 
Package contents:

    * 'adaptive_curvelet.m' and 'adaptive_icurvelet' implement the forward and inverse adaptive curvelet transforms.
	* 'find_optimal_curvelet_max' and 'find_optimal_curvelet_min' find the adaptive curvelet
	   tiling that better represents a given image by either maximizing or minimizing a cost function.
	* 'findoptimalnumscales': Finds the optimal number of scale decompositoins according to the centered impulse
           response criteria given in [1].
	* 'img_denoise_adaptive_curvelet' and 'img_denoise_default_curvelet': Denoise an image using default and 
	   adaptive curvelets. Default curvelets are assumed to use the 'real' option and periodic extensions.
	*  'fdct_wrapping' and 'ifdct_wrapping': Default forward and invese curvelet tranfsorms.
	*  'fdct_warpping_window': The smoothing function is used in default and adaptive curvelets.
	*  'demo': Denoising and partial reconstruction error examples.


Hasan AlMarzouqi, hasalmarzouqi@pi.ac.ae
Ghassan AlRegib, alregib@gatech.edu
6/22/2017



    

