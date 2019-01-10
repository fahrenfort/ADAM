CubeHelix Colormaps
===================

This is the only colorscheme-generator you will ever need for MATLAB: Cubehelix are simply the best colormaps for published or distributed documents as they are very attractive in full color and yet are suitable for grayscale conversion. This submssion consists of an extremely versatile colormap generator and visualization tool.

This submission allows you to create different colormaps using just a few parameters. The standard Cubehelix algorithm offer very attractive colorschemes for online and electronic documents (e.g. PDF), and yet when printed in grayscale they keep exactly the sequence information of the original data. This submission also includes two extra controls over the range and domain of the Cubehelix scheme, giving a practically unlimited number of colormaps with many different styles: maximally distinct, multi or single hue, suitable for grayscale printing or even simple grayscale.

This submission includes three functions for working with Cubehelix colormaps:

* CUBEHELIX returns a colormap created using Dave Green's Cubehelix colorscheme function. 
* CUBEHELIX_VIEW creates a figure for creating Cubehelix colorschemes with real-time interactive adjustment of the scheme's parameter values, plus a 'random' demonstration mode and the ability to control other figures' or axes' colormaps.
* CUBEHELIX_FIND can be used to retrieve the parameters from an existing Cubehelix colormap, or to find the best Cubehelix colorscheme that matches any selection of colors (e.g. a document or corporate colorscheme).

### Cubehelix ###

Cubehelix colorschemes consist of nodes along a tapered helix in the RGB color cube, with a continuous increase in perceived intensity (e.g. black->white). Thus the scheme defines attractive colormaps with a huge choice of hue, saturation and brightness, and yet printing a figure (or image) in Black-and-White (e.g. postscript) results in a monotonically increasing grayscale that retains the brightness order of the original colormap. The sequence information of the colormap is retained even in grayscale, which means an attractive colored image can be printed in grayscale and still be informative to the end-user.

The scheme is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf

For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/

Note: The original specification (the links above) misnamed the saturation option as "hue". In this submission the saturation option is named "sat".

### Examples ###

	% New colors for the "colormap" example: 
	load spine 
	image(X) 
	colormap(cubehelix)
	
	% New colors for the "surf" example: 
	[X,Y,Z] = peaks(30); 
	surfc(X,Y,Z) 
	colormap(cubehelix([],0.5,-1.5,1,1,[0.29,0.92])) 
	axis([-3,3,-3,3,-10,5])

### Examples of Viewing Cubehelix Colormaps ###

	% Interactive colorscheme parameter viewer: 
	cubehelix_view 
	% Set/reset the viewer with new parameter values: 
	cubehelix_view([],0.5,-1.5,1,1)
	
	# Control external axes/figure colormaps:
	load spine
	image(X)
	cubehelix_view({gca})

### Examples of Retrieving Colormap Parameters ###

	cubehelix_find(cubehelix(10))
	 ans = [0.5,-1.5,1,1]
	
	map = cubehelix(10, 1.4,-0.7,0.9,1.2, [0.05,0.97]); 
	[vec,irg,dmn] = cubehelix_find(map)
	 vec = [1.4,-0.7,0.9,1.2]
	 irg = [0.05,0.97]
	 dmn = [0,1]

	map = cubehelix(64, [2.3,0.4,0.5,0.6], [0.05,0.24], [0.19,0.85]);
	[vec,irg,dmn] = cubehelix_find(map)
	 vec = [2.3,0.4,0.5,0.6]
	 irg = [0.05,0.24]
	 dmn = [0.19,0.85]
