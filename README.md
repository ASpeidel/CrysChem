# CrysChem
MATLAB scripts to acquire Euler angles from surface topography data

Overview: CrysChemMain.m is a MATLAB script using inbuilt and bespoke functions (detailed below) that returns Euler angles from etch topography data (e.g. acquired from an interferometer). This is achieved by extracting characteristic peak data from the topography data and fitting an orthogonal vector set that simulates (100) etch response.

Requirements: The script requires the MATLAB Image Processing Toolbox and Phased Array System Toolbox.

Methodology: To acquire the characteristic peaks CrysChemMain accesses a topography data file, splits this into separate fields of view, and performs a form removal subtraction. The resulting data (surfFOV) is input into topo.m, which plots azimuthal and elevation data from surfFOV to output Euler angles. 

Characteristic peaks are determined by peakfinder.m, which uses an image thresholding approach if there is no low intensity peak (typical of etched Al surfaces vicinal to (100)) or finds the low intensity peak if one exists. This decision is taken where >20% of the surface intensity falls <20° elevation. peakfinder.m uses Imageproc.m to open and close the topography intensity image prior to thresholding. To deal with the situation where one peak is split over the azimuthal origin (peak 1 <30° && peak n >330°), splitpeaks.m is called, which returns one peak with the mean azimuth and elevation of the parent peaks (peak 1 and peak n).

The output characteristic peaks are then compared with the model to output Euler angles. Where one low angle, high intensity peak is determined (>20% of the surface intensity falls <20° elevation), vectormatching.m is called, which fixes the spherical coordinates of the low angle peak as the Z’’-axis about which the vector set is rotated to find the effective in-plane rotation of the crystal. This is the minimised angular difference between the characteristic peaks and the vector set (AngularDiff.m).

Where there is no low angle, high intensity peak, bestfit1.m is called that finds the minimum separation using AngularDiff.m against all rotations using a coarse fitting pass (6° increments) to define an approximate Euler angle set, prior to a fine pass (1° increments).
  
Inputs: Three column surface height data (for example .xyz files).

Outputs: Either an Euler angle set that corresponds to the underlying crystalline orientation for a specified sampling area, or a non-indexed point (NaN), where <2 characteristic peaks can be determined.
