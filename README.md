# p-WIA for virtual population
 Code for Imperial college London Msc Research Project by Xiaoqi Lin: 

Comparison of arterial Wave Intensity Analysis by pressure only and by pressure-velocity methods using computational modelling

Code was modified from # Sphygmocor-reservoir
Reservoir Analysis using a matlab batch process for multiple Sphygmocor© files (bRes_sp)

For background see *Alun Hughes, Kim Parker. The modified arterial reservoir: an update with consideration of asymptotic pressure (P∞) and zero-flow pressure (Pzf). Proceedings of the Institution of Mechanical Engineers, Part H: Journal of Engineering in Medicine in press.* https://doi.org/10.1177/0954411920917557 and Kim Parker's website pages on Reservoir/excess pressure (http://www.bg.ic.ac.uk/research/k.parker/res_press_web/rp_home.html).

  Main_v06.m  - Code for batch analysis of p-WIA on healthy virtual population of n=240.
  Main_v06_single.m  - Code for single subject analysis of p-WIA on healthy virtual population.
  Main_v06_Svi.m   - Code for analysis of p-WIA simulated impaired stroke volume patients of n=2000.
  Main_v06_tp.m - Code for applying p-WIA on ASCOT dataset from the CAFE Study.

## Using the script

Put dataset on directory:

C:\\Spdata

Open matlab script and ensure that the working directory is the one that contains the relevant script and function files

Open matlab and ensure that the working directory is the one that contains the relevant script and function files.

Run the scrpit and after some time (depending on how many files are analysed the run should complete, returning to command prompt. It will show progress on command window while it is running.

A new folder should now exist in  C:\\Spdata\\results.


