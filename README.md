Cartesian package consists of three different modules: cartesian.py, cartesian_ana.py, cartesian_batch.py

###############################################################################################
Cartesian.py: parses coordinates from GROMACS .xvg (gmx traj) and .pdb files of MD trajectories.
  options:
  -h, --help     show this help message and exit
  --f F          Input traj (.xvg from GROMACS gmx traj or .pdb of trajectory) (aligned for rotation+translation or for rotation in
                 conjuction with --s COM file)
  --s S          OPTIONAL Input traj -com -ox filepath
  --o O          Output.cart (parquet format)
  --pdbnp PDBNP  OPTIONAL number of asynchronous processes for .pdb parsing. =1 turns on serial backend, >1 parallel. =8 is
                 recommended

    > outputs parquet-type files <name>.cart
    > can be used to recalculate vectors using an additional .xvg (gmx traj -com) file, so all the vectors start at the
      COM position
    > these output vectors are then used for either pair-wise analysis (cartesian_ana.py) or batch analysis (cartesian_batch.py)
    > .pdb parsing is slower than .xvg, so using an .xvg format is recommended. To speed up reading of a .pdb, use --pdbnp.
###############################################################################################
Cartesian_ana.py: uses .cart vectors to do a pair-wise analysis of two MD trajectories.
  options:
  -h, --help            show this help message and exit
  --f F                 Input vector.json for trajectory 1
  --s S                 OPTIONAL Input vector.json for trajectory 2
  --o O                 Output directory, defaults to names of trajectories used separated by _
  --pdbs PDBS           OPTIONAL Input structure in .pdb format of the second trajectory
  --method {grid,grid_scan,violin,volume,all}
                        Method - grid, grid_scan, violin, volume, all; defaults to all
  --grid GRID           Grid size in nm, defaults to 0.1 nm
  --gridbackend {0,1}   Backend to use for unique grid assignment
  --p P                 Disable all visualization and plotting, overrides all other plot parameters
  --pdbshow PDBSHOW     Creates a PyMol visualization of dynamicity and position change, requires --pdbs. Default=True
  --plotrscore [PLOTRSCORE]
                        Plot R-scores from grid, defaults to True
  --plot_3d [PLOT_3D]   OPTIONAL plot results in a 3D-plot
  --plot_diff [PLOT_DIFF]
                        OPTIONAL plot explored volume by atom
  --plot_violin [PLOT_VIOLIN]
                        OPTIONAL plot violin plots of spatial positions

   @@@ Volume method
   For each atom's positions XYZ as a function of time in the trajectory, all the points are connected into a polyhedron using
   Convex Hull method, which is then used to approximate EXPLORED VOLUMES through the whole trajectory for the atom, and it's
   per-step value (independent of trajectory length or the amount of datapoints).
   These values are output in a .csv format and plotted (--plot_diff) to compare the two trajectories. If a specific atom shows higher
   value of an explored volume than the other per-step, it means that it explored a larger volume overall (i.e. information about its dynamicity)

   @@@ Grid
   This method identifies regions of conformational change (change in distributions positions). It uses a custom scoring functions, which is
   independent on the distributions widths - it's a score valid for differentiating two multidimensional distributions and how much they differ.

   First, according to a parameter --grid, all the points are binned into a 3D histogram and only the unique bins (and their populations) are output into
   a .cart format. This substantially decreases the amount of datapoints for further calculations without much loss in accuracy (and given the approximate nature of MD and forcefields,
   we don't get accuracy on the order of <0.1 nm anyway). Each of these bins is then used to calculate the R-score. This is then output in a normalized (-1, 1) form into a .csv file.

   This R-score is formulated by (i, j are bins of trajectory 1 and trajectory 2, we're summing over all of them):
   R-pair(atomN) = sqrt( (1/samples) * SUM(i) * SUM(j) [ {deltaX_ij}^2 + {deltaY_ij}^2 + {deltaZ_ij}^2 ] * [Population_of_bin_i * Population_of_bin_j]   )
   R-int1(atomN) = R-pair only for trajectory 1
   R-int2(atomN) = R-pair only for trajectory 2

   >> R-score(atomN) = 2x(R-pair(atomN)) - (R-int1(atomN) + R-int2(atomN)) <<

   The higher the value of the score, the more the two distributions differ for the N-th atom.

   @@@ Grid_scan
   In order to find the optimal --grid parameter, a scan is performed with values of 3, 2, 1.5, 1, 0.5, 0.25, 0.10, 0.05 nm.
   R-scores are output for all the grid values, the one which is the largest but "converges" with the smaller ones is usually the optimal one.

   @@@ Violin
   Plots violin plots of XYZ distributions for either one or both of the trajectories, this doesn't rely on binning and is useful for visualizing distribution differences
   between atoms and trajectories in all axes.

   @@@ All
   Combines grid, volume and violin and produces a comprehensive plot.
###############################################################################################
