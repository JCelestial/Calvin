# CALVIN v.0.4.8
### *For when you are tired of writing new scripts for every Molecular Simulation....*
CALVIN will allow you to go from simulation to plot-ready data in minutes!

<img src="https://github.com/JCelestial/Calvin/blob/master/tio2_il_system.png" width="275" height="300">  <img src="https://github.com/JCelestial/Calvin/blob/master/old_algo_ordpar.gif" width="350" height="300">

And more!

------------------------------------------------------------------------------------------------------------------------------------
### What is CALVIN?
Calvin stands for **C**omputes **A**bsurdly **L**arge, **V**apid, **I**ndexed **N**umbers and is a utility that is used to process large scale trajectory data from Molecular Dynamics simulation programs such as [VMD](https://www.ks.uiuc.edu/Research/vmd/). In programs like VMD, simulation trajectory coordinates can be exported under multiple formats, e.g. *.xyz*, *.pdb*, or *.dcd*; as of this current version, CALVIN is most compatible with *.xyz* and currently looking to build classes to expand upon CALVIN's usability for more extensions.

------------------------------------------------------------------------------------------------------------------------------------
#### Current methods available:
* Print coordinates (no header lines, just pure trajectories)
* Center of Mass
* Angular Conformation
* [Order parameter](https://en.wikipedia.org/wiki/Phase_transition)
* Energy Average and standard deviation

#### Planned methods:
* [Radial distribution](https://en.wikipedia.org/wiki/Radial_distribution_function)
* Protein surface area

------------------------------------------------------------------------------------------------------------------------------------
### Why do we need CALVIN?
Simulation coordinate files are very cumbersome to parse as it is littered with unecessary headers and columns, which is why this version begins with one of the least complex extensions, *.xyz*. A typical *.xyz* file can appear as such:
```
36288
 generated by VMD
  C1        31.869669      -20.711391      -34.581379
  N1        30.844368      -20.086567      -33.779266
  H1        31.610401      -21.606871      -35.055351
  H2        32.800102      -20.909853      -34.014919
  H3        32.144650      -20.023277      -35.347126
  C2        30.857540      -18.913921      -33.099659
  C3        29.650688      -18.775717      -32.479210
  N2        28.933077      -19.915432      -32.859100
  C4        29.690832      -20.693037      -33.642181
```
Not to mention the first two lines shows up intermittently to mark the beginning of every simulation frame. In short, a large simulation containing hundreds of thousands of frames can have hundreds of these headers, which can cause potential parsing errors. In addition, making shell scripts that utilizes **grep** commands and regex leads to messy outcomes. So, how does CALVIN deal with this?

CALVIN takes advantage of these headers by first appending the *.xyz* files with the *terminate_xyz* script that will attach a termination sequence at the very end of the file and then using those lines as checkpoints for allocation and deallocation of memory space, preventing memory leaks.
```bash
$ ./terminate_xyz trajectories.xyz > terminatedTrajectories.txt
```
Once appended, CALVIN can properly parse these files and due to the nature of how MD programs arrange their data, this allows CALVIN to analyze the simulation metrics of the simulation, containing the following:
* Total number of simulation frames
* Total number of molecules within simulation frames
* Number of atoms per molecule

```
========================================================
CALVIN : The data has the following array dimensions... 

Simulation Frames: 500
Molecules per Frame: 672 
Atoms per Molecule: 56
========================================================
```

------------------------------------------------------------------------------------------------------------------------------------
### Why Fortran?

One reason, for speed. Typically, multiple simulated systems have to be analyzed, and due to their large size, it's difficult to analyze all of them concurrently without using [High-Performance Clusters](https://insidehpc.com/hpc101/intro-to-hpc-whats-a-cluster/), so the next best thing is to quickly analyze them one by one. To put it in perspective, it takes minutes for Interpreted languages like Python or R to analyze files that are roughly half a GB in size, whereas CALVIN takes about 30 seconds to give the user simulation metrics and analyses.

With that said, there are future plans on overhauling CALVIN into a different language, such as C++ to support more Object Orientation and the utilization of more Data Structures in order to perform more complex methods. Another proposed alternative is to strip CALVIN of its main interface and leave its methods alone and turn the Fortran components into a dynamic library that can be invoked in memory.

------------------------------------------------------------------------------------------------------------------------------------
### Bugs and Future Plans

#### Bugs
* Fix the garray subroutine to reset the simulation frame metric

#### Future Plans
* Incorporate methods to analyze protein files such surface area
* Expand methods on Energy module to include other statistical metrics
* Overhaul CALVIN into C++ and use Linked Lists instead of Arrays to remove the need for the termination shell script
* Greater support for multithreading to allow concurrent execution of multiple analyses
* Development of a GUI to make CALVIN more user friendly
