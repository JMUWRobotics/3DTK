# 3DTK - The 3D Toolkit

[![Build Status](https://travis-ci.org/3DTK/3DTK.svg?branch=master)](https://travis-ci.org/3DTK/3DTK)
[![Build status](https://ci.appveyor.com/api/projects/status/7f8n86ed859vw8j4/branch/master?svg=true)](https://ci.appveyor.com/project/josch/3dtk/branch/master)

## About
_The 3D Toolkit_ provides algorithms and methods to process 3D point clouds. It includes automatic high-accurate registration (6D simultaneous localization and mapping, 6D SLAM) and other tools, e.g., a fast 3D viewer, plane extraction software, etc. Several file formats for the point clouds are natively supported, new formats can be implemented easily. Check out our [website (threedtk.de)](http://threedtk.de) for more information.

## Build Instructions
For installation instructions on **Ubuntu**, **Debian** and **macOS**, please have a look at the [**INSTALL file**](INSTALL.md).

To compile the project simply call `make`. This will configure slam6d using the default settings. If you wish to configure the project using custom settings do: `make config`. This command requires _ccmake_ be installed on your system.  Alternatively you may change into the build directory `.build` and configure the project with your preferred cmake configurator, i.e.:
```
mkdir .build
cd .build && cmake ..
```
For **Microsoft Windows**, use the _cmake-gui_ application provided by _cmake_ to configure and generate project files for the appropriate version of _Microsoft Visual Studio C++_ of your system. Use the `INSTALL`  target to built the entire project.  Executables (and .dll's) will then reside in the `windows` folder. For running the binaries you need to install the proper redistributable package.
Alternatively on **Windows**, execute the batch script `build.cmd` either from the command line or from the Explorer. It will download and extract the right versions of _boost_ and _opencv_ and then run _cmake_ with the right configuration options.

Some Boost libraries (_graph, regex, serialization, filesystem, interprocess_) are needed to compile the slam6D program. **_3DTK_ is not compatible with boost version 1.49.** Please make sure to use a different version.

### 3DConnexion Space Mouse Support
The Applications `show` and `qtshow` supports navigation with the _3DConnexion SpaceMouse_. The left key on the _SpaceMouse Wireless/Compact_ fixes the rotation and the right key fixes the translation.
The keyset for other models can be different.

#### Linux
No further steps are required for the basic functions. For full support the `libspnav` and the `spacenavd` daemon ([show on GitHub](https://github.com/FreeSpacenav)) is required.
Then build _3DTK_ with the `WITH_3DMOUSE` option.

#### macOS
It's required to install the official driver and build _3DTK_ with the `WITH_3DMOUSE` option.

<!--#### Windows
Should be the same as for macOS, but is not tested!-->

## Getting Started
For a detailed explanation of the programm, its usage, etc., please refer to the [high level documentation](doc/documentation_HL.pdf) (esp. sections 4-6, for starters) and [here](doc/3d_video.md) for creating videos with _3DTK_. Further a reference manual can be found in `doc/html/index.html` after building the doxygen domcumentation (type in `make docu` to get the HTML files).

##### Important:
Take care to register scans first (`bin/slam6D`) before trying to display them (`bin/show`), and think about using the point reduction ([see section 6](doc/documentation_HL.pdf)) for a much faster scan matching result. Extremely large scans might need to be reduced (using `bin/scan_red`) before registration. This will write reduced scans in the `uos` format into a directory `reduced` in the data directory.

### Datasets:
Three example scans are included in the dat directory, several larger data sets can be downloaded from the [data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/).

### Examples:
#### slam6d and show
###### Using the data set in the slam6d repository:
```
bin/slam6D -m 500 -R 5 -d 25.0 --metascan dat
bin/show dat
```
```
bin/slam6D --max=500 -r 10.2 -i 20 --metascan dat
bin/show dat
```

###### Using `hannover1.tgz` [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
```
bin/slam6D -s 1 -e 65 -r 10 -i 100 -d 75 -D 250 --epsICP=0.00001
           -I 50 --cldist=750 -L 0 -G 1 /YOUR/DOWNLOAD/DIRECTORY/dat_hannover1
bin/show -s 1 -e 65 /YOUR/DOWNLOAD/DIRECTORY/dat_hannover1
```

###### Using `hannover2.tgz` [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
```
bin/slam6D -q -r 10 -f rts -s 23 -d 75 -L 4 --cldist=1500 -G 1 -D -1
           --DlastSLAM 250 --graphDist 200 -I 50 /YOUR/DOWNLOAD/DIRECTORY/hannover2
bin/show -f rts -s 23 /YOUR/DOWNLOAD/DIRECTORY/hannover2
```

###### Using `kvarntorp_mine.tgz` (`dat_mine1`) [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
```
bin/slam6D -s 1 -e 76 -r 10 -m 3000 -d 50 -i 1000 --epsICP=0.000001
           -I 50 -D 75 --clpairs=5000 -f old /YOUR/DOWNLOAD/DIRECTORY/dat_mine1/
bin/show -s 1 -e 76 -m 3000 -f old /YOUR/DOWNLOAD/DIRECTORY/dat_mine1/
```

###### Using `bremen_city.zip` [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
```
bin/scan_red -s 0 -e 12 -r OCTREE -v 10 --octree 0 -f RIEGL_TXT /YOUR/DOWNLOAD/DIRECTORY/bremen_city
bin/slam6D -a 2 -q /YOUR/DOWNLOAD/DIRECTORY/bremen_city/reduced -f uos -d 150
            -s 0 -e 12 --anim=1 -n /YOUR/DOWNLOAD/DIRECTORY/bremen.net
		  -G 1 -D 100 -i 0 -I 50 -p --epsSLAM=0.0
bin/show -s 0 -e 12 /YOUR/DOWNLOAD/DIRECTORY/bremen_city/reduced
```

###### Using `UniKoblenz_CampusTour3_OsnabrueckFormat.tar.gz` [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
```
bin/slam6D -s 1 -e 320 -r 20 -i 300 --epsICP=0.000001 -d 45 -D 45
           -f uos --algo=2 -l 10 -L 4 -I 100 -G 1
           /YOUR/DOWNLOAD/DIRECTORY/UniKoblenz_CampusTour3_OsnabrueckFormat/
bin/show -s 1 -e 320 -f uos /YOUR/DOWNLOAD/DIRECTORY/UniKoblenz_CampusTour3_OsnabrueckFormat/
```

#### Detecting Planes
For detecting planes compile with the `WITH_PLANE` option. Adapt the settings in `bin/hough.cfg` for your data set.
Using the data set in the slam6d repository, no modification of `bin/hough.cfg` necessary:
```
bin/planes -s 0 dat
bin/show -s 0 -e 0 dat -l dat/planes/planes.list
```

###### Using `bremen_city.zip` [from data repository](http://kos.informatik.uni-osnabrueck.de/3Dscans/):
Adapt these settings in `bin/hough.cfg`:
```
RhoNum              500
RhoMax              5000
MaxPointPlaneDist   50.0
MinPlaneSize        50
PointDist           100.0
```
```
bin/planes -f riegl_txt -s 0 /YOUR/DOWNLOAD/DIRECTORY/bremen_city/ -r 50 -O 1 -m 5000
bin/show -s 0 -e 0 /YOUR/DOWNLOAD/DIRECTORY/bremen_city/ -f riegl_txt -l dat/planes/planes.list -r 10 -O 1 -m 5000
```

## Frequent Use

If you find yourself using 3DTK frequently, read about [Configuration Files](doc/Configuration%20Files.md) to find out how not to have to type as much on the command line.
