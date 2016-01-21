#FAST v1.0: Fast Automated Spud Trekker Manual#

Please cite [***Aksel T, Yu EC, Sutton S, Ruppel KM, Spudich JA. Cell Reports. 2015. Ensemble force changes that result from human cardiac myosin mutations and a small molecule effector.***][1]

[1]: http://www.cell.com/cell-reports/abstract/S2211-1247(15)00381-2

**Tural Aksel**

**Last updated on 06/13/2015**

##Installation on Ubuntu 14.04##

- **IMPORTANT: If you previously installed fast following the old instructions, before proceeding delete the FAST codes on your system typing the lines below.** 
    
    ```
    cd /usr/bin
    sudo rm -f fast lima stack2tifs motility.py plotparams.py
    ```

- Install numpy, matplotlib, scipy, texlive, scikit-image. To install all the necessary packages, type the line below in the terminal.
   
     ```
     sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose texlive libav-tools python-skimage python-opencv
     ```
- To display fonts properly in FAST output, install MS fonts.
  
    ``` 
    sudo apt-get install ttf-mscorefonts-installer
    ```
- After installing MD fonts, remove font cache file for matplotlib in your home directory.
    
    ```
    rm -f ~/.cache/matplotlib/fontList.cache
    ```

- Move **FAST_linux** directory to the desired location. After you move the folder, on terminal change move to FAST_linux and run the installation script.     
    ```
    chmod +x install.sh
    ./install.sh
    ```
- After running the script don't move your FAST directory to some other location. If you move to another location, run the install.sh script in the new directory.

##Installation on Mac OS X Yosemite##

- **IMPORTANT: If you previously installed fast following the old instructions, before proceeding delete the FAST codes on your system typing the lines below.** 
    
    ```
    cd /usr/local/bin
    sudo rm -f fast lima stack2tifs
    cd ~/anaconda/lib/python2.7/site-packages
    rm -f motility.py plotparams.py
    ```

- MacOS X comes with python installed, however I personally prefer the [Anaconda](http://continuum.io/downloads) distribution. Download the 64 bit installer for Python2.7 and install it.

- FAST depends on some functions from [FreeImage](http://freeimage.sourceforge.net/) library. To install FreeImage, first install [macports](https://www.macports.org/). 

- Download the package for [Mac OS X 10.10 Yosemite](https://distfiles.macports.org/MacPorts/MacPorts-2.3.3-10.10-Yosemite.pkg).

- After installation of macports, install free image on terminal.

    ``` 
    sudo port install freeimage
    ```
- Install python package for openCV on terminal.

    ``` 
    conda install opencv
    ```
- Numpy package will be downgraded. Install again the updated version of numpy on terminal.

    ```
    conda install numpy
    ```
- Move **FAST_mac** directory to the desired location. After you move the folder, on terminal change move to FAST_mac and run the installation script. 
    
    ```
    chmod +x install.sh
    ./install.sh
    ```

- After running the script don't move your FAST directory to some other location. If you move to another location, run the install.sh script in the new directory.

###Installation instructions for generating tracking movies on Mac OSX###

- To create filament tracking movies as a visual report for quality of filament tracking, two extra packages need to be installed.

- Install ImageMagick by downloading and installing [ImageMagick package](http://cactuslab.com/imagemagick/assets/ImageMagick-6.9.1-0.pkg.zip).

- Install ffmpeg by downloading the binaries for your system from [ffmpeg page](http://ffmpegmac.net/).

- After downloading and opening the zip file, move the contents to ```/usr/local/bin``` on a terminal.

- First ```cd``` to the package directory and move contents to destination.
    
    ```
    sudo mv ffmpeg ffprobe ffserver /usr/local/bin
    ```

##Preparation of movie files##

- **fast** only analyzes movie tif files recorded using  [micro-manager](https://www.micro-manager.org/). For movies, recorded using other software, first save the movie as tiff stacks and convert the stacks to micro-manager output format using **stack2tiffs**.
   
     ```     stack2tiffs -d DIRECTORY -f FRAMERATE -s SIZELOWERBOUND
     ```
    
    - **DIRECTORY** is the top directory in which tiff stacks are stored.
    - **FRAMERATE** is the frame rate of the movies in frame per second **(Default: 1)**. Process movies with different frame rates separately.
    - **SIZELOWERBOUND** is the lower bound for the size (Mbytes) of the tiff stacks to be converted into individual tiffs **(Default: 6)**. Only tiffstacks bigger in size than SIZELOWERBOUND are processed.

##Analysis of movies using FAST##
- Although not necessary, it is recommended to organize the movies to be analyzed in a hierarchical order.
   
     - LEVEL1 (e.g. date)
        - LEVEL2 (e.g. slide number)
            - LEVEL3 (e.g. experimental condition)
                - LEVEL4 (e.g. replicates)   
 
- All **fast** needs is the top directory the movie folders are located at.

    ```    fast -d LEVEL1
    ```

- **FAST** first finds the lowest LEVEL directories that have movie folders under **LEVEL1**, and analyzes them in order. The lowest level movies (folders) under the same directory are treated as replicates. The results from replicates are combined to determine the average results. Therefore, it is important that the replicates have identical frame rates. Please check example movie files in the **examples/unloaded_motility** directory.

- **FAST** accepts various parameters for comprehensive analysis of filament velocities and for display of results.
    - ``` -n  WINDOWSIZE ``` : Number of consecutive frames for velocity averaging **(Default:5)**.
    - ``` -p  PATHLENGTH ``` : Minimum length for the tracked filament paths in the analysis **(Default:5)**.
    - ``` -pt TOLERANCE``` : Percent tolerance parameter to filter fluctuating velocities **(Default:None)**.
    - ```-cl COLOR```: Color of the data points in velocity scatter plot **(Default:blue)**.
    - ```-fx FUNCTION```: Function to be fitted to maximal velocity data. **exp** for single exponential decay, **uyeda** for Uyeda equation and **none** for no curve fitting **(Default:none)**.
    - ```-px PIXEL```: Pixel size in nm **(Default:80.65)**.  
    - ```-ymax YMAX```: Maximum velocity in nm/s for the scatter plot **(Default:1500)**.
    - ```-xmax XMAX```: Maximum filament length in nm for the scatter plot **(Default:10000)**.

- To estimate maximum velocities TOP5% and PLATEAU, I recommend the following parameter set.
    - ``` fast -n 5 -p 10 -pt 20 -d LEVEL1```
- For loaded motility experiments, I recommend the following parameter set.
    - ``` fast -n 5 -p 10 -d LEVEL1 ```
- Analysis results are stored in **outputs** folder in the path FAST is executed. Analysis results with different parameter sets are stored in different folders. For example, the results for **LEVEL1** analyzed using the parameters ``` -n 5 -p 10 and -pt 20``` are stored in **outputs/LEVEL1_n_5_p_10_pt_20**. Combined results from replicates at the lowest level (LEVEL4) are stored in a subfolder called **combined**.

- To analyze the movies with a new parameter set, use ```-r``` flag for speedy analysis.
    - ```fast -r -n 10 -p 10 -pt 20 -d LEVEL1```
- To force re-analyze the movies by processing through individual images, use ``` -f ``` flag.
   -  ```fast -f -n 10 -p 10 -pt 20 -d LEVEL1```
 
 - To make tracking movies, use ``` -m ``` flag.
 
     - ``` fast -m -n 10 -p 10 -pt 20 -d LEVEL1```   

- To abort execution, press ```CTRL+C``` on terminal.

- Please check the examples in **examples/unloaded_motility** to get familiar with **stack2tiffs** and **fast**.

##Result descriptions##

- **fast** plots velocities as png files and prints velocity data as text files. Complete list of unfiltered velocity points are saved with the extension ```*_full_length_velocity.txt```. Maximum path velocities, which are colored in the scatter plot, are saved with the extension ```*_max_length_velocity.txt```. The plots are saved with the extension ```*_length_velocity.png```. Combined results are saved in ```combined``` folder in ```outputs``` directory.   

- First column in ```*_length_velocity.txt``` files is the filament length in nm. Second column is the mean velocity over the ```n``` frame window (see above -n WINDOWSIZE). Third column is the standard deviation of velocities within ```n``` frame window. Fourth column is the length of the track from which the velocity is measured. 
- For description of ```*_length_velocity.png``` and the algorithms of **fast**, see [**Aksel et al. 2015**][1] 

- ```*_paths_2D.png``` shows the tracks for each filament ad the number is the average velocity for each filament track in nm/s.

- Tracking movies are saved as ```*_filament_tracks.avi``` if ```-m``` is used in fast execution. Please remember that movies will be generated, if only the packages required for movie generation are installed.

[1]: http://www.cell.com/cell-reports/abstract/S2211-1247(15)00381-2

- In addition, mean and standard error of mean (SEM) for the velocity parameters are stored in **MEAN_values.txt** and **SEM_values.txt** in **combined** folder.

##Loaded in vitro motility analysis##

- FAST is designed for high throughput analysis of loaded in vitro motility movies. For the experimental setup and the details of the loaded motility analysis please read through [our paper][1].

[1]: http://www.cell.com/cell-reports/abstract/S2211-1247(15)00381-2

- To extract the "force" parameter from a set of data collected at different utrophin (or any other actin binding protein) concentrations, I wrote a python script called **lima**. LIMA stands for Loaded In vitro Motility Analysis.

- To use lima for loaded motility analysis, user has to name the movie files in a specific format.

- **LEVEL3** (described above) should be minimally named in the following way ```PROTEINNAME_XnM_utr```. ```X``` is the utrophin concentration. For example, for a movie recorded at 0.5 nM utrophin for a myosin called **alpha**, I would name the LEVEL3 folder as ```alpha_0.5nM_utr```. For LEVEL3 and hierarchical organization of the movie folders, see above. For an example set of loaded motility data, check under ```examples/loaded_motility``` directory.

- To run **lima**, on a set of movies processed by **fast**, first go to outputs directory where the results for the complete data set are stored. For example, if user is in ```examples/loaded_motility```, enter in terminal ```cd outputs``` to change directory to outputs.

- To perform a loaded motility analysis for a **FOLDER** in ```outputs``` directory, enter in terminal,

    ```lima -d FOLDER```
- Analysis results will be stored in ```FOLDER/combined/lima```.

- For the analysis of an example data set, check ```examples/loaded_motility```.
    - First, analyze the movies:
    ```fast -r -d 032714```
    - Move to outputs folder:
    ```cd outputs```
    - Process the only directory in **outputs**: 
    ```lima -d 032714__pt_none__n_5__ymax_1500__p_5__fx_none```
    - Check the analysis results under
    ```032714__pt_none__n_5__ymax_1500__p_5__fx_none/combined/lima```. 
   
- For different analysis options, enter ```lima -h```.

        
##FAQ##
- For questions and to report bugs, please contact me by turalaksel[at]gmail.com.
