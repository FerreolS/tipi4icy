tipi4icy
========

Master: ![Travis build status](https://travis-ci.org/FerreolS/tipi4icy.svg?branch=master)  [![DOI](https://zenodo.org/badge/25529468.svg)](https://zenodo.org/badge/latestdoi/25529468)


tipi4icy is a collection of classes to link  [Icy](http://icy.bioimageanalysis.org/)  and  [TiPi](https://github.com/emmt/TiPi), a Java *Toolkit for Inverse Problems and Imaging*. Along with [microTiPi](https://github.com/FerreolS/microTiPi) is mainly used to implements the DEconvolution MICroscopy Suite.

DEconvolution MICroscopy Suite (DEMICS)
========

DEMICS is a collection of Icy plugins for deconvolution microscopy. It  provides user friendly and state of the art blind deconvolution methods for  fluorescence microscopy.

**SimpleDEMIC**

SimpleDEMIC is a simple non blind deconvolution plugin that uses total variation regularization. As it is non blind, it requires the knowledge of the PSF.

**EpiDEMIC**

EpiDEMIC stands for Epifluorescence DEconvolution MICroscopy. It is a blind deconvolution plugin for Epifluorescence / Widefield fluorescence microscopy. It requires only basic knowledge about the data (numerical aperture, wavelength, pixel size,...)

Automated batch processing
========

**Protocols**

Both plugins are compatible with Icy's protocols and can be easily added to any processing workflow.

**DEMICS on server**


Deconvolution is a memory and computation intensive task and it may not fit in end-users laptop. To process batches of large datasets, both plugins can be used without the Icy GUI on any server with [Docker](https://www.docker.com/). The docker images of [SimpleDEMIC](https://hub.docker.com/r/ferreol/simpledemic/) and [EpiDEMIC](https://hub.docker.com/r/ferreol/epidemic/) can be pulled from docker hub and executed using (for epidemic):
```
docker run --rm -v FOLDER:/data -e NAME='DATAFILE'   -e XMS=MEMORY ferreol/epidemic'
 ```
  where:
  - `FOLDER` is the folder were the dataset is,
  - `DATAFILE` is the name of the file and 
  - `MEMORY` is the memory size of the JVM (default: 10G). 
  
 All the parameters should be stored in an XML file `datafile-param.xml` in the same folder. This parameter file can be generated easily using the plugin in the Icy GUI. The result will be stored in the file `datafile-dec.tif`.


Developers corner
=======

**Eclipse**

This project can be used in Eclipse, just copy/clone it into Eclipse workspace, then from Eclipse: 

```
File -> import -> Existing project in workspace and choose TiPi folder.
```

Then Eclipse should automatically accept it as a known project.

To add the needed libraries to TiPi Eclipse project, move all the jars into a folder. Then in eclipse:

```
Right clic on TiPi project in package explorer -> Properties -> Java Build Path -> Libraries -> Add External jar and select all Library jar
```

**Icy in Eclipse**

There is a plugin to launch quickly and easily Icy from eclipse [here](http://icy.bioimageanalysis.org/index.php?display=startDevWithIcy).


## Javadoc

The javadoc is auto generated at each push and is [HERE](http://ferreols.github.io/tipi4icy)

