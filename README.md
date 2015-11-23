tipi4icy
========

![Travis build status](https://travis-ci.org/emmt/tipi4icy.svg?branch=master)

tipi4icy is a collection of [Icy](http://icy.bioimageanalysis.org/) plugins based on
[TiPi](https://github.com/emmt/TiPi), a Java *Toolkit for Inverse Problems and Imaging*.

Eclipse
=======

This project can be used in Eclipse, just copy/clone it into Eclipse workspace, then from Eclipse: 

```
File -> import -> Existing project in workspace and choose TiPi folder.
```

Then Eclipse should automatically accept it as a known project.

To add the needed libraries to TiPi Eclipse project, move all the jars into a folder. Then in eclipse:

```
Right clic on TiPi project in package explorer -> Properties -> Java Build Path -> Libraries -> Add External jar and select all Library jar
```

Icy in Eclipse
==============

There is a plugin to launch quickly and easily Icy from eclipse [here](http://icy.bioimageanalysis.org/index.php?display=startDevWithIcy).

Plugins
========

**MitivDeconvolution**

MitivDeconvolution is the plugin that use linear algorithms to deconvolve 2D or 3D images.

The algorithms are wiener and linear conjugate gradients.

In ICY this plugin is available as a protocol.

**MitivTotalVariation**

MitivTotalVartion use the Total variation (non linear) to deconvolve 2D or 3D images.

In ICY this plugin is available as a protocol.

**MitivGlobalDeconv**

MitivGlobalDeconv (temporary name) is the plugin that is able to do blind deconvolution for microscopy images.

In progress...

## Javadoc

The javadoc is auto generated at each push and is [HERE](http://emmt.github.io/tipi4icy)

## Credits

[TiPi](https://github.com/emmt/TiPi) and [tipi4icy](https://github.com/emmt/tipi4icy)
have been developed as part of the [MiTiV](http://mitiv.univ-lyon1.fr) project.

## Links

[Main site web](http://mitiv.univ-lyon1.fr/)

[Plugins presentation](http://mitiv.univ-lyon1.fr/software/available-plugins.html)

[Tutorial](http://mitiv.univ-lyon1.fr/software/tutorial.html)
