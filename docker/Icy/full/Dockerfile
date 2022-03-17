FROM ferreol/icy:latest

MAINTAINER Ferreol Soulez ferreol.soulez@epfl.ch

ENV  XMS=10G 

# do all in one step
RUN  	cd /usr/src/icy && \
	java -jar icy.jar -hl -x plugins.ferreol.icyhlplugininstaller.IcyHLPluginInstaller --all && \  
	cd /usr/src/icy/lib/ && \
        rm -rf  mac64 unix32 win32 win64 
