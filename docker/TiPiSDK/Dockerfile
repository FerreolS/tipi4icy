FROM ferreol/icy:bare

MAINTAINER Ferreol Soulez ferreol.soulez@univ-lyon1.fr

# do all in one step
RUN     apk add  --update openssl && \
        mkdir -p /usr/src && \
        cd /usr/src && \
        wget https://github.com/FerreolS/tipi4icy/archive/jars.zip && \
        unzip jars.zip  && \
        mv tipi4icy-jars/TiPiSDK.jar icy/plugins/. &&\ 
        rm -r tipi4icy-jars jars.zip &&\
	cd icy && \ 
	java -jar icy.jar -hl  -x plugins.ferreol.icyhlplugininstaller.IcyHLPluginInstaller \
	     plugins.adufour.ezplug.EzPlug \
	     plugins.adufour.blocks.Blocks \
	     plugins.mitiv.JTransforms3 \ 
	     plugins.mitiv.JLargeArrays &&\
        cd .. && \
	apk del openssl 
