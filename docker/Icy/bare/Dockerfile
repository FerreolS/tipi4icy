FROM anapsix/alpine-java

MAINTAINER Ferreol Soulez ferreol.soulez@epfl.ch

ENV  XMS=10G 

# do all in one step
RUN  	apk add  --update openssl && \
   	mkdir /usr/src && \ 
	cd /usr/src && \
	wget https://github.com/FerreolS/Icy-App/archive/master.zip && \
	unzip master.zip && \
	rm master.zip && \
	mv Icy-App-master icy && \ 
	wget https://github.com/FerreolS/IcyHLPluginInstaller/archive/jar.zip &&\ 
	unzip jar.zip && \
	rm -rf /usr/src/icy/plugins/* && \
	mv IcyHLPluginInstaller-jar/IcyHLPluginInstaller.jar icy/plugins/. &&\ 
	rm -r IcyHLPluginInstaller-jar jar.zip &&\
	cd /usr/src/icy/lib/ && \
	rm -rf  mac64 unix32 unix64 win32 win64 xuggler.jar && \
	apk del openssl 
