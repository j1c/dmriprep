# Use Ubuntu 16.04 LTS
FROM ubuntu:xenial-20161213

# Pre-cache neurodebian key
COPY docker/files/neurodebian.gpg /usr/local/etc/neurodebian.gpg

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    curl \
                    bzip2 \
                    ca-certificates \
                    xvfb \
                    cython3 \
                    build-essential \
                    autoconf \
                    libtool \
                    pkg-config \
		    vim \
		    zip \
		    unzip \
		    wget \
                    git && \
    curl -sL https://deb.nodesource.com/setup_10.x | bash - && \
    apt-get install -y --no-install-recommends \
                    nodejs && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    curl -o pandoc-2.2.2.1-1-amd64.deb -sSL "https://github.com/jgm/pandoc/releases/download/2.2.2.1/pandoc-2.2.2.1-1-amd64.deb" && \
    dpkg -i pandoc-2.2.2.1-1-amd64.deb && \
    rm pandoc-2.2.2.1-1-amd64.deb

# Installing Neurodebian packages (FSL, AFNI, git)
RUN curl -sSL "http://neuro.debian.net/lists/$( lsb_release -c | cut -f2 ).us-ca.full" >> /etc/apt/sources.list.d/neurodebian.sources.list && \
    apt-key add /usr/local/etc/neurodebian.gpg && \
    (apt-key adv --refresh-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9 || true) && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
                    afni=16.2.07~dfsg.1-5~nd16.04+1 \
                    git-annex-standalone && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV FSLDIR="/opt/fsl-6.0.1" \
    PATH="/opt/fsl-6.0.1/bin:$PATH" \
    FSLOUTPUTTYPE="NIFTI_GZ"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           bc \
           dc \
           file \
           libfontconfig1 \
           libfreetype6 \
           libgl1-mesa-dev \
           libglu1-mesa-dev \
           libgomp1 \
           libice6 \
           libxcursor1 \
           libxft2 \
           libxinerama1 \
           libxrandr2 \
           libxrender1 \
           libxt6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Downloading FSL ..." \
    && wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py \
    && chmod 775 fslinstaller.py \
    && /fslinstaller.py -d /opt/fsl-6.0.1 -V 6.0.1 -q \
    && rm -rf /opt/fsl-6.0.1/data \
    && rm -rf /opt/fsl-6.0.1/bin/FSLeyes* \
    && rm -rf /opt/fsl-6.0.1/src \
    && rm -rf /opt/fsl-6.0.1/extras/src \
    && rm -rf /opt/fsl-6.0.1/doc \
    && rm -rf /opt/fsl-6.0.1/bin/fslview.app \
    && rm -rf /opt/fsl-6.0.1/data/atlases \
    && rm -rf /opt/fsl-6.0.1/data/first \
    && rm -rf /opt/fsl-6.0.1/data/mist \
    && rm -rf /opt/fsl-6.0.1/data/possum

# Installing ANTs 2.2.0 (NeuroDocker build)
ENV ANTSPATH=/usr/lib/ants
RUN mkdir -p $ANTSPATH && \
    curl -sSL "https://dl.dropbox.com/s/2f4sui1z6lcgyek/ANTs-Linux-centos5_x86_64-v2.2.0-0740f91.tar.gz" \
    | tar -xzC $ANTSPATH --strip-components 1

ENV AFNI_INSTALLDIR=/usr/lib/afni \
    PATH=${PATH}:/usr/lib/afni/bin \
    AFNI_PLUGINPATH=/usr/lib/afni/plugins \
    AFNI_MODELPATH=/usr/lib/afni/models \
    AFNI_TTATLAS_DATASET=/usr/share/afni/atlases \
    AFNI_IMSAVE_WARNINGS=NO \
    FSLOUTPUTTYPE=NIFTI_GZ \
    PATH=$ANTSPATH:$PATH \
    ANTS_VERSION=2.2.0

# Installing SVGO
RUN curl -sL https://deb.nodesource.com/setup_10.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g svgo \
    && npm install -g bids-validator@1.2.3 \
    && curl -sSLO https://repo.continuum.io/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh && \
    bash Miniconda3-4.5.11-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-4.5.11-Linux-x86_64.sh

# Set CPATH for packages relying on compiled libs (e.g. indexed_gzip)
ENV PATH="/usr/local/miniconda/bin:$PATH" \
    CPATH="/usr/local/miniconda/include/:$CPATH" \
    LANG="C.UTF-8" \
    LC_ALL="C.UTF-8" \
    PYTHONNOUSERSITE=1

# Installing precomputed python packages
RUN conda install -y python=3.7.1 \
                     pip=19.1 \
                     mkl=2018.0.3 \
                     mkl-service \
                     numpy=1.15.4 \
                     scipy=1.1.0 \
                     scikit-learn=0.19.1 \
                     matplotlib=2.2.2 \
                     pandas=0.23.4 \
                     libxml2=2.9.8 \
                     libxslt=1.1.32 \
                     graphviz=2.40.1 \
                     traits=4.6.0 \
                     zlib; sync && \
		     cython && \
    chmod -R a+rX /usr/local/miniconda; sync && \
    chmod 777 /usr/local/miniconda/bin/*; sync && \
    conda build purge-all; sync && \
    conda clean -tipsy && sync

# Unless otherwise specified each process should only use one thread - nipype
# will handle parallelization
ENV MKL_NUM_THREADS=4 \
    OMP_NUM_THREADS=4

# Create a shared $HOME directory
RUN useradd --no-user-group --create-home --shell /bin/bash dmriprep
WORKDIR /home/dmriprep
ENV HOME="/home/dmriprep"

# Precaching fonts, set 'Agg' as default backend for matplotlib
RUN python -c "from matplotlib import font_manager" \
    && sed -i 's/\(backend *: \).*$/\1Agg/g' $( python -c "import matplotlib; print(matplotlib.matplotlib_fname())" ) \
    && pip install --upgrade pip \
    && apt-get update && apt-get install -y sudo gfortran liblapack-dev libopenblas-dev \
    && pip install ipython cython parse \
    && git clone -b homecooked https://github.com/dPys/dmriprep.git dmriprep \
    && cd dmriprep \
    && python setup.py install \
    && pip install ipython cython parse \
    && pip install --no-cache-dir https://github.com/samuelstjean/nlsam/archive/master.zip

RUN chown -R dmriprep:dmriprep /usr/local/miniconda/lib/python3.7 \
    && chown -R dmriprep:dmriprep /home/dmriprep \
    && chmod -R 777 /home/dmriprep \
    && find $HOME -type d -exec chmod go=u {} + \
    && find $HOME -type f -exec chmod go=u {} + \
    && mkdir /inputs \
    && chmod -R 777 /inputs \
    && mkdir /outputs \
    && chmod -R 777 /outputs \
    && echo "dmriprep ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/user \
    && chmod a+s -R /usr/local \
    && chmod -R 775 /usr/local/miniconda/lib/python3.7/site-packages \
    && chmod -R 777 /usr/local/miniconda/lib/python3.7/site-packages/dmriprep*; sync \
    && apt-get purge -y --auto-remove \
       git \
       wget \
       curl \
       build-essential \
       ca-certificates \
       cython3

ENV IS_DOCKER_8395080871=1

RUN ldconfig
WORKDIR /tmp/
USER dmriprep
ENTRYPOINT ["/usr/local/miniconda/bin/dmriprep"]

ARG BUILD_DATE
ARG VCS_REF
ARG VERSION
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="dMRIPrep" \
      org.label-schema.description="dMRIPrep - robust dMRI preprocessing tool" \
      org.label-schema.url="http://dmriprep.org" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/nipreps/dmriprep" \
      org.label-schema.version=$VERSION \
      org.label-schema.schema-version="1.0"
