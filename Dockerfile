FROM crazyhsu/mrbigr_env:v1.1
MAINTAINER xufeng <crazyhsu9527@gmail.com>

ENV CONDA_DEFAULT_ENV=mrbigr
ENV CONDA_PREFIX=/root/conda/envs/$CONDA_DEFAULT_ENV
ENV CONDA_AUTO_UPDATE_CONDA=false

RUN echo ". /root/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate mrbigr" >> ~/.bashrc

RUN $CONDA_PREFIX/bin/pip install setuptools==67
# RUN git clone https://github.com/CrazyHsu/MRBIGR.git /root/MRBIGR && \
RUN git clone https://gitee.com/crazyhsu/MRBIGR.git /root/MRBIGR && \
    cd /root/MRBIGR && \
    $CONDA_PREFIX/bin/python setup.py build && \
    $CONDA_PREFIX/bin/python setup.py install

# RUN git clone https://github.com/liusy-jz/MODAS /root/MODAS && \
RUN git clone https://gitee.com/crazyhsu/modas.git /root/MODAS && \
    cd /root/MODAS && \
    $CONDA_PREFIX/bin/python setup.py build && \
    $CONDA_PREFIX/bin/python setup.py install

WORKDIR /root/MRBIGR

RUN chmod +x MRBIGR.py

ENV PATH=$CONDA_PREFIX/bin:/root/MRBIGR/utils:$PATH
ENV LD_LIBRARY_PATH=/root/MRBIGR/utils/libs:$LD_LIBRARY_PATH

