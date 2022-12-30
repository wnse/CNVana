# 建立conda环境

FROM centos:latest

MAINTAINER yangkai07@gmail.com

## SET WORKING DIRECTORY
WORKDIR /cnv_ana



COPY . /cnv_ana
ENV PATH /opt/conda/bin:$PATH

RUN cd /etc/yum.repos.d/ \
&& sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-* \
&& sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-* \
&& yum install wget -y \
&& yum install which -y \
&& yum install libXrender-0.9.10-7.el8.i686 -y \
&& yum install cairo-devel -y \
&& wget -O /etc/yum.repos.d/CentOS-Base.repo https://mirrors.aliyun.com/repo/Centos-vault-8.5.2111.repo \
&& yum makecache

## TIMEZONE
RUN cp /usr/share/zoneinfo/Asia/Shanghai /etc/localtime
RUN echo Asia/Shanghai > /etc/timezone

RUN cd /cnv_ana \
&& sh Miniconda3-py38_4.12.0-Linux-x86_64.sh -b -p /opt/conda \
&& rm -rf Miniconda3-py38_4.12.0-Linux-x86_64.sh \
&& ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
	&& echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
&& /opt/conda/bin/conda clean -afy \
&& source /opt/conda/bin/activate \
&& conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ \
&& conda env create -f env_config.yml


ENV PATH /opt/conda/envs/pgt/bin:$PATH
ENV CONDA_DEFAULT_ENV pgt

# SHELL ["conda", "run", "-n", "pgt", "/bin/bash", "-c"]

RUN R -e "install.packages('BiocManager',repos='http://cran.us.r-project.org')" \
&& R -e "BiocManager::install('DNAcopy')" \
&& R -e "install.packages(c('ggplot2','cowplot','hash','gridExtra'),repos='http://cran.us.r-project.org')" \
&& rm -rf /tmp/downloaded_packages/

#CMD ["/bin/bash"]
#CMD ["source /opt/conda/bin/activate"]
CMD ["/opt/conda/envs/pgt/bin/gunicorn", "-c", "./gunicorn.conf.py", "wsgi:app"]
