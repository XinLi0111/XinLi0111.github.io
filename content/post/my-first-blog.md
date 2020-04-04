---
title: "My first blog"
date: 2020-03-26T05:27:53Z
draft: false
---
idba_ud 是一个可以针对不同测序深度的短reads的基于交互式De Bruijin作图的从头拼接软件。他从小的k-mer开始到大的的k-mer步步前进，设定阈值，短的和低深度的contigs被删掉，这样来完成低深度和高深度的拼接。
拼接之前需要对原始下机数据进行清洗质控。
拼接后需要使用诱饵序列调取最终的拼接结果，所以这里将拼接过程和调取脚本合并于同一个脚本中。


1. fastp质控

通常来说，二代测序公司返回的测序数据为clean data，即经过去除低质量测序数据后的数据。目前有公司已经不再返回clean data，而是直接返回下机数据raw data。这里将这一步也加入处理流程中。
生物信息中，一旦涉及到 fastq 数据的处理，通常都要消耗大量内存、CPU 和时间。从 fastq 文件中随机抽取数据。使用的工具是 seqtk，详细安装教程：http://www.genek.tv/article/28

使用fastp软件。

它可以仅仅扫描 FASTQ 文件一次，就完成比FASTQC + cutadapt + Trimmomatic 这三个软件加起来还多很多的功能，而且速度上比仅仅使用 Trimmomatic 一个软件还要快 3 倍左右，因为它使用 C++开发，处处使用了高效算法，而且完美支持多线程！

安装

安装 fastp 十分简单，如果你使用的是 Linux 系统，可以直接使用网站上提供的预编译好的版本，下载地址是http://opengene.org/fastp/fastp，可以使用 wget等命令进行下载，下载了之后记得使用 chmod a+x ./fastp 增加该文件的可执行权限，然后就可以使用了。

也可以从源代码进行编译，需要使用 git 工具或者直接在网页上下载 release 的 源代码，以 git 下载最新的代码为例：

git clone https://github.com/OpenGene/fastp.git
cd fastp
make
sudo make install

使用
