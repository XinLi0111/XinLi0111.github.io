---
title: "Idba_ud拼接昆虫线粒体全基因组脚本"
date: 2020-04-04T18:38:14+08:00
draft: false
---

>  [idba_ud](http://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/) 是一个可以针对不同测序深度的短reads的基于交互式De Bruijin作图的从头拼接软件。他从小的k-mer开始到大的的k-mer步步前进，设定阈值，短的和低深度的contigs被删掉，这样来完成低深度和高深度的拼接。
> 拼接之前需要对原始下机数据进行清洗质控。
> 拼接后需要使用诱饵序列调取最终的拼接结果，所以这里将拼接过程和调取脚本合并于同一个脚本中。


<br />


# 1. fastp质控


> 通常来说，二代测序公司返回的测序数据为clean data，即经过去除低质量测序数据后的数据。目前有公司已经不再返回clean data，而是直接返回下机数据raw data。这里将这一步也加入处理流程中。
> 生物信息中，一旦涉及到 fastq 数据的处理，通常都要消耗大量内存、CPU 和时间。从 fastq 文件中随机抽取数据。使用的工具是 [seqtk](https://github.com/lh3/seqtk)，详细安装教程：[http://www.genek.tv/article/28](http://www.genek.tv/article/28)



- 使用fastp软件。



> 它可以仅仅扫描 FASTQ 文件一次，就完成比FASTQC + cutadapt + Trimmomatic 这三个软件加起来还多很多的功能，而且速度上比仅仅使用 Trimmomatic 一个软件还要快 3 倍左右，因为它使用 C++开发，处处使用了高效算法，而且完美支持多线程！




## 安装

<br />安装 fastp 十分简单，如果你使用的是 Linux 系统，可以直接使用网站上提供的预编译好的版本，下载地址是http://opengene.org/fastp/fastp，可以使用 wget等命令进行下载，下载了之后记得使用 chmod a+x ./fastp 增加该文件的可执行权限，然后就可以使用了。<br />
<br />也可以从源代码进行编译，需要使用 git 工具或者直接在网页上下载 release 的 源代码，以 git 下载最新的代码为例：<br />

```shell
git clone https://github.com/OpenGene/fastp.git
cd fastp
make
sudo make install
```



## 使用


```shell
fastp -i in.R1.fq -o out.R1.fq -I in.R2.fq -O out.R2.fq -q 20 -w 30
# 可以看到，-i 和-o 还是用来指定 read1 的输入了输出，而大写的-I 和-O(不是零)则是用于指定read2的输入和输出，其他都保持不变。而且fastp软件最初的研发就是为了更好地处理 PE 数据，所以对于 PE 数据开发了更多的算法，比如基于 overlap 分析进行碱基校正等功能，就是只有 PE 数据独享的。fastp 对于输入和输出都支持 gzip 压缩，使用方法也很简单，只要文件名的末尾带有.gz，就会被认为是 gzip 压缩文件，会启用 gzip 对输入输出进行压缩和解压处理，例如以上 PE 的例子，如果是压缩的，就可以是以下命令：
# fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz -q 20 -w 30
```

---


# 2. idba_ud组装


- 线粒体使用idba_ud, celera，idba，ray和_de novo_ TRANS。如果是单种建库，还可以用mitoZ，一键获取线粒体组装注释上传结果。




## 安装



##### 第一步：下载安装idba

<br />这里提供在GitHub中的下载地址：[点这里下载idba](https://github.com/loneknightpy/idba/releases)<br />
<br />Ubuntu系统中使用如下命令安装：<br />

```shell
# 进入idba安装文件所在文件夹内操作
./configure
make
    
# 如果需要清除configure文件
make clean
```

<br />   idba软件初始设计是支持kmer值为128的，也就是说目前二代测序结果读长都超过其设计阈值，因此我们需要做一些更改，以适应150和250的测序数据拼接。共需要改两处：<br />
<br />

```perl
1. idba安装文件夹下的src/sequence/short_sequence.h文件中，改为kMaxShortSequence = 512；
2. basic文件夹中的kmer.h文件中，更改kNumUnit64 = 16.
```


> 作者原文：Note that IDBA assemblers are designed for short reads (around 100bp). If you want to assemble paired-end reads with longer read length, please modify the constant **kMaxShortSequence** in src/sequence/short_sequence.h to support longer read length.



> 



## 使用
在执行完质控过程后，在数据所在文件夹中运行idba_ud-assemble-0.3.5.py<br />

> 通常需要一天左右时间。该种方法极耗计算资源，要求200g以上内存，500g以上硬盘空间。


<br />

```shell
# 如果你对上述配置和编程脚本无爱，请在安装完成blast之后，手动执行以下命令。在组装之后，使用BioEdit执行blast调取。

#首先执行质控，得到解压后的fq文件。
mv *1.fq R1.fq
mv *2.fq R2.fq

fq2fa --merge R1.fq R2.fq merged.fa --filter 

# 本参数适用于250测序，如果是150测序，将--mink 40 --maxk 140用为参数
idba_ud -r merged.fa --mink 80 --maxk 240 --step 20 -o OutFolderName --num_threads 30 --min_contig 500 --similar 0.98
#在此步之后，将组装结果和你自己的诱饵序列baits，手动进行blast（BioEdit）

#如果idbassemble.py脚本出错，seqspicker.py正常，可以正常使用本套以下命令
makeblastdb -in COI.fa -dbtype nucl -parse_seqids -out COIndex

blastn -query scaffold.fa -out seq.blast -db COIndex -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30

python seqspicker.py scaffold.fa seq.blast
```



---

# 3. SequencePicker调取

<br />这一步主要是从已经组装好的scaffold中，找到你需要的那一个。<br />
<br />一般线粒体全基因组处理方法为使用blast，使用COI, 16S之类的片段，从scaffold中匹配目标序列，选取blast结果最佳的那一个。<br />


## 安装

### 第一步：下载安装blastn+

<br />NCBI网站上下载安装blastn+：[点击这里下载blastn+]()<br />
<br />下载完成后，解压，添加到环境变量中就可以使用了。<br />


### 第二步：调取target序列脚本

<br />使用Python脚本调取target序列 [点这下载此脚本](https://github.com/XinLi0111/Python/blob/master/seqspicker.py)<br />

```python
#!/usr/bin/env python
# -*- ecoding: utf-8 -*-

import pandas as pd
import sys
    
    
'''
Author:  Eric Lee
date:    2018-04-24
fileName:seqspicker.py
version: 0.2
    
this script only used for simplly pick up sequence blast results 
   from IDBA_UD scaffold fasta file.
    
usage: python seqspicker.py scaffold.fa seqblastout
'''



def parse_blastn(blastn6):
   '''parse blastn tabular results'''
   with open(blastn6) as fi:
       blast_out_list = pd.read_csv(fi, sep='\t',header = None)
       scaffold_name = list(blast_out_list[0])
       index_name = list(blast_out_list[1])
       index_dict = dict(zip(scaffold_name,index_name))
       
   return index_dict



def read_fasta(scaffold_fa):
   '''read fasta file as dict'''
   fasta = {}
   with open(scaffold_fa) as file_one:
       for line in file_one:
           line = line.strip()
           if not line:
               continue
           if line.startswith(">"):
               active_sequence_name = line[1:]
               if active_sequence_name not in fasta:
                  fasta[active_sequence_name] = []
               continue
           sequence = line
           fasta[active_sequence_name].append(sequence)

   return fasta



def get_hits(scaffold, seq_blast):
   '''get seq targets'''
   index = parse_blastn(seq_blast)
   fasta = read_fasta(scaffold)
   hit = {}
   with open('hits.fas','w') as hits:
       for key, value in index.items():
           newkey = key+'_'+value
           hit[newkey] = ''.join(fasta[key])

       for newkey,hit_value in hit.items():
           hits.write('>{}\n{}\n'.format(newkey,hit_value))



if __name__ == '__main__':

   get_hits(sys.argv[1],sys.argv[2])
   print('\nCompleted!\n')
```


### 或者：使用Python整合拼接流程


```python
#!/usr/bin/env python
# -*-coding: utf-8-*-

'''

This script used for idba_ud NGS sequences assemble.

# Author:   Eric Lee
# Date:     2020-01-16
# Version:  0.3.5

'''
import time
import os
import base64
import getpass
import argparse
import subprocess
from glob import glob
from multiprocessing import cpu_count
try:
    import yagmail
except ImportError:
    subprocess.call('pip install yagmail', shell=True)
    import yagmail

try:
    from tqdm import tqdm
except ImportError:
    subprocess.call('pip install tqdm', shell=True)
    from tqdm import tqdm

try:
    from colorama import init, Fore, Back
except ImportError:
    subprocess.call('pip install colorama', shell=True)
    from colorama import init, Fore, Back

try:
    import pandas as pd
except ImportError:
    subprocess.call('pip install pandas', shell=True)
    import pandas as pd


'''
History:
0.1   create
0.1.1 add unzip function
      add choice of read length (150/250)
0.2   rewrite with def funtion, print color symbol notice
0.2.1 add prcessing state
0.2.2 add judge function to process the baits input filenames
0.2.3 fix the fq2fa parameters mistakes R1 R2 is not global variations
0.2.4 support extract more type compressed input_files
0.2.5 integrated seqspicker function in one script
0.2.6 add colorama to show colors in windows console
0.2.7 add end color
0.2.8 make import more robust
0.3   TODO：用click获取参数，tqdm创建进度条，最后建立一个各步骤名称的元组set，循环，
      依据参数判定执行哪一步
0.3.1 去除click，逻辑太难弄，用argparse重写
0.3.2 增加任务完成后可选邮件通知
0.3.3 排除了tqdm进度条末尾残留字符的问题。中文字符导致，英文无问题
0.3.4 调整tqdm显示设置
0.3.5 增加文件检测，避免运行到中间因缺失文件而中断
'''

parser = argparse.ArgumentParser(prog='idba-assemble.py', usage='python %(prog)s ' +
                                 'fastq -i R1.fq R2.fq -r 150 -b COI.fa',
                                 description='原始数据需为clean data')
parser.add_argument('-v', '--version', action='version', version='%(prog)s v 0.3.5')
parser.add_argument('filetype', default='fastq', metavar='bz2|targz|gz|fastq',
                    choices=['bz2', 'targz', 'gz', 'fastq'],
                    help='输入文件格式，如为压缩文件将自动解压')
parser.add_argument('-i', '--inputs', nargs=2, required=True,
                    metavar='*-clean.fq', help='输入文件名称[双端测序结果文件]')
parser.add_argument('-r', '--read', choices=[150, 250], type=int,
                    metavar='150|250', required=True,
                    help='输入双端测序读长[150/250]')
parser.add_argument('-b', '--baits', nargs='+', metavar='COI.fa',
                    required=True, help='输入诱饵序列文件名称，至少一个，可多个')
parser.add_argument('-m', '--mail', nargs='?', default=False,
                    const='xxxx@126.com', help='接收任务完成的邮箱地址') 

args = parser.parse_args()

# init(autoreset=True)  # colorama set


def green(log):
    'print logs with green color'
    return Fore.GREEN + log + Fore.RESET


def red(error):
    'print error with red color'
    return Fore.RED + error + Fore.RESET
    

def yellow(notice):
    'print notice with yellow color'
    return Fore.YELLOW + notice + Fore.RESET


def strong_notice(notice=['strong']):
    'notice in center and red background'
    def cprint(repeat, char='#'):
        print(Back.GREEN + char * repeat + Back.RESET)

    width = os.get_terminal_size()[0]
    cprint(width)
    for line in notice: 
        print(Back.GREEN + line.center(width) + Back.RESET)
    cprint(width)


cpu_number = cpu_count()  # get & check cpu numbers
def check_cpu():
    if cpu_number < 8:
        print(red(f'NOTICE: The CPU number is {cpu_number}, it is not ' +
                  'recommanded to run following steps!'))
        os._exit(0)


def check_file_existence():
    'check whether the files are existed in current work folder'
    file_names = args.baits + args.inputs
    exit = False
    for file_name in file_names:
        if not os.path.exists(file_name):
            print(red(f'{file_name} can not be found!'))
            exit = True
    if exit:
        os._exit(0)


# mail sender
def triple_try_mail(passwd_salt):
    '''
    try send mail three times
    :passwd_salt: ''
    '''
    
    passwd_encrypt = 'xxxxxxxxxxxxxx'
    for i in passwd_salt:
        passwd_encrypt = passwd_encrypt.replace(i, '', 1)
    passwd_decrypt = base64.decodebytes(passwd_encrypt.encode()).decode()
    yag = yagmail.SMTP(user='xxxxx@qq.com',
                       password=passwd_decrypt,
                       host='smtp.qq.com', port='465')  # 验证邮件服务器
    attempts = 0
    success = False
    while attempts < 3 and not success:
        try:
            if os.path.exists('hits.fas'):
                yag.send(args.mail, subject = "idba_ud task completed",
                         contents = "<h2>idba_ud assemble mission completed!</h2>",
                         attachments = 'hits.fas')
            else:
                yag.send(args.mail, subject = "idba_ud task completed",
                         contents = "<h2>idba_ud assemble mission completed!</h2>")
            success = True
            print(green('邮件发送成功'))
        except:
            attempts += 1
            if attempts == 3:
                print(red('邮件发送失败'))
                os._exit(0)


# seqspicker
class SequencePicker(object):
    '''
    extract target sequences from the completed scaffold assembled by idba_ud
    Usage:
    SequencePicker.get_hits('./OutFolder/scaffold.fa', blast_out_filename)
    '''
    def __init__(self):
        pass

    def _parse_blastn(self, blastn6):
        '''parse blastn tabular results'''
        with open(blastn6) as fi:
            blast_out_list = pd.read_csv(fi, sep='\t', header=None)
            scaffold_name = list(blast_out_list[0])
            index_name = list(blast_out_list[1])
            index_dict = dict(zip(scaffold_name, index_name))

        return index_dict

    def _read_fasta(self, scaffold_fa):
        '''read fasta file as dict'''
        fasta = {}
        with open(scaffold_fa) as file_one:
            for line in file_one:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    active_sequence_name = line[1:]
                    if active_sequence_name not in fasta:
                        fasta[active_sequence_name] = []
                    continue
                sequence = line
                fasta[active_sequence_name].append(sequence)

        return fasta

    @classmethod
    def get_hits(cls, scaffold, seq_blast):
        '''get seq targets'''
        index = cls()._parse_blastn(seq_blast)
        fasta = cls()._read_fasta(scaffold)
        hit = {}
        with open('hits.fas', 'w') as hits:
            for key, value in index.items():
                newkey = key + '_' + value
                hit[newkey] = ''.join(fasta[key])

            for newkey, hit_value in hit.items():
                hits.write(f'>{newkey}\n{hit_value}\n')


# 第一步 输入
class InFiles(object):
    '''
    获取fastq文件
    obj = InFiles.getfastq(args.filetype, *args.inputs)
    '''
    #def __init__(self, filetype, *input_files):
    #    self.filetype = filetype
    #    self.input_files = input_files

    def _bz2(self, *input_files):
        'two tar.bz2 input_files name <*.tar.bz2>'
        for bz2_file in input_files:
            subprocess.call(f'tar -jxvf {bz2_file}', shell=True)
        fq_files = glob(r'*.fq')
        fq_files.extend(glob(r'*.fastq'))

        return fq_files

    def _targz(self, *input_files):
        'two tar.gz input_files name <*.tar.gz>'
        for tar_file in input_files:
            subprocess.call(f'tar -zxvf {tar_file}', shell=True)
        fq_files = glob(r'*.fq')
        fq_files.extend(glob(r'*.fastq'))

        return fq_files

    def _gz(self, *input_files):
        'two gz input_files name <*.gz>'
        for gz_file in input_files:
            subprocess.call(f'gzip -d {gz_file}', shell=True)
        fq_files = glob(r'*.fq')
        fq_files.extend(glob(r'*.fastq'))

        return fq_files

    def _fastq(self, *fq_files):
        'two fastq input_files name <*.fq>/<*.fastq>'

        return fq_files

    @classmethod
    def get_fastqs(cls, filetype, *input_files):
        'input func'
        param_dict = {'bz2': cls()._bz2,
                      'targz': cls()._targz,
                      'gz': cls()._gz,
                      'fastq': cls()._fastq}
        tqdm.write(f'Deal with {filetype} files: {yellow(input_files[0])} {yellow(input_files[1])}')
        fastq_files = param_dict.get(filetype, -1)(*input_files)

        return fastq_files


# 第二步 设定参数
def set_kmers():
    'get and check read length parameters'
    kmer_parameters = ''  # set kmers
    if args.read > 151:
        kmer_parameters = '--mink 80 --maxk 240'
        tqdm.write(f'Read length is {args.read}, and kmers setted as "{yellow(kmer_parameters)}".')
    else:
        kmer_parameters = '--mink 40 --maxk 140'
        tqdm.write(f'Read length is {args.read}, and kmers setted as "{yellow(kmer_parameters)}".')

    return kmer_parameters


# 第三步 组装
def assemble(kmer_parameters, *fq_files):
    'idba_ud组装过程'
   
    merge = f'fq2fa --merge {fq_files[0]} {fq_files[1]} merged.fsa --filter'
    assembly = f'idba_ud -r merged.fsa {kmer_parameters} --step 20 -o OutFolder ' + \
               f'--num_threads {cpu_number - 1} --min_contig 500 --similar 0.98'
    assemble_parameters = [merge, assembly]
    assemble_tqdm = tqdm(assemble_parameters, desc='  Assembling', position=0,
                         bar_format='{desc}: {n_fmt}/{total_fmt}|{bar}|')
    for step in assemble_tqdm:
        tqdm.write(step)  # test
        subprocess.call(step)


# 第四步 调取
def blast_hits(baits):
    baits_tqdm = tqdm(baits, desc='  Pickup', position=1,
                      bar_format='{desc}: {n_fmt}/{total_fmt}|{bar}|')    
    for bait in baits_tqdm:
        blast_out_filename = f'{bait}.blastout'
        makeblastdb = f'makeblastdb -in {bait} -dbtype nucl -parse_seqids -out {bait}'
        blastn = f'blastn -query ./OutFolder/scaffold.fa -out {blast_out_filename} -db {bait} -outfmt 6 \
-evalue 1e-5 -max_target_seqs 1 -num_threads {cpu_number - 1}'

        blast_parameters = [makeblastdb, blastn, 'get_hits']
        for step in tqdm(blast_parameters, desc=f'    {bait} selecting', position=0,
                         bar_format='{desc}: {n_fmt}/{total_fmt}|{bar}|'):
            if step != 'get_hits':
                tqdm.write(step)   
                subprocess.call(step)
            else:
                tqdm.write(f'SequencePicker.get_hits(./OutFolder/scaffold.fa, {blast_out_filename})')  # test
                SequencePicker.get_hits('./OutFolder/scaffold.fa', blast_out_filename)



if __name__ == '__main__':

    check_cpu()  # check whether CPU is enough for analysis
    check_file_existence()
    
    if args.mail:
        passwd_salt=getpass.getpass()

    with tqdm(total=4, desc='Main progress', position=2, bar_format='{desc}: {n_fmt}/{total_fmt}|{bar}|') as pbar:
        tqdm.write(green('Getting fastqs:')) # desc参数为中文会导致显示不完美
        fq_file_names = InFiles.get_fastqs(args.filetype, *args.inputs)
        pbar.update(1)
        
        tqdm.write(green('Set kmers:'))
        setted_kmers = set_kmers()  # set kmer start, end and step
        pbar.update(1)
        
        tqdm.write(green('Start assembly:'))
        assemble(setted_kmers, *fq_file_names)
        pbar.update(1)
        
        tqdm.write(green('Start pickup:'))
        blast_hits(args.baits)
        pbar.update(1)


    if args.mail:
        print(green(f'Send mail to {yellow(args.mail)}'))
        triple_try_mail(passwd_salt)  # send mail

    # 结束
    notice_list = ['ASSEMBLE COMPLETED!',
                   'The result file is stored in `./hits.fas`']
    strong_notice(notice_list)

```

将上述软件安装成功后并配置好环境变量即可使用，若出现python报错没有依赖包，将脚本中import后的包都安装下。python版本要求大于3.


> linux环境变量一般在~/.bashrc中，使用文本编辑器（如vim）更改相关配置就好：

```shell
# add idbassemble to path
export PATH="$PATH:/path/to/idbassemble.pyfolder/"
```
> 如果出现脚本不能运行，请更改脚本文件的可执行权限``chmod +x idba-ud-assemble-0.3.5.py`
