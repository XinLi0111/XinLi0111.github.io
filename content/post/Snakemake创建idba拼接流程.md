---
title: "Snakemake创建idba拼接流程"
date: 2020-04-04T18:14:11+08:00
draft: false
---


> 上一篇文章中说到了自己写脚本创建的idba流程，但是在工作流程方面是有一些专门的软件存在。比如Next flow以及snakemake等，因为snakemake是以Python语言为主要书写语言，因此学习了snakemake。

# snakemake简介
snakemake是一个创建流程的工具，其最主要的部分就是一个个**"rule"**，一个rule就是一系列流程中的最基础的一步。每个rule简单来说就相当于一个shell脚本，主要由 `input` 、`output` 、`shell` 3个部分构成。多个rule之间由它们的 `input` 、 `output` 相互关系决定执行顺序。上一步的结果文件是下一步运行的条件。

```python
#这是一个简单的snakemake rule示例，将两个文件合并为一个文件
rule concat:
    input:
        expand("{file}.txt", file=["hello", "world"])
    output:
        "merge.txt"
    shell:
        "cat {input} > {output}"
```

# snakemake安装
之前尝试过自己安装snakemake，但是遇到了不知名错误，在手机上使用termux安装并无异常。查阅报错信息和官方指导后，建议使用[官方提供的方法](https://snakemake.readthedocs.io/en/stable/index.html)：使用conda直接创建新环境安装 `conda create -c bioconda -c conda-forge -n snakemake snakemake` ，自己新建环境再安装snakemake遭遇了依赖无法解决的情况，不知道是最近版本问题还是一直都有问题。

# 激活snakemake环境
上一步已经安装好了新的conda环境——snakemake，在终端中使用 `conda activate snakemake` 激活snakemake环境。
# snakemake参数简介
snakemake的参数还挺多的，具体的请参阅官方文档，这里仅作简介。

```python
--snakefile, -s 指定Snakefile，否则是当前目录下的Snakefile
--dryrun, -n  不真正执行，一般用来查看Snakefile是否有错
--printshellcmds, -p   输出要执行的shell命令
--reason, -r  输出每条rule执行的原因,默认FALSE
--cores, --jobs, -j  指定运行的核数，若不指定，则使用最大的核数
--force, -f 重新运行第一条rule或指定的rule
--forceall, -F 重新运行所有的rule，不管是否已经有输出结果
--forcerun, -R 重新执行Snakefile，当更新了rule时候使用此命令

#一些可视化命令
$ snakemake --dag | dot -Tsvg > dag.svg
```

# snakemake创建的idba拼接流程

```python
# filename: Snakemake.py
BAIT=['coi', 'cob']
rule all:
    input: 
        expand('bait/{bait}.hits.fas', bait=BAIT),  # 注意多个文件之间有逗号
        'idba_assembled/scaffold.fa'           # 需要最终保存的文件可以写在这里，或者使用protect保护起来


rule fq2fa:
    input:
        expand('r{index}.fq', index=['1','2'])
    output:
        'merged/merged.fsa'
    shell:
        # 'touch {output}'
        'fq2fa --merge {input} {output} --filter'

rule idba_ud:
    input:
        rules.fq2fa.output
    output:
        'idba_assembled/scaffold.fa'
    params: '--mink 40 --maxk 140'
    threads: 7          # 计算机线程数
    shell:
        # 'touch {output} coi.fa cob.fa'
        'idba_ud -r {input} {params} --step 20 -o idba_assembled \
        --num_threads {threads} --min_contig 500 --similar 0.98'


rule makeblastdb:
    input: 
        'bait/{bait}.fa'  # {bait}变量需要在涉及rule的input和output一直传承，否则在解决文件依赖时可能会出错
    output:
        'blast/{bait}.blastdb'
    shell:
        # 'touch {output}'
        'makeblastdb -in {input} -dbtype nucl -parse_seqids -out {output}'

rule blastn:
    input:    # 这里是涉及bait的rule和源文件合并rule的交点
        scaffold='idba_assembled/scaffold.fa',
        blastdb='blast/{bait}.blastdb'
    output:
        'blast/{bait}.blastout'
    params: 7
    shell:
        # 'touch {output}'
        'blastn -query {input.scaffold} -out {output} -db {input.blastdb} -outfmt 6 \
        -evalue 1e-5 -max_target_seqs 1 -num_threads {params}'

rule pick_up:
    input:
        scaffold = 'idba_assembled/scaffold.fa',
        blastout = 'blast/{bait}.blastout'
    output: 'bait/{bait}.hits.fas'
    run:         # `shell`执行shell命令，`run`直接执行Python语言，`script`执行Python脚本文件【注意传参格式】
        # 'snakemake_idba_picker.py'
        # 'python seqspicker({input.scaffold} {input.blastout})'
        
        def parse_blastn(blastn6):
            '''parse blastn tabular results'''
            with open(blastn6) as fi:
                blast_out_list = pd.read_csv(fi, sep='\t', header=None)
                scaffold_name = list(blast_out_list[0])
                index_name = list(blast_out_list[1])
                index_dict = dict(zip(scaffold_name, index_name))
                
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


        def get_hits(scaffold, seq_blast, resultsfilename):
            '''get seq targets'''
            index = parse_blastn(seq_blast)
            fasta = read_fasta(scaffold)
            hit = {}
            with open(resultsfilename,'w') as hits:
                for key, value in index.items():
                    newkey = key+'_'+value
                    hit[newkey] = ''.join(fasta[key])

                for newkey, hit_value in hit.items():
                    hits.write('>{}\n{}\n'.format(newkey, hit_value))


        get_hits(input.scaffold, input.blastout, output[0])

# python script should write like this example
# def do_something(data_path, out_path, threads, myparam):
#     # python code

# do_something(snakemake.input[0], snakemake.output[0], snakemake.config["myparam"])

```
# 运行snakemake

```shell
#参数含义见上文
#仅查看是否可以运行
snakemake -s Snakefile.py -np
#或者
snakemake -s Snakemake.py -nr
#创建流程图普
snakemake -s Snakemake.py --dag | dot -Tsvg > dag.svg
#正式运行
snakemake -s Snakemake.py -p
```

