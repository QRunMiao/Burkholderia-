```shell

cd /mnt/c/shengxin/data

cat Amycolatopsis/assembly/Amycolatopsis.assembly.tsv Streptomyces/assembly/Streptomyces.assembly.tsv Burkholderia/assembly/Burkholderia.assembly.tsv > all/all.assembly.tsv

nwr template /mnt/c/shengxin/data/all/all.assembly.tsv \
    --pro  #生成Protein/species.tsv ，加上Ros_intestinalis_BSD2780061689_150309_G12_GCF_015554785_1和Si_sun_NSJ_8_GCF_014337175_1	Simiaoa_sunii，共11085行


rsync -avP \
/mnt/c/shengxin/data/Protein/species.tsv \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all

#ref 的
rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/ASSEMBLY/Simiaoa_sunii/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/ASSEMBLY

rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/ASSEMBLY/Roseburia_intestinalis/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/ASSEMBLY

rsync -avP \
/mnt/c/shengxin/data/Protein/collect.sh \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/

bsub -n 24 -q mpi -J "DddA-like" "bash ~/qyl/data/all/collect.sh"
bsub -n 24 -q mpi -J "DddA-like" "gzip -d all.replace.fa.gz"

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/2replace.tsv \
    /mnt/c/shengxin/data/all

bsub -n 24 -q mpi -J "DddA-like1" "E_VALUE=1e-10 ; cat ~/data/Bacteria/Protein/species.tsv |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 20 \"
            if [[ ! -d ~/qyl/data/all/ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ~/qyl/data/all/ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E \${E_VALUE} --domE \${E_VALUE} --noali --notextw ~/qyl/data/Burkholderia/HMM/DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \\\$1, {1}; '
        \" \
        > ~/qyl/data/Burkholderia/2replace.tsv" #394

bsub -n 24 -q mpi -J "DddA-like1" "E_VALUE=1e-10 ; cat ~/qyl/data/all/species0.tsv |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 20 \"
            if [[ ! -d ~/qyl/data/all/ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ~/qyl/data/all/ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E \${E_VALUE} --domE \${E_VALUE} --noali --notextw ~/qyl/data/Burkholderia/HMM/DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \\\$1, {1}; '
        \" \
        >> ~/qyl/data/Burkholderia/2replace.tsv" #加2个

cat ~/qyl/data/Burkholderia/2replace.tsv | tsv-select -f 2,1 | tr "\t" "_"> ~/qyl/data/Burkholderia/3.tsv
rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/3.tsv \
    /mnt/c/shengxin/data/Burkholderia #在vscode中将3.tsv文件中的后缀.1去掉

rsync -avP \
/mnt/c/shengxin/data/Burkholderia/3.tsv \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/3.tsv

bsub -n 24 -q mpi -J "DddA-like1" "faops some ~/qyl/data/all/all.replace.fa ~/qyl/data/Burkholderia/3.tsv ~/qyl/data/all/Protein/all-DddA-like.fa" #396个序列

rsync -avP \
/mnt/c/shengxin/data/all/Protein/all-DddA-like.fa \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein/
    
muscle -in ~/qyl/data/all/Protein/all-DddA-like.fa -out ~/qyl/data/all/Protein/all-DddA-aln.fa

FastTree ~/qyl/data/all/Protein/all-DddA-aln.fa > ~/qyl/data/all/Protein/all-DddA-aln.newick
rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein/all-DddA-aln.newick \
    /mnt/c/shengxin/data/Burkholderia 

#系统发育树
cat ~/qyl/data/Burkholderia/2replace.tsv |cut -f 2 > ~/qyl/data/Burkholderia/1.tsv
cat ~/data/Bacteria/Protein/species.tsv |tsv-join -f ~/qyl/data/Burkholderia/1.tsv -k 1 -d 1 > ~/qyl/data/Burkholderia/2.tsv
#376个物种

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/2.tsv \
    /mnt/c/shengxin/data/Burkholderia #376物种
rsync -avP \
/mnt/c/shengxin/data/Burkholderia/2.tsv \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/2.tsv

cat ~/qyl/data/Burkholderia/2.tsv | tsv-select -f 2,1 > ~/qyl/data/Burkholderia/5.tsv #378物种

rsync -avP \
    /mnt/c/shengxin/data/bac120.sh  \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia

#bac120.sh内容: 
E_VALUE=1e-10 ; for marker in PF02576.12 PF01025.14 PF03726.9 PF00466.15 PF00410.14 PF00380.14 TIGR00006 TIGR00019 TIGR00020 TIGR00029 TIGR00043 TIGR00054 TIGR00059 TIGR00061 TIGR00064 TIGR00065 TIGR00082 TIGR00083 TIGR00084 TIGR00086 TIGR00088 TIGR00090 TIGR00092 TIGR00095 TIGR00115 TIGR00116 TIGR00138 TIGR00158 TIGR00166 TIGR00168 TIGR00186 TIGR00194 TIGR00250 TIGR00337 TIGR00344 TIGR00362 TIGR00382 TIGR00392 TIGR00396 TIGR00398 TIGR00414 TIGR00416 TIGR00420 TIGR00431 TIGR00435 TIGR00436 TIGR00442 TIGR00445 TIGR00456 TIGR00459 TIGR00460 TIGR00468 TIGR00472 TIGR00487 TIGR00496 TIGR00539 TIGR00580 TIGR00593 TIGR00615 TIGR00631 TIGR00634 TIGR00635 TIGR00643 TIGR00663 TIGR00717 TIGR00755 TIGR00810 TIGR00922 TIGR00928 TIGR00959 TIGR00963 TIGR00964 TIGR00967 TIGR01009 TIGR01011 TIGR01017 TIGR01021 TIGR01029 TIGR01032 TIGR01039 TIGR01044 TIGR01059 TIGR01063 TIGR01066 TIGR01071 TIGR01079 TIGR01082 TIGR01087 TIGR01128 TIGR01146 TIGR01164 TIGR01169 TIGR01171 TIGR01302 TIGR01391 TIGR01393 TIGR01394 TIGR01510 TIGR01632 TIGR01951 TIGR01953 TIGR02012 TIGR02013 TIGR02027 TIGR02075 TIGR02191 TIGR02273 TIGR02350 TIGR02386 TIGR02397 TIGR02432 TIGR02729 TIGR03263 TIGR03594 TIGR03625 TIGR03632 TIGR03654 TIGR03723 TIGR03725 TIGR03953; do
    >&2 echo "==> marker [${marker}]"

mkdir -p ~/qyl/data/Burkholderia/${marker}

cat ~/qyl/data/Burkholderia/5.tsv |
        parallel --colsep '\t' --no-run-if-empty  --linebuffer -k -j 8 "
            gzip -dcf ~/data/Bacteria/ASSEMBLY/{1}/{2}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/Bacteria/HMM/hmm/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {2}; '
        "         > ~/qyl/data/Burkholderia/${marker}/replace.tsv

    >&2 echo

done

rsync -avP \
    /mnt/c/shengxin/data/bac120.sh \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/bac120.sh 

bsub -n 24 -q mpi -J "DddA-like" "bash ~/qyl/data/Burkholderia/bac120.sh"

du -sh * 
#du 是表示"disk usage"的缩写，用于显示文件或目录的磁盘使用情况。
#-s 选项用于仅显示总计信息而不显示每个子目录的详细信息，这样可以更清晰地查看当前目录下各个文件和文件夹的大小。
#-h 选项用于以人类可读的方式显示文件和目录的大小，以便更易于理解。
#* 表示通配符，用于匹配当前目录下的所有文件和文件夹。
# 执行该命令后，会列出当前目录下每个文件和文件夹的磁盘使用情况总结，并以人类可读的格式显示出来。

### Align and concat marker genes to create species tree #比对和合并标记基因，创建物种树


cat ~/data/Bacteria/HMM/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat ~/qyl/data/Burkholderia/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#378     381.5   779.5

cat ~/data/Bacteria/HMM/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat ~/qyl/data/Burkholderia/{}/replace.tsv |
            wc -l

    ' |
    paste - - |
    tsv-filter --invert --ge 2:360 --le 2:390 |
    cut -f 1 \
    > ~/qyl/data/Burkholderia/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process

rsync -avP \
    /mnt/c/shengxin/data/bac120.aln.sh \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/bac120.aln.sh

#bac120.aln.sh内容：
cat ~/data/Bacteria/HMM/bac120.tsv |  sed '1d' | cut -f 1 |
    grep -v -Fx -f ~/qyl/data/Burkholderia/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        cat ~/qyl/data/Burkholderia/{}/replace.tsv \
            > ~/qyl/data/Burkholderia/{}/{}.replace.tsv

        faops some ~/qyl/data/all/all.uniq.fa.gz <(
            cat ~/qyl/data/Burkholderia/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > ~/qyl/data/Burkholderia/{}/{}.pro.fa
    '

bsub -n 24 -q mpi -J "DddA-like" "bash ~/qyl/data/Burkholderia/bac120.aln.sh"

rsync -avP \
    /mnt/c/shengxin/data/tree.sh \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/tree.sh

bsub -n 24 -q mpi -J "DddA-like" "bash ~/qyl/data/Burkholderia/tree.sh"

#tree.sh内容：

# Align each markers with muscle
cat ~/data/Bacteria/HMM/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s ~/qyl/data/Burkholderia/{}/{}.pro.fa ]; then
            exit
        fi

        muscle -in ~/qyl/data/Burkholderia/{}/{}.pro.fa -out ~/qyl/data/Burkholderia/{}/{}.aln.fa
    '


for marker in $(cat ~/data/Bacteria/HMM/bac120.tsv | sed '1d' |cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s ~/qyl/data/Burkholderia/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s ~/qyl/data/Burkholderia/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat ~/qyl/data/Burkholderia/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s ~/qyl/data/Burkholderia/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > ~/qyl/data/Burkholderia/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat ~/data/Bacteria/HMM/bac120.tsv | sed '1d' | cut -f 1 ); do
    if [ ! -s ~/qyl/data/Burkholderia/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s ~/qyl/data/Burkholderia/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 ~/qyl/data/Burkholderia/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > ~/qyl/data/Burkholderia/bac120.aln.fas

cat ~/qyl/data/Burkholderia/2replace.tsv | cut -f 2 | sort | uniq |
    fasops concat ~/qyl/data/Burkholderia/bac120.aln.fas stdin -o ~/qyl/data/Burkholderia/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in ~/qyl/data/Burkholderia/bac120.aln.fa -out ~/qyl/data/Burkholderia/bac120.trim.fa -automated1

faops size ~/qyl/data/Burkholderia/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#24152
19021

# To make it faster
FastTree -fastest -noml ~/qyl/data/Burkholderia/bac120.trim.fa > ~/qyl/data/Burkholderia/bac120.trim.newick
bsub -n 24 -q mpi -J "tree" "FastTree -fastest -noml ~/qyl/data/Burkholderia/bac120.trim.fa > ~/qyl/data/Burkholderia/bac120.trim.newick"

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/bac120.trim.newick \
    /mnt/c/shengxin/data/Burkholderia 

#三个科的

rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/ASSEMBLY/Simiaoa_sunii/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/ASSEMBLY1

rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/ASSEMBLY/Roseburia_intestinalis/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/ASSEMBLY1

rsync -avP \
/mnt/c/shengxin/data/Protein/collect1.sh \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein1

bsub -n 24 -q mpi -J "DddA-like1" "bash ~/qyl/data/all/Protein1/collect1.sh"
bsub -n 24 -q mpi -J "DddA-like" "gzip -d all.replace.fa.gz"

cat ~/qyl/data/all/Protein1/all.replace.fa ~/qyl/data/all/all.replace.fa > ~/qyl/data/all/Protein1/All.replace.fa

bsub -n 24 -q mpi -J "DddA-like1" "E_VALUE=1e-10 ; cat ~/qyl/data/all/species.tsv |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 20 \"
            if [[ ! -d ~/qyl/data/all/ASSEMBLY1/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ~/qyl/data/all/ASSEMBLY1/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E \${E_VALUE} --domE \${E_VALUE} --noali --notextw ~/qyl/data/Burkholderia/HMM/DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \\\$1, {1}; '
        \" \
        > ~/qyl/data/all/Protein1/ALL-replace.tsv" 

rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/Protein/species.tsv \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/species0.tsv

bsub -n 24 -q mpi -J "DddA-like1" "E_VALUE=1e-10 ; cat ~/qyl/data/all/species0.tsv |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 20 \"
            if [[ ! -d ~/qyl/data/all/ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ~/qyl/data/all/ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E \${E_VALUE} --domE \${E_VALUE} --noali --notextw ~/qyl/data/Burkholderia/HMM/DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \\\$1, {1}; '
        \" \
        >> ~/qyl/data/all/Protein1/ALL-replace.tsv" #找所有三个科中有domain的，共2339个

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein1/ALL-replace.tsv \
    /mnt/c/shengxin/data/all

cat ~/qyl/data/all/Protein1/ALL-replace.tsv | tsv-select -f 2,1 | tr "\t" "_"> ~/qyl/data/all/Protein1/3.tsv
rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein1/3.tsv \
    /mnt/c/shengxin/data/all

rsync -avP \
    /mnt/c/shengxin/data/all/3.tsv \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein1/3.tsv  #在vscode中将3.tsv文件中的后缀.1去掉

bsub -n 24 -q mpi -J "DddA-like" "perl ~/qyl/data/all/fasplit.pl All.replace.fa input 3000000"

ls file{1..10}.fasta > list1
ls file{11..21}.fasta > list2
ls file{21..31}.fasta > list3
ls file{32..43}.fasta > list4

for i in {1..4}; do
  bsub -q largemem -n 24 -J "DddA-like" "bash job${i}.sh"
done

cat output/*.output > ~/qyl/data/all/Protein1/ALL-DddA-like.fa

faops filter -u ~/qyl/data/all/Protein1/ALL-DddA-like.fa ~/qyl/data/all/Protein1/All-unique.fa

bsub -n 24 -q mpi -J "DddA-like" "mafft --auto --thread 4 ~/qyl/data/all/Protein1/All-unique.fa > ~/qyl/data/all/Protein1/ALL-DddA-aln.fa"

FastTree ~/qyl/data/all/Protein1/ALL-DddA-aln.fa > ~/qyl/data/all/Protein1/ALL-DddA-aln.newick

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/all/Protein1/ALL-DddA-aln.newick \
    /mnt/c/shengxin/data/all
