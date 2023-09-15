### List all ranks 列出所有等级

```shell
nwr member Actinospica |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'
```
| rank | count |
| --- | ---:|
| genus | 1 |
| species | 4 |
| no rank | 1 |
| strain | 1 |
```shell
nwr lineage Actinospica |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tBacteria/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

 #查看丛生放线菌属属的分类上的基本情况   
```
| #rank | sci_name | tax_id |
| --- | --- | --- |
| superkingdom | Bacteria | 2 |
| phylum | Actinomycetota | 201174 |
| class | Actinomycetes | 1760 |
| order | Catenulisporales | 414714 |
| family | Actinospicaceae | 414713 |
| *genus* | Actinospica | 414715 |
### Species with assemblies 具有组装的物种

Actinospicaceae是细菌界中的一个科

```shell
cd /mnt/c/shengxin
mkdir -p data/Actinospica/summary
cd /mnt/c/shengxin/data/Actinospica/summary

#找出和丛生放线菌属同一科(Actinospicaceae)的所有属

# should have a valid name of genus
nwr member Actinospicaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
# 2 genus.list.tsv

#丛生放线菌属同属的所有参考物种基因组信息
cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %' -- Crossbreeding of two species
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

#丛生放线菌属同属的所有物种信息(genbank)
cat genus.list.tsv | cut -f 1 |
while read RANK_ID; do
    echo "
        SELECT
            species_id,
            species,
            COUNT(*) AS count
        FROM ar
        WHERE 1=1
            AND genus_id = ${RANK_ID}
            AND species NOT LIKE '% sp.%'
            AND species NOT LIKE '% x %'
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
done |
    tsv-sort -k2,2 \
    > GB1.tsv

wc -l RS*.tsv GB*.tsv
# 4 RS1.tsv
  5 GB1.tsv
  9 total

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     4
#GB1     9

```

## Download all assemblies

### Create .assembly.tsv

```shell
cd /mnt/c/shengxin/data/Actinospica/summary

# Reference genome
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharothrix')
            " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv
  #菌属的参考菌株的基因组信息

# RS
SPECIES=$(
    cat RS1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        assembly_accession
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv
 # 丛生放线菌属同科的各属的参考菌株基因组信息

# Preference for refseq
cat raw.tsv |
    tsv-select -H -f "assembly_accession" \
    > rs.acc.tsv

# GB
SPECIES=$(
    cat GB1.tsv |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        species || ' ' || infraspecific_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, biosample, assembly_level,
        gbrs_paired_asm
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND genome_rep IN ('Full')
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite |
    tsv-join -f rs.acc.tsv -k 1 -d 7 -e \
    >> raw.tsv
 #丛生放线菌属同科的各属的菌株基因组信息

cat raw.tsv |
    tsv-uniq |
    datamash check
#10 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-6 |
    perl ../../abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\tbiosample\tspecies\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[6]}++; # abbr_name
        $seen{$F[6]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\t%s\n}, $F[6], $F[3], $F[4], $F[1], $F[5];
        ' |
    tsv-filter --or --str-in-fld 2:ftp --str-in-fld 2:http |
    keep-header -- tsv-sort -k4,4 -k1,1 \
    > Actinospica.assembly.tsv
 #创建简写名称的丛生放线菌属属水平的各菌株基因组下载文件

 datamash check < Actinospica.assembly.tsv
 #10 lines, 5 fields

# find potential duplicate strains or assemblies 
cat Actinospica.assembly.tsv |
    tsv-uniq -f 1 --repeated
#检查有没有重复

cat Actinospica.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp
    # 检查下载链接是否正确

# Save the file to another directory to prevent accidentally changing it
cp Actinospica.assembly.tsv /mnt/c/shengxin/data/Actinospica/assembly/

# Cleaning
rm raw*.*sv
```
### Count before download 下载前计数

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class
Strains.taxon.tsv - 分类信息：物种、属、科、目和类

```shell
cd /mnt/c/shengxin/data/Actinospica

nwr template summary/Actinospica.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh #strains.taxon.tsv共6列，是为了统计各个菌株的数量以及所在的物种，属，科，目，纲的数量

#taxa.tsv内容
item    count
strain  9
species 5
genus   2
family  1
order   1
class   1

# genus.lst and genus.count.tsv
bash Count/rank.sh

#genus.lst 文件只有一列，提取了属的名称
#genus.count.tsv 统计丛生放线菌属属以及各属的species和strains去重后的数量

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'
```
| genus | #species | #strains |
|---|--:|--:|
| Actinocrinis | 1 | 1 |
| Actinospica | 4 | 8 |


### Download and check 下载并检查
```shell
cd /mnt/c/shengxin/data/Actinospica
nwr template ../Actinospica/assembly/Actinospica.assembly.tsv\
    --ass

# --ass: ASSEMBLY/
#     * One TSV file
#         * url.tsv #三列：菌株简写名称；下载网址；物种名
#     * And five Bash scripts
#         * rsync.sh
#         * check.sh
#         * n50.sh [LEN_N50] [N_CONTIG] [LEN_SUM]

# Run
bash ASSEMBLY/rsync.sh


# Check md5; create check.lst
# rm ASSEMBLY/check.lst
bash ASSEMBLY/check.sh

# N50 C S; create n50.tsv and n50.pass.tsv
bash ASSEMBLY/n50.sh 100000 1000 1000000
# 10000：最小n50长度
# 1000： contig数量小于1000
# 1000000: 总长度大于1000000.

# Adjust parameters passed to `n50.sh`#调整传递给' n50.sh '的参数
cat ASSEMBLY/n50.tsv |
    tsv-filter -H --str-in-fld "name:_GCF_" | 
    tsv-summarize -H --min "N50,S" --max "C"

#N50_min S_min   C_max
#26256   6765010 944

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0,0.5" --quantile "N50:0,0.5"  --quantile "C:0.5,0.9"
    #计算数据的分位数
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#6765010 7519806 11223   45308   430     1014.4

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh
#收集已经下载的基因组的详细信息

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

| #item | fields | lines |
|---|--:|--:|
| url.tsv | 3 | 9 |
| check.lst | 1 | 9 |
| collect.tsv | 20 | 10 |
| n50.tsv | 4 | 10 |
| n50.pass.tsv | 4 | 4 |
| collect.pass.tsv | 23 | 4 |
| pass.lst | 1 | 3 |
| omit.lst | 1 | 1 |
| rep.lst | 1 | 2 |

# 统计下载的文件行数：
# collect.pass.tsv ：下载的基因组通过了n50的检验 -5308行（有一行标题）；pass.lst 5307行（不含标题），取collect.pass.tsv的第一列名称。
# omit.lst ：没有蛋白序列信息，文件只有一列，基因组名 225个基因组
# rep.lst : 含有参考序列的基因组 268个基因组
```
### Rsync to hpcc 把本地文件与超算同步

```bash
rsync -avP \
    /mnt/c/shengxin/data/Actinospica/assembly/ \
    wangq@20219.37.251:qyl/data/Actinospica/ASSEMBLY
#本地运行

rsync -avP \
    /mnt/c/shengxin/data/Actinospica/summary/ \
    wangq@20219.37.251:qyl/data/Actinospica/summary
```

## BioSample 生物样本
```shell
cd /mnt/c/shengxin/data/Actinospica

ulimit -n `ulimit -Hn`

cp ASSEMBLY/Actinospica.assembly.tsv summary/
nwr template /mnt/c/shengxin/data/Actinospica/assembly/Actinospica.assembly.tsv \
    --bs

head BioSample/sample.tsv
#SAMN02597168    Actinos_rob_DSM_44927_GCF_000504285_1   Actinospica_robiniae
SAMN02645325    Actinos_acidip_NRRL_B_24431_GCA_000717315_1     Actinospica_acidiphila

bash BioSample/download.sh
#在本地下载

# Ignore rare attributes#忽略稀少属性
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#6473 lines, 107 fields

cp Biosample/attributes.lst summary/ 
cp Biosample/biosample.tsv summary/


rsync -avP \
    /mnt/c/shengxin/data/Actinospica/Biosample/ \
    wangq@20219.37.251:qyl/data/Actinospica/Biosample
```

## MinHash
```shell
cd /mnt/c/shengxin/data/Actinospica

nwr template /mnt/c/shengxin/data/Actinospica/summary/Actinospica.assembly.tsv \
    --mh \
    --parallel 16 \
    --in ASSEMBLY/pass.lst \
    --ani-ab 0.05 \
    --ani-nr 0.005 \
    --height 0.4 

# Compute assembly sketches
bash MinHash/compute.sh

# Distances within species
bash MinHash/species.sh

# Abnormal strains
bash MinHash/abnormal.sh

cat MinHash/abnormal.lst


cat MinHash/abnormal.lst |wc -l


# Non-redundant strains within species
bash MinHash/nr.sh

# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh

rsync -avP \
    wangq@20219.37.251:/share/home/wangq/qyl/data/Actinospica/Biosample/ \
    /mnt/c/shengxin/data/Actinospica/Biosample

rsync -avP \
    wangq@20219.37.251:/share/home/wangq/qyl/data/Actinospica/MinHash/ \
    /mnt/c/shengxin/data/Actinospica/MinHash
    #本地运行   

```

### Condense branches in the minhash tree
```shell
mkdir -p cd /mnt/c/shengxin/data/Actinospica/tree
cd /mnt/c/shengxin/data/Actinospica/tree

nw_reroot ../MinHash/tree.nwk S_vio |
    nw_order -cn - \
    > minhash.reroot.newick

# rank::col
ARRAY=(
#    'order::5'
#    'family::4'
#    'genus::3'
    'species::2'
)

rm minhash.condensed.map
CUR_TREE=minhash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ../condense_tree.sh ${CUR_TREE} ../Count/strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick minhash.${GROUP_NAME}.newick
    cat condense.map >> minhash.condensed.map

    CUR_TREE=minhash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 minhash.species.newick |
    rsvg-convert -o Actinospica.minhash.png
#报错：超过了 32767 像素，这是 rsvg-convert 当前无法处理的最大尺寸限制。
nw_display -s -b 'visibility:hidden' -w 800 -v 15 minhash.species.newick | rsvg-convert -o Actinospica.minhash.png
#还是过大

nw_display -s -b 'visibility:hidden' -w 600 -v 10 minhash.species.newick |
    rsvg-convert -o Actinospica.minhash.png

```
## Count valid species and strains #计算有效物种数和菌株数

### For *genomic alignments* 用于*基因组比对*

```shell
cd /mnt/c/shengxin/data/Actinospica
nwr template /mnt/c/shengxin/data/Actinospica/assembly/Actinospica.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --rank genus \
    --lineage family --lineage genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# Can accept N_COUNT
bash Count/lineage.sh

cat Count/lineage.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/genome.taxon.tsv
```
### For *protein families* #用于*蛋白质家族*

```shell
cd /mnt/c/shengxin/data/Actinospica/

nwr template /mnt/c/shengxin/data/Actinospica/assembly/Actinospica.assembly.tsv \
    --count \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst \
    --rank genus

# strains.taxon.tsv
bash Count/strains.sh

# .lst and .count.tsv
bash Count/rank.sh

cat Count/genus.count.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

# copy to summary/
cp Count/strains.taxon.tsv summary/protein.taxon.tsv

```


## Collect Protein

```shell
cd /mnt/c/shengxin/data/Actinospica/

nwr template /mnt/c/shengxin/data/Actinospica/assembly/Actinospica.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# * --pro: Protein/
#     * One TSV file
#         * species.tsv
#     * collect.sh

# collect Protein
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat
```

## Phylogenetics with bac120 #系统发育

### Find corresponding Protein by `hmmsearch` #通过`hmmsearch`查找相应的蛋白质

```shell
cd /mnt/c/shengxin/data/Actinospica/
mkdir -p HMM
E_VALUE=1e-20

# Find all genes
    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ../DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/replace.tsv

```
# 建立蛋白树
```shell
cd /mnt/c/shengxin/data/Actinospica/
cat Protein/replace.tsv | tsv-select -f 2,1 > Protein/3.tsv
faops some Protein/all.replace.fa.gz <(tsv-select -f 1 Protein/3.tsv) Protein/DddA-like.fa
#Protein/3.tsv格式是B_pseudoma_UMC107_GCF_002921075_1_WP_004533223，
#zcat Protein/all.replace.fa.gz |head 显示格式是>B_aenigmatica_AU17325_GCF_002223275_1_WP_000522996
所以要将Protein/3.tsv中变成B_pseudoma_UMC107_GCF_002921075_1_WP_004533223，可以直接在vscode中进行替换 #要转义.1，是\.1 

muscle -in Protein/DddA-like.fa -out Protein/DddA-like.aln.fa

FastTree Protein/DddA-like.aln.fa > Protein/DddA-like.aln.newick

nw_reroot Protein/DddA-like.aln.newick $(nw_labels Protein/DddA-like.aln.newick | grep -E "S_vio") |
    nw_order -c n - \
    > Tree/DddA-like.reoot.newick

nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/DddA-like.reoot.newick | rsvg-convert -o tree/DddA-like.reoot.png
```
## Phylogenetics with bac120 #系统发育

### Find corresponding Protein by `hmmsearch` #通过`hmmsearch`查找相应的蛋白质
```shell
#TIGRFAM
mkdir -p mnt/shengxin/data/Actinospica/HMM/TIGRFAM
cd mnt/shengxin/data/Actinospica/HMM/TIGRFAM
wget -N --content-disposition ftp://ftp.jcvi.org/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz

mkdir -p HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR02013.HMM
#tar：tar 命令用于操作 tar 归档文件。
#--directory HMM：指定目标目录为 HMM，在解压缩时将文件提取到该目录中。
#-x：解压缩模式，表示从归档文件中提取文件。
#-z：使用 gzip 解压缩算法进行解压缩。
#-v：显示详细的解压缩过程。
#-f TIGRFAMs_14.0_HMM.tar.gz：指定要解压缩的归档文件为 TIGRFAMs_14.0_HMM.tar.gz。
#TIGR02013.HMM：指定要提取的文件名
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR00485.HMM

mkdir -p /mnt/c/shengxin/data/Actinospica/HMM/bac120
cd /mnt/c/shengxin/data/Actinospica/HMM/bac120


#下载bac120.tsv文件
mnt/shengxin/data/Actinospica/HMM/bac120/

mkdir -p HMM

cat /mnt/c/shengxin/data/Actinospica/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory HMM -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat /mnt/c/shengxin/data/Actinospica/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L http://pfam.xfam.org/family/{}/hmm > HMM/{}.HMM
    '
#bac120进行系统发育
E_VALUE=1e-20

# Find all genes
cd /mnt/c/shengxin/data/Actinospica/
cat Protein/replace.tsv |cut -f 2 >Protein/1.tsv
cat Protein/species.tsv |tsv-join -f Protein/1.tsv -k 1 -d 1 >Protein/2.tsv
cat Protein/2.tsv | tsv-select -f 2,1 > temp.tsv

for marker in $(cat HMM/bac120/bac120.tsv | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

mkdir -p Protein/${marker}

cat temp.tsv |
        parallel --colsep '\t' --no-run-if-empty  --linebuffer -k -j 8 "
            gzip -dcf ASSEMBLY/{1}/{2}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/bac120/HMM/${marker}.HMM - |
                grep '>>' |
                perl -nl -e ' m{>>\s+(\S+)} and printf qq{%s\t%s\n}, \$1, {2}; '
        " \
        > Protein/${marker}/replace.tsv

    >&2 echo

done

```
### Align and concat marker genes to create species tree #比对和合并标记基因，创建物种树
```shell
cd /mnt/c/shengxin/data/Actinospica/

cat HMM/bac120/bac120.tsv | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75
#1414    1415    1815.25

cat HMM/bac120/bac120.tsv | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l

    ' |
    paste - - |
    tsv-filter --invert --ge 2:1300 --le 2:1500 |
    cut -f 1 \
    > Protein/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/bac120/bac120.tsv | cut -f 1 |
    grep -v -Fx -f Protein/bac120.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        cat Protein/{}/replace.tsv \
            > Protein/{}/{}.replace.tsv

        faops some Protein/all.uniq.fa.gz <(
            cat Protein/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > Protein/{}/{}.pro.fa
    '

# Align each markers with muscle
cat HMM/bac120/bac120.tsv | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '


for marker in $(cat HMM/bac120/bac120.tsv |cut -f 1); do
    >&2 echo "==> marker [${marker}]"
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi

    # sometimes `muscle` can not produce alignments
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # 1 name to many names
    cat Protein/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s Protein/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > Protein/${marker}/${marker}.replace.fa
done

# Concat marker genes
for marker in $(cat HMM/bac120/bac120.tsv |cut -f 1); do
    if [ ! -s Protein/${marker}/${marker}.pro.fa ]; then
        continue
    fi
    if [ ! -s Protein/${marker}/${marker}.aln.fa ]; then
        continue
    fi

    # sequences in one line
    faops filter -l 0 Protein/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > Protein/bac120.aln.fas

cat Protein/replace.tsv | cut -f 2 | sort | uniq |
    fasops concat Protein/bac120.aln.fas stdin -o Protein/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/bac120.aln.fa -out Protein/bac120.trim.fa -automated1

faops size Protein/bac120.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#31965
#27534

# To make it faster
FastTree -fastest -noml Protein/bac120.trim.fa > Protein/bac120.trim.newick

nw_reroot Protein/bac120.trim.newick S_vio |
    nw_order -c n - \
    > tree/bac120.reroot.newick

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/bac120.reroot.newick |
    rsvg-convert -o tree/Actinospica.marker.png

rsync -avP \
    /mnt/c/shengxin/data/Actinospica/HMM/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Actinospica/HMM
```