### List all ranks 列出所有等级

```shell
nwr member Streptomyces |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'
```
| rank             | count |
| ---------------- | ----: |
| genus            |     1 |
| species          |  1127 |
| species group    |    15 |
| strain           |    91 |
| subspecies       |   111 |
| no rank          |    37 |
| species subgroup |    10 |
```shell
nwr lineage Streptomyces |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tBacteria/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

 #查看链霉属的分类上的基本情况   
```
| #rank        | sci_name          | tax_id |
| ------------ | ----------------- | ------ |
| superkingdom | Bacteria          | 2      |
| phylum       | Actinomycetota    | 201174 |
| class        | Actinomycetes     | 1760   |
| order        | Kitasatosporales  | 85011  |
| family       | Streptomycetaceae | 2062   |
| *genus*      | Streptomyces      | 1883   |

### Species with assemblies 具有组装的物种

Streptomycetaceae 是细菌界中的一个科

```shell
cd /mnt/c/shengxin
mkdir -p data/Streptomyces/summary
cd /mnt/c/shengxin/data/Streptomyces/summary

#找出和链霉同一科(Streptomycetaceae)的所有属

# should have a valid name of genus
nwr member Streptomycetaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
# 16 genus.list.tsv

#链霉同属的所有参考物种基因组信息
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
            AND species NOT LIKE '% x %' -- Crossbreeding of two species
            AND genome_rep IN ('Full')
        GROUP BY species_id
        HAVING count >= 1
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
done |
    tsv-sort -k2,2 \
    > RS1.tsv

#链霉同属的所有物种信息(genbank)
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
#  2083 RS1.tsv
  2282 GB1.tsv
  4365 total

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     2921
#GB1     3320

```

## Download all assemblies

### Create .assembly.tsv

```shell
cd /mnt/c/shengxin/data/Streptomyces/summary

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
 # 链霉同科的各属的参考菌株基因组信息

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
 #链霉同科的各属的菌株基因组信息

cat raw.tsv |
    tsv-uniq |
    datamash check
#3321 lines, 7 fields

# Create abbr.
cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    tsv-select -f 1-6 |
    perl abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
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
    > Streptomyces.assembly.tsv
 #创建简写名称的链霉属水平的各菌株基因组下载文件

 datamash check < Streptomyces.assembly.tsv
 #3321 lines, 5 fields

# find potential duplicate strains or assemblies 
cat Streptomyces.assembly.tsv |
    tsv-uniq -f 1 --repeated
#检查有没有重复

cat Streptomyces.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp
    # 检查下载链接是否正确

# Save the file to another directory to prevent accidentally changing it
mkdir -p ../assembly/
cp Streptomyces.assembly.tsv /mnt/c/shengxin/data/Streptomyces/assembly/

# Cleaning
rm raw*.*sv
```
### Count before download 下载前计数

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class
Strains.taxon.tsv - 分类信息：物种、属、科、目和类

```shell
cd /mnt/c/shengxin/data/Streptomyces

nwr template summary/Streptomyces.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh #strains.taxon.tsv共6列，是为了统计各个菌株的数量以及所在的物种，属，科，目，纲的数量

#taxa.tsv内容
item    count
strain  1687
species 700
genus   12
family  1
order   1
class   1

#strains.taxon.tsv内容
Ac_acidid_KK5PA1_GCF_016918855_1        Actinacidiphila acididurans     Actinacidiphila Streptomycetaceae       Kitasatosporales        Actinomycetes
Ac_alni_CGMCC_4_3510_GCF_900112845_1    Actinacidiphila alni    Actinacidiphila Streptomycetaceae       KitasatosporaleActinomycetes

# genus.lst and genus.count.tsv
bash Count/rank.sh

#genus.lst 文件只有一列，提取了属的名称
#genus.count.tsv 统计链霉属以及各属的species和strains去重后的数量

mv Count/genus.count.tsv Count/genus.before.tsv

cat Count/genus.before.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'
```
| genus              | #species | #strains |
| ------------------ | -------: | -------: |
| Actinacidiphila    |       12 |       17 |
| Allostreptomyces   |        1 |        1 |
| Embleya            |        2 |        3 |
| Kitasatospora      |       23 |       41 |
| Mangrovactinospora |        1 |        1 |
| Peterkaempfera     |        2 |        2 |
| Phaeacidiphilus    |        1 |        1 |
| Streptacidiphilus  |       10 |       12 |
| Streptantibioticus |        1 |        3 |
| Streptomyces       |      644 |     1602 |
| Wenjunlia          |        2 |        3 |
| Yinghuangia        |        1 |        1 |


### Download and check 下载并检查
```shell
cd /mnt/c/shengxin/data/Streptomyces
nwr template ../Streptomyces/assembly/Streptomyces.assembly.tsv\
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

#这步速度太慢了，复制超算中ASSEMBLY的部分文件到我的超算中
cd ~/qyl/data/
mkdir -p Streptomyces/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Streptantibioticus_* ~/qyl/data/Streptomyces/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Streptomyces_* ~/qyl/data/Streptomyces/ASSEMBLY

#下载到本地
rsync -avP \
wangq@202.119.37.251:qyl/data/Streptomyces/ASSEMBLY \
/mnt/c/shengxin/data/Streptomyces


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
#5805    1550368 1888

cat ASSEMBLY/n50.tsv |
    tsv-summarize -H --quantile "S:0,0.5" --quantile "N50:0,0.5"  --quantile "C:0.5,0.9"
    #计算数据的分位数
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#1027029 7123830 1203    201283  96      234

# Collect; create collect.tsv
bash ASSEMBLY/collect.sh
#收集已经下载的基因组的详细信息

# After all completed
bash ASSEMBLY/finish.sh

cp ASSEMBLY/collect.pass.tsv summary/

cat ASSEMBLY/counts.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'

| #item            | fields | lines |
| ---------------- | -----: | ----: |
| url.tsv          |      3 | 6,475 |
| check.lst        |      1 | 6,475 |
| collect.tsv      |     20 | 6,476 |
| n50.tsv          |      4 | 6,476 |
| n50.pass.tsv     |      4 | 5,308 |
| collect.pass.tsv |     23 | 5,308 |
| pass.lst         |      1 | 5,307 |
| omit.lst         |      1 |   225 |
| rep.lst          |      1 |   268 |

# 统计下载的文件行数：
# collect.pass.tsv ：下载的基因组通过了n50的检验 -5308行（有一行标题）；pass.lst 5307行（不含标题），取collect.pass.tsv的第一列名称。
# omit.lst ：没有蛋白序列信息，文件只有一列，基因组名 225个基因组
# rep.lst : 含有参考序列的基因组 268个基因组
```
### Rsync to hpcc 把本地文件与超算同步

```bash
rsync -avP \
    /mnt/c/shengxin/data/Streptomyces/Count/ \
    wangq@202.119.37.251:qyl/data/Streptomyces/Count
#本地运行

rsync -avP \
    /mnt/c/shengxin/data/Streptomyces/summary/ \
    wangq@202.119.37.251:qyl/data/Streptomyces/summary
```

## BioSample 生物样本
```shell
cd /mnt/c/shengxin/data/Streptomyces

ulimit -n `ulimit -Hn`

cp ASSEMBLY/Streptomyces.assembly.tsv summary/
nwr template /mnt/c/shengxin/data/Streptomyces/assembly/Streptomyces.assembly.tsv \
    --bs

head BioSample/sample.tsv
#SAMD00000356    Parab_fer_NBRC_106233_GCF_000685035_1   ParaStreptomyces_ferrariae
#SAMD00000357    Parab_mim_NBRC_106338_GCF_000739815_1   ParaStreptomyces_mimosarum

bash BioSample/download.sh
#在本地下载

# Ignore rare attributes#忽略稀少属性
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#6473 lines, 107 fields

cp Biosample/attributes.lst summary/ 
cp Biosample/biosample.tsv summary/

#biosample.tsv 文件内容
#name   BioSample       sample name     observed biotic relationship    collection date environmental medium    geographic location     isolation and growth condition  latitude and longitude  locus_tag_prefix        number of replicons    project name     reference for biomaterial       strain  trophic level   source material identifiers     isolation sourchost     External Id     Submitter Id    culture collection      host health state       anonymized name common name    serovar  subject id      supplier_name   investigation type      sequencing method       assembly quality        assembly software       binning parameters      binning software        completeness score      completeness software   contamination score     metagenomic source      sample derived from     taxonomic identity marker       scientific_name local environmental context     broker name     isolate Alias   SRA accession   Title   alias   strain_synonym  substrain      anonymized_name  relationship to oxygen  collected by    passage history GOLD Stamp ID   Gram Staining   disease environment     Cell Shape      Motility        Temperature Range       Sporulation     Phenotypes      Isolation Site  Temperature Optimum     host disease    sample type     environmental package   biomaterial provider    depth   alternate strain name   organism modifier note  host age        note    host sex        host disease outcome    genotype        host tissue sampled     alternate_ID    elevation       identified by   host description        serotype        pathotype      subgroup subtype specimen_category       provider        collection method       uFORGE_Sample_ID        pH      altitudhost disease stage       host subject id isolate name alias      identification method   biovar  phylotype       sample size     metagenomic     derived from    host taxonomy ID        source type     Phylotype       temperature     sample storage location sample storage temperature      host infra specific name
Parab_fer_NBRC_106233_GCF_000685035_1   SAMD00000356    NBRC 106233     free living             Iron Mine       Brazil:Minas Gerais     17012573                BFE01S          Streptomyces ferrariae NBRC 106233 genome sequencing project   NBRC 106233      heterotroph     NBRC 106233

rsync -avP \
    /mnt/c/shengxin/data/Streptomyces/Biosample/ \
    wangq@202.119.37.251:qyl/data/Streptomyces/Biosample
```

## MinHash
```shell
cd /mnt/c/shengxin/data/Streptomyces

nwr template /mnt/c/shengxin/data/Streptomyces/summary/Streptomyces.assembly.tsv \
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
#B_amb_RZ2MS16_GCF_001443045_1
#B_anth_3723STDY6437370_GCA_900631665_1

cat MinHash/abnormal.lst |wc -l
493

# Non-redundant strains within species
bash MinHash/nr.sh

# Distances between all selected sketches, then hierarchical clustering
bash MinHash/dist.sh

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Streptomyces/Biosample/ \
    /mnt/c/shengxin/data/Streptomyces/Biosample

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Streptomyces/MinHash/ \
    /mnt/c/shengxin/data/Streptomyces/MinHash
    #本地运行   

```

### Condense branches in the minhash tree
```shell
mkdir -p cd /mnt/c/shengxin/data/Streptomyces/tree
cd /mnt/c/shengxin/data/Streptomyces/tree

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
    rsvg-convert -o Streptomyces.minhash.png
#报错：超过了 32767 像素，这是 rsvg-convert 当前无法处理的最大尺寸限制。
nw_display -s -b 'visibility:hidden' -w 800 -v 15 minhash.species.newick | rsvg-convert -o Streptomyces.minhash.png
#还是过大

nw_display -s -b 'visibility:hidden' -w 600 -v 10 minhash.species.newick |
    rsvg-convert -o Streptomyces.minhash.png

```
## Count valid species and strains #计算有效物种数和菌株数

### For *genomic alignments* 用于*基因组比对*

```shell
cd /mnt/c/shengxin/data/Streptomyces
nwr template /mnt/c/shengxin/data/Streptomyces/assembly/Streptomyces.assembly.tsv \
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
| genus               | #species | #strains |
| ------------------- | -------: | -------: |
| Streptomyces        |       41 |     3969 |
| Caballeronia        |       23 |       36 |
| Chitinasiproducens  |        1 |        1 |
| Chitinimonas        |        3 |        4 |
| Cupriavidus         |       21 |      117 |
| Ephemeroptericola   |        1 |        1 |
| Formosimonas        |        1 |        1 |
| Lautropia           |        2 |        5 |
| Limnobacter         |        3 |        3 |
| Mycetohabitans      |        2 |        2 |
| Mycoavidus          |        1 |        1 |
| Pandoraea           |       28 |       64 |
| ParaStreptomyces    |       83 |      218 |
| Pararobbsia         |        2 |        2 |
| Polynucleobacter    |       29 |      157 |
| Quisquiliibacterium |        1 |        1 |
| Ralstonia           |       11 |      211 |
| Robbsia             |        1 |        4 |
| Saccharothrix       |        1 |        1 |
| Trinickia           |        7 |       14 |
| Zeimonas            |        2 |        2 |


| genus               | #species | #strains |
| ------------------- | -------: | -------: |
| Streptomyces        |       41 |     3969 |
| Caballeronia        |       23 |       36 |
| Chitinasiproducens  |        1 |        1 |
| Chitinimonas        |        3 |        4 |
| Cupriavidus         |       21 |      117 |
| Ephemeroptericola   |        1 |        1 |
| Formosimonas        |        1 |        1 |
| Lautropia           |        2 |        5 |
| Limnobacter         |        3 |        3 |
| Mycetohabitans      |        2 |        2 |
| Mycoavidus          |        1 |        1 |
| Pandoraea           |       28 |       64 |
| ParaStreptomyces    |       83 |      218 |
| Pararobbsia         |        2 |        2 |
| Polynucleobacter    |       29 |      157 |
| Quisquiliibacterium |        1 |        1 |
| Ralstonia           |       11 |      211 |
| Robbsia             |        1 |        4 |
| Saccharothrix       |        1 |        1 |
| Trinickia           |        7 |       14 |
| Zeimonas            |        2 |        2 |
qin@Qin:/mnt/c/shengxin/data/Streptomyces$ bash Count/lineage.sh
==> Count/lineage.sh <==
==> Done.
qin@Qin:/mnt/c/shengxin/data/Streptomyces$ cat Count/lineage.count.tsv |
>     mlr --itsv --omd cat |
>     perl -nl -e 's/-\s*\|$/-:|/; print'
| #family            | genus               | species                             | count |
| ------------------ | ------------------- | ----------------------------------- | ----: |
| Streptomycetaceae  | Streptomyces        | Streptomyces aenigmatica            |    20 |
|                    |                     | Streptomyces ambifaria              |    89 |
|                    |                     | Streptomyces anthina                |    27 |
|                    |                     | Streptomyces arboris                |     6 |
|                    |                     | Streptomyces catarinensis           |     1 |
|                    |                     | Streptomyces cenocepacia            |   395 |
|                    |                     | Streptomyces cepacia                |   234 |
|                    |                     | Streptomyces contaminans            |    86 |
|                    |                     | Streptomyces diffusa                |    16 |
|                    |                     | Streptomyces dolosa                 |    20 |
|                    |                     | Streptomyces gladioli               |   260 |
|                    |                     | Streptomyces glumae                 |    63 |
|                    |                     | Streptomyces guangdongensis         |     1 |
|                    |                     | Streptomyces humptydooensis         |     3 |
|                    |                     | Streptomyces lata                   |    16 |
|                    |                     | Streptomyces latens                 |     6 |
|                    |                     | Streptomyces mallei                 |    58 |
|                    |                     | Streptomyces mayonis                |     2 |
|                    |                     | Streptomyces metallica              |     7 |
|                    |                     | Streptomyces multivorans            |   490 |
|                    |                     | Streptomyces oklahomensis           |     9 |
|                    |                     | Streptomyces orbicola               |     3 |
|                    |                     | Streptomyces paludis                |     3 |
|                    |                     | Streptomyces perseverans            |     1 |
|                    |                     | Streptomyces plantarii              |     5 |
|                    |                     | Streptomyces pseudomallei           |  1504 |
|                    |                     | Streptomyces pseudomultivorans      |    10 |
|                    |                     | Streptomyces puraquae               |     2 |
|                    |                     | Streptomyces pyrrocinia             |     3 |
|                    |                     | Streptomyces reimsis                |     1 |
|                    |                     | Streptomyces savannae               |     4 |
|                    |                     | Streptomyces semiarida              |     4 |
|                    |                     | Streptomyces seminalis              |    15 |
|                    |                     | Streptomyces singularis             |     1 |
|                    |                     | Streptomyces sola                   |     1 |
|                    |                     | Streptomyces stabilis               |     2 |
|                    |                     | Streptomyces stagnalis              |   101 |
|                    |                     | Streptomyces territorii             |    37 |
|                    |                     | Streptomyces thailandensis          |    27 |
|                    |                     | Streptomyces ubonensis              |   297 |
|                    |                     | Streptomyces vietnamiensis          |   139 |
|                    | Caballeronia        | Caballeronia arationis              |     1 |
|                    |                     | Caballeronia arvi                   |     1 |
|                    |                     | Caballeronia calidae                |     1 |
|                    |                     | Caballeronia catudaia               |     1 |
|                    |                     | Caballeronia concitans              |     1 |
|                    |                     | Caballeronia cordobensis            |     2 |
|                    |                     | Caballeronia fortuita               |     1 |
|                    |                     | Caballeronia glathei                |     2 |
|                    |                     | Caballeronia glebae                 |     1 |
|                    |                     | Caballeronia grimmiae               |     3 |
|                    |                     | Caballeronia hypogeia               |     1 |
|                    |                     | Caballeronia insecticola            |     1 |
|                    |                     | Caballeronia jiangsuensis           |     1 |
|                    |                     | Caballeronia novacaledonica         |     3 |
|                    |                     | Caballeronia pedi                   |     1 |
|                    |                     | Caballeronia peredens               |     1 |
|                    |                     | Caballeronia ptereochthonis         |     1 |
|                    |                     | Caballeronia sordidicola            |     3 |
|                    |                     | Caballeronia telluris               |     1 |
|                    |                     | Caballeronia temeraria              |     1 |
|                    |                     | Caballeronia terrestris             |     1 |
|                    |                     | Caballeronia turbans                |     1 |
|                    |                     | Caballeronia zhejiangensis          |     6 |
|                    | Chitinasiproducens  | Chitinasiproducens palmae           |     1 |
|                    | Chitinimonas        | Chitinimonas arctica                |     1 |
|                    |                     | Chitinimonas koreensis              |     2 |
|                    |                     | Chitinimonas taiwanensis            |     1 |
|                    | Cupriavidus         | Cupriavidus agavae                  |     1 |
|                    |                     | Cupriavidus alkaliphilus            |     7 |
|                    |                     | Cupriavidus basilensis              |     5 |
|                    |                     | Cupriavidus campinensis             |     3 |
|                    |                     | Cupriavidus cauae                   |     2 |
|                    |                     | Cupriavidus gilardii                |    17 |
|                    |                     | Cupriavidus laharis                 |     1 |
|                    |                     | Cupriavidus malaysiensis            |     1 |
|                    |                     | Cupriavidus metallidurans           |    15 |
|                    |                     | Cupriavidus nantongensis            |     2 |
|                    |                     | Cupriavidus necator                 |     5 |
|                    |                     | Cupriavidus neocaledonicus          |     3 |
|                    |                     | Cupriavidus numazuensis             |     1 |
|                    |                     | Cupriavidus oxalaticus              |     3 |
|                    |                     | Cupriavidus pampae                  |     1 |
|                    |                     | Cupriavidus pauculus                |     2 |
|                    |                     | Cupriavidus pinatubonensis          |     2 |
|                    |                     | Cupriavidus plantarum               |     6 |
|                    |                     | Cupriavidus respiraculi             |     2 |
|                    |                     | Cupriavidus taiwanensis             |    37 |
|                    |                     | Cupriavidus yeoncheonensis          |     1 |
|                    | Ephemeroptericola   | Ephemeroptericola cinctiostellae    |     1 |
|                    | Formosimonas        | Formosimonas limnophila             |     1 |
|                    | Lautropia           | Lautropia dentalis                  |     1 |
|                    |                     | Lautropia mirabilis                 |     4 |
|                    | Limnobacter         | Limnobacter alexandrii              |     1 |
|                    |                     | Limnobacter humi                    |     1 |
|                    |                     | Limnobacter thiooxidans             |     1 |
|                    | Mycetohabitans      | Mycetohabitans endofungorum         |     1 |
|                    |                     | Mycetohabitans rhizoxinica          |     1 |
|                    | Mycoavidus          | Mycoavidus cysteinexigens           |     1 |
|                    | Pandoraea           | Pandoraea anapnoica                 |     1 |
|                    |                     | Pandoraea anhela                    |     1 |
|                    |                     | Pandoraea apista                    |    18 |
|                    |                     | Pandoraea aquatica                  |     1 |
|                    |                     | Pandoraea bronchicola               |     1 |
|                    |                     | Pandoraea capi                      |     1 |
|                    |                     | Pandoraea captiosa                  |     1 |
|                    |                     | Pandoraea cepalis                   |     2 |
|                    |                     | Pandoraea commovens                 |     2 |
|                    |                     | Pandoraea communis                  |     2 |
|                    |                     | Pandoraea eparura                   |     1 |
|                    |                     | Pandoraea faecigallinarum           |     1 |
|                    |                     | Pandoraea fibrosis                  |     2 |
|                    |                     | Pandoraea horticolens               |     1 |
|                    |                     | Pandoraea iniqua                    |     2 |
|                    |                     | Pandoraea morbifera                 |     1 |
|                    |                     | Pandoraea norimbergensis            |     1 |
|                    |                     | Pandoraea nosoerga                  |     3 |
|                    |                     | Pandoraea oxalativorans             |     1 |
|                    |                     | Pandoraea pneumonica                |     1 |
|                    |                     | Pandoraea pnomenusa                 |    10 |
|                    |                     | Pandoraea pulmonicola               |     2 |
|                    |                     | Pandoraea soli                      |     1 |
|                    |                     | Pandoraea sputorum                  |     2 |
|                    |                     | Pandoraea terrae                    |     1 |
|                    |                     | Pandoraea terrigena                 |     1 |
|                    |                     | Pandoraea thiooxydans               |     2 |
|                    |                     | Pandoraea vervacti                  |     1 |
|                    | ParaStreptomyces    | ParaStreptomyces acidicola          |     1 |
|                    |                     | ParaStreptomyces acidipaludis       |     1 |
|                    |                     | ParaStreptomyces acidiphila         |     1 |
|                    |                     | ParaStreptomyces acidisoli          |     1 |
|                    |                     | ParaStreptomyces agricolaris        |     1 |
|                    |                     | ParaStreptomyces antibiotica        |     1 |
|                    |                     | ParaStreptomyces aspalathi          |     7 |
|                    |                     | ParaStreptomyces atlantica          |     7 |
|                    |                     | ParaStreptomyces azotifigens        |     1 |
|                    |                     | ParaStreptomyces bonniea            |     1 |
|                    |                     | ParaStreptomyces bryophila          |     4 |
|                    |                     | ParaStreptomyces caballeronis       |     9 |
|                    |                     | ParaStreptomyces caffeinilytica     |     3 |
|                    |                     | ParaStreptomyces caffeinitolerans   |     1 |
|                    |                     | ParaStreptomyces caledonica         |     5 |
|                    |                     | ParaStreptomyces caribensis         |     9 |
|                    |                     | ParaStreptomyces diazotrophica      |     1 |
|                    |                     | ParaStreptomyces dilworthii         |     1 |
|                    |                     | ParaStreptomyces dinghuensis        |     1 |
|                    |                     | ParaStreptomyces dioscoreae         |     1 |
|                    |                     | ParaStreptomyces dipogonis          |     1 |
|                    |                     | ParaStreptomyces dokdonensis        |     1 |
|                    |                     | ParaStreptomyces domus              |     7 |
|                    |                     | ParaStreptomyces eburnea            |     2 |
|                    |                     | ParaStreptomyces edwinii            |     1 |
|                    |                     | ParaStreptomyces ferrariae          |     1 |
|                    |                     | ParaStreptomyces flava              |     1 |
|                    |                     | ParaStreptomyces franconis          |     1 |
|                    |                     | ParaStreptomyces fungorum           |    17 |
|                    |                     | ParaStreptomyces fynbosensis        |     1 |
|                    |                     | ParaStreptomyces gardini            |     2 |
|                    |                     | ParaStreptomyces ginsengisoli       |     2 |
|                    |                     | ParaStreptomyces ginsengiterrae     |     2 |
|                    |                     | ParaStreptomyces graminis           |     3 |
|                    |                     | ParaStreptomyces guartelaensis      |     1 |
|                    |                     | ParaStreptomyces haematera          |     1 |
|                    |                     | ParaStreptomyces hayleyella         |     1 |
|                    |                     | ParaStreptomyces heleia             |     1 |
|                    |                     | ParaStreptomyces hiiakae            |     1 |
|                    |                     | ParaStreptomyces hospita            |     7 |
|                    |                     | ParaStreptomyces humisilvae         |     1 |
|                    |                     | ParaStreptomyces kirstenboschensis  |     1 |
|                    |                     | ParaStreptomyces kururiensis        |     4 |
|                    |                     | ParaStreptomyces lacunae            |     1 |
|                    |                     | ParaStreptomyces lycopersici        |     1 |
|                    |                     | ParaStreptomyces madseniana         |     3 |
|                    |                     | ParaStreptomyces megapolitana       |     3 |
|                    |                     | ParaStreptomyces metrosideri        |     1 |
|                    |                     | ParaStreptomyces mimosarum          |     4 |
|                    |                     | ParaStreptomyces monticola          |     1 |
|                    |                     | ParaStreptomyces nemoris            |     6 |
|                    |                     | ParaStreptomyces pallida            |     1 |
|                    |                     | ParaStreptomyces panacisoli         |     1 |
|                    |                     | ParaStreptomyces phenazinium        |     3 |
|                    |                     | ParaStreptomyces phenoliruptrix     |     7 |
|                    |                     | ParaStreptomyces phosphatilytica    |     1 |
|                    |                     | ParaStreptomyces phymatum           |     2 |
|                    |                     | ParaStreptomyces piptadeniae        |     1 |
|                    |                     | ParaStreptomyces podalyriae         |     1 |
|                    |                     | ParaStreptomyces polaris            |     1 |
|                    |                     | ParaStreptomyces rhizosphaerae      |     1 |
|                    |                     | ParaStreptomyces rhynchosiae        |     2 |
|                    |                     | ParaStreptomyces ribeironis         |     1 |
|                    |                     | ParaStreptomyces sabiae             |     2 |
|                    |                     | ParaStreptomyces sacchari           |     3 |
|                    |                     | ParaStreptomyces saeva              |     3 |
|                    |                     | ParaStreptomyces sartisoli          |     1 |
|                    |                     | ParaStreptomyces silvatlantica      |     5 |
|                    |                     | ParaStreptomyces silviterrae        |     1 |
|                    |                     | ParaStreptomyces solisilvae         |     1 |
|                    |                     | ParaStreptomyces sprentiae          |     2 |
|                    |                     | ParaStreptomyces steynii            |     1 |
|                    |                     | ParaStreptomyces strydomiana        |     2 |
|                    |                     | ParaStreptomyces susongensis        |     1 |
|                    |                     | ParaStreptomyces tagetis            |     1 |
|                    |                     | ParaStreptomyces terrae             |     4 |
|                    |                     | ParaStreptomyces terricola          |     4 |
|                    |                     | ParaStreptomyces tropica            |    14 |
|                    |                     | ParaStreptomyces tuberum            |     3 |
|                    |                     | ParaStreptomyces ultramafica        |     1 |
|                    |                     | ParaStreptomyces unamae             |     3 |
|                    |                     | ParaStreptomyces xenovorans         |     2 |
|                    |                     | ParaStreptomyces youngii            |     5 |
|                    | Pararobbsia         | Pararobbsia alpina                  |     1 |
|                    |                     | Pararobbsia silviterrae             |     1 |
|                    | Polynucleobacter    | Polynucleobacter acidiphobus        |     1 |
|                    |                     | Polynucleobacter aenigmaticus       |     1 |
|                    |                     | Polynucleobacter alcilacus          |     1 |
|                    |                     | Polynucleobacter antarcticus        |     1 |
|                    |                     | Polynucleobacter arcticus           |     1 |
|                    |                     | Polynucleobacter asymbioticus       |     4 |
|                    |                     | Polynucleobacter bastaniensis       |     1 |
|                    |                     | Polynucleobacter brandtiae          |     1 |
|                    |                     | Polynucleobacter campilacus         |     1 |
|                    |                     | Polynucleobacter corsicus           |     1 |
|                    |                     | Polynucleobacter cosmopolitanus     |     1 |
|                    |                     | Polynucleobacter difficilis         |     1 |
|                    |                     | Polynucleobacter duraquae           |     1 |
|                    |                     | Polynucleobacter finlandensis       |     1 |
|                    |                     | Polynucleobacter hallstattensis     |     1 |
|                    |                     | Polynucleobacter hirudinilacicola   |     1 |
|                    |                     | Polynucleobacter ibericus           |     1 |
|                    |                     | Polynucleobacter kasalickyi         |     1 |
|                    |                     | Polynucleobacter nymphae            |     1 |
|                    |                     | Polynucleobacter paludilacus        |     1 |
|                    |                     | Polynucleobacter paneuropaeus       |   113 |
|                    |                     | Polynucleobacter parvulilacunae     |     1 |
|                    |                     | Polynucleobacter rarus              |     1 |
|                    |                     | Polynucleobacter sinensis           |     1 |
|                    |                     | Polynucleobacter sphagniphilus      |    11 |
|                    |                     | Polynucleobacter tropicus           |     1 |
|                    |                     | Polynucleobacter victoriensis       |     1 |
|                    |                     | Polynucleobacter wuianus            |     2 |
|                    |                     | Polynucleobacter yangtzensis        |     3 |
|                    | Quisquiliibacterium | Quisquiliibacterium transsilvanicum |     1 |
|                    | Ralstonia           | Ralstonia chuxiongensis             |     1 |
|                    |                     | Ralstonia insidiosa                 |    23 |
|                    |                     | Ralstonia mannitolilytica           |    14 |
|                    |                     | Ralstonia mojiangensis              |     5 |
|                    |                     | Ralstonia nicotianae                |     1 |
|                    |                     | Ralstonia pickettii                 |     3 |
|                    |                     | Ralstonia pseudosolanacearum        |    33 |
|                    |                     | Ralstonia solanacearum              |   121 |
|                    |                     | Ralstonia soli                      |     1 |
|                    |                     | Ralstonia syzygii                   |     7 |
|                    |                     | Ralstonia wenshanensis              |     2 |
|                    | Robbsia             | Robbsia andropogonis                |     4 |
|                    | Trinickia           | Trinickia caryophylli               |     5 |
|                    |                     | Trinickia dabaoshanensis            |     1 |
|                    |                     | Trinickia diaoshuihuensis           |     1 |
|                    |                     | Trinickia dinghuensis               |     1 |
|                    |                     | Trinickia fusca                     |     1 |
|                    |                     | Trinickia soli                      |     2 |
|                    |                     | Trinickia symbiotica                |     3 |
|                    | Zeimonas            | Zeimonas arvi                       |     1 |
|                    |                     | Zeimonas sediminis                  |     1 |
| Pseudonocardiaceae | Saccharothrix       | Saccharothrix violaceirubra         |     1 |

### For *protein families* #用于*蛋白质家族*

```shell
cd /mnt/c/shengxin/data/Streptomyces/

nwr template /mnt/c/shengxin/data/Streptomyces/assembly/Streptomyces.assembly.tsv \
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
| genus               | #species | #strains |
| ------------------- | -------: | -------: |
| Streptomyces        |       41 |     3962 |
| Caballeronia        |       23 |       36 |
| Chitinasiproducens  |        1 |        1 |
| Chitinimonas        |        3 |        4 |
| Cupriavidus         |       21 |      117 |
| Ephemeroptericola   |        1 |        1 |
| Formosimonas        |        1 |        1 |
| Lautropia           |        2 |        5 |
| Limnobacter         |        3 |        3 |
| Mycetohabitans      |        2 |        2 |
| Mycoavidus          |        1 |        1 |
| Pandoraea           |       28 |       63 |
| ParaStreptomyces    |       83 |      217 |
| Pararobbsia         |        2 |        2 |
| Polynucleobacter    |       29 |      157 |
| Quisquiliibacterium |        1 |        1 |
| Ralstonia           |       11 |      207 |
| Robbsia             |        1 |        4 |
| Saccharothrix       |        1 |        1 |
| Trinickia           |        7 |       14 |
| Zeimonas            |        2 |        2 |

## Collect Protein

```shell
cd /mnt/c/shengxin/data/Streptomyces/

nwr template /mnt/c/shengxin/data/Streptomyces/assembly/Streptomyces.assembly.tsv \
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
| #item                          | count      |
| ------------------------------ | ---------- |
| Proteins                       | 29,719,431 |
| Unique headers and annotations | 6,971,708  |
| Unique proteins                | 6,889,042  |
| all.replace.fa                 | 29,719,431 |
| all.annotation.tsv             | 29,719,432 |
| all.info.tsv                   | 29,719,432 |
## Phylogenetics with bac120 #系统发育

### Find corresponding Protein by `hmmsearch` #通过`hmmsearch`查找相应的蛋白质

```shell
cd /mnt/c/shengxin/data/Streptomyces/
mkdir -p HMM
E_VALUE=1e-1

# Find all genes
    cat Protein/species.tsv |
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
        > HMM/replace.tsv

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/DYW.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/DYWreplace.tsv

    cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Toxin.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Toxinreplace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Adeamin.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Adeamin-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/FdhD-NarQ.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/FdhD-NarQ-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/MafB19.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/MafB19-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Pput2613.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Pput2613-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/XOO.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/XOO-replace.tsv

     cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/YwqJ.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/YwqJ-replace.tsv

     cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/AICARFT_IMPCHas.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/AICARFT_IMPCHas-replace.tsv


        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Bd3614.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Bd3614-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/dCMP_cyt_deam_1.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/dCMP_cyt_deam-replace.tsv


        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Inv-AAD.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Inv-AAD-replace.tsv


        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/LmjF365940-deam.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/LmjF365940-deam-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/LpxI C.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/LpxIC-replace.tsv


        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/MafB19.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/MafB19-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/OTT_1508.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/OTT_1508-replace.tsv


        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Pput2613.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Pput2613-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/TM1506.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/TM1506-replace.tsv

        cat Protein/species.tsv |
        tsv-join -f ASSEMBLY/pass.lst -k 1 |
        tsv-join -e -f MinHash/abnormal.lst -k 1 |
        tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            if [[ ! -d ASSEMBLY/{2}/{1} ]]; then
                exit
            fi

            gzip -dcf ASSEMBLY/{2}/{1}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/Toxin.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/Toxin-replace.tsv
```
# 建立蛋白树
```shell
cd /mnt/c/shengxin/data/Streptomyces/
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
mkdir -p mnt/shengxin/data/Streptomyces/HMM/TIGRFAM
cd mnt/shengxin/data/Streptomyces/HMM/TIGRFAM
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

mkdir -p /mnt/c/shengxin/data/Streptomyces/HMM/bac120
cd /mnt/c/shengxin/data/Streptomyces/HMM/bac120


#下载bac120.tsv文件
mnt/shengxin/data/Streptomyces/HMM/bac120/

mkdir -p HMM

cat /mnt/c/shengxin/data/Streptomyces/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory HMM -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat /mnt/c/shengxin/data/Streptomyces/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L http://pfam.xfam.org/family/{}/hmm > HMM/{}.HMM
    '
#bac120进行系统发育
E_VALUE=1e-20

# Find all genes
cd /mnt/c/shengxin/data/Streptomyces/
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
cd /mnt/c/shengxin/data/Streptomyces/

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
    rsvg-convert -o tree/Streptomyces.marker.png

rsync -avP \
    /mnt/c/shengxin/data/Streptomyces/HMM/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Streptomyces/HMM
```

rsync -avP \
    /mnt/c/shengxin/data/Streptomyces/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Streptomyces

rsync -avP \
    /mnt/c/shengxin/data/Clostridium/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Clostridium

rsync -avP \
    /mnt/c/shengxin/data/Enterobacter/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Enterobacter

rsync -avP \
    /mnt/c/shengxin/data/Treponema/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Treponema

rsync -avP \
    /mnt/c/shengxin/data/Lachnospira/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Lachnospira

rsync -avP \
    /mnt/c/shengxin/data/Frankiaceae/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Frankiaceae

rsync -avP \
    /mnt/c/shengxin/data/Actinospica/ \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Actinospica