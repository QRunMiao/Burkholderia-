### List all ranks 列出所有等级

```shell
nwr member Burkholderia |
    grep -v " sp." |
    tsv-summarize -H -g rank --count |
    mlr --itsv --omd cat |
    perl -nl -e 's/-\s*\|$/-:|/; print'
```
| rank | count |
| --- | ---:|
| genus | 1 |
| species | 57 |
| no rank | 15 |
| species group | 2 |
| strain | 299 |
| subspecies | 1 |
```shell
nwr lineage Burkholderia |
    tsv-filter --str-ne 1:clade |
    tsv-filter --str-ne "1:no rank" |
    sed -n '/kingdom\tBacteria/,$p' |
    sed -E "s/\b(genus)\b/*\1*/"| # Highlight genus
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

 #查看伯克霍尔德属的分类上的基本情况   
```
| #rank | sci_name | tax_id |
| --- | --- | --- |
| superkingdom | Bacteria | 2 |
| phylum | Pseudomonadota | 1224 |
| class | Betaproteobacteria | 28216 |
| order | Burkholderiales | 80840 |
| family | Burkholderiaceae | 119060 |
| *genus* | Burkholderia | 32008 |

### Species with assemblies 具有组装的物种

Burkholderiaceae 是细菌界中的一个科

```shell
cd /mnt/c/shengxin
mkdir -p data/Burkholderia/summary
cd /mnt/c/shengxin/data/Burkholderia/summary

#找出和伯克霍尔德同一科(Burkholderiaceae)的所有属

# should have a valid name of genus
nwr member Burkholderiaceae -r genus |
    grep -v -i "Candidatus " |
    grep -v -i "candidate " |
    grep -v " sp." |
    grep -v " spp." |
    sed '1d' |
    sort -n -k1,1 \
    > genus.list.tsv

wc -l genus.list.tsv
# 23 genus.list.tsv

#伯克霍尔德同属的所有参考物种基因组信息
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

#伯克霍尔德同属的所有物种信息(genbank)
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
# 276 RS1.tsv
# 281 GB1.tsv
# 557 total

for C in RS GB; do
    for N in $(seq 1 1 10); do
        if [ -e "${C}${N}.tsv" ]; then
            printf "${C}${N}\t"
            cat ${C}${N}.tsv |
                tsv-summarize --sum 3
        fi
    done
done
#RS1     5695
#GB1     6474

```

## Download all assemblies

### Create .assembly.tsv

```shell
cd /mnt/c/shengxin/data/Burkholderia/summary

# Reference genome
echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus IN ('Saccharomyces')
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-select -H -f organism_name,species,genus,ftp_path,biosample,assembly_level,assembly_accession \
    > raw.tsv
  #酵母菌属的参考菌株的基因组信息

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
 # 伯克霍尔德同科的各属的参考菌株基因组信息

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
 #伯克霍尔德同科的各属的菌株基因组信息

cat raw.tsv |
    tsv-uniq |
    datamash check
#6476 lines, 7 fields

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
    > Burkholderia.assembly.tsv
 #创建简写名称的木霉属水平的各菌株基因组下载文件

 datamash check < Burkholderia.assembly.tsv
 #6476 lines, 5 fields

# find potential duplicate strains or assemblies 
cat Burkholderia.assembly.tsv |
    tsv-uniq -f 1 --repeated
#检查有没有重复

cat Burkholderia.assembly.tsv |
    tsv-filter --str-not-in-fld 2:ftp
    # 检查下载链接是否正确

# Save the file to another directory to prevent accidentally changing it
# cp Burkholderia.assembly.tsv /mnt/c/shengxin/data/Burkholderia/assembly/

# Cleaning
rm raw*.*sv
```
### Count before download 下载前计数

* `strains.taxon.tsv` - taxonomy info: species, genus, family, order, and class
Strains.taxon.tsv - 分类信息：物种、属、科、目和类

```shell
cd /mnt/c/shengxin/data/Burkholderia

nwr template assembly/Burkholderia.assembly.tsv \
    --count \
    --rank genus

# strains.taxon.tsv and taxa.tsv
bash Count/strains.sh #strains.taxon.tsv共6列，是为了统计各个菌株的数量以及所在的物种，属，科，目，纲的数量

#taxa.tsv内容
item	count
strain	6475
species	282
genus	23
family	2
order	2
class	2

#strains.taxon.tsv内容
B_aenigmatica_AU17325_GCF_002223275_1	Burkholderia aenigmatica	Burkholderia	Burkholderiaceae	Burkholderiales	Betaproteobacteria
B_aenigmatica_BCC0217_GCF_902500525_1	Burkholderia aenigmatica	Burkholderia	Burkholderiaceae	Burkholderiales	Betaproteobacteria

# genus.lst and genus.count.tsv
bash Count/rank.sh

#genus.lst 文件只有一列，提取了属的名称
#genus.count.tsv 统计伯克霍尔德属以及各属的species和strains去重后的数量

cat Count/genus.before.tsv |
    mlr --itsv --omd cat |
    perl -nl -e 'm/^\|\s*---/ and print qq(|---|--:|--:|) and next; print'
```
| genus | #species | #strains |
|---|--:|--:|
| Burkholderia | 41 | 4965 |
| Caballeronia | 28 | 50 |
| Chitinasiproducens | 1 | 1 |
| Chitinimonas | 3 | 4 |
| Cupriavidus | 22 | 155 |
| Ephemeroptericola | 1 | 1 |
| Formosimonas | 1 | 1 |
| Hydromonas | 1 | 1 |
| Lautropia | 2 | 19 |
| Limnobacter | 3 | 3 |
| Mycetohabitans | 2 | 2 |
| Mycoavidus | 1 | 1 |
| Pandoraea | 28 | 79 |
| Paraburkholderia | 91 | 266 |
| Pararobbsia | 2 | 2 |
| Paucimonas | 1 | 2 |
| Polynucleobacter | 31 | 210 |
| Quisquiliibacterium | 1 | 1 |
| Ralstonia | 11 | 689 |
| Robbsia | 1 | 5 |
| Saccharomyces | 1 | 1 |
| Trinickia | 7 | 15 |
| Zeimonas | 2 | 2 |

### Download and check 下载并检查
```shell
cd /mnt/c/shengxin/data/Burkholderia
nwr template ../Burkholderia/assembly/Burkholderia.assembly.tsv\
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
cp -r ~/data/Bacteria/ASSEMBLY/Burkholderia_* ~/qyl/data/Burkholderia
cp -r ~/data/Bacteria/ASSEMBLY/Cupriavidus_* ~/qyl/data/Burkholderia/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Pandoraea_* ~/qyl/data/Burkholderia/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Paraburkholderia_* ~/qyl/data/Burkholderia/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Polynucleobacter_* ~/qyl/data/Burkholderia/ASSEMBLY
cp -r ~/data/Bacteria/ASSEMBLY/Ralstonia_* ~/qyl/data/Burkholderia/ASSEMBLY

#下载到本地
rsync -avP \
wangq@202.119.37.251:qyl/data/Burkholderia/ASSEMBLY \
/mnt/c/shengxin/data/Burkholderia
#超算上运行此命令

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
    tsv-summarize -H --quantile "S:0.1,0.5" --quantile "N50:0.1,0.5"  --quantile "C:0.5,0.9"
    #计算数据的分位数
#S_pct10 S_pct50 N50_pct10       N50_pct50       C_pct50 C_pct90
#5555984 7123830 73229.2 201283  96      234

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
| url.tsv | 3 | 6,475 |
| check.lst | 1 | 6,475 |
| collect.tsv | 20 | 6,476 |
| n50.tsv | 4 | 6,476 |
| n50.pass.tsv | 4 | 5,308 |
| collect.pass.tsv | 23 | 5,308 |
| pass.lst | 1 | 5,307 |
| omit.lst | 1 | 225 |
| rep.lst | 1 | 268 |

# 统计下载的文件行数：
# collect.pass.tsv ：下载的基因组通过了n50的检验 -5308行（有一行标题）；pass.lst 5307行（不含标题），取collect.pass.tsv的第一列名称。
# omit.lst ：没有蛋白序列信息，文件只有一列，基因组名 225个基因组
# rep.lst : 含有参考序列的基因组 268个基因组
```
### Rsync to hpcc 把本地文件与超算同步

```bash
rsync -avP \
    /mnt/c/shengxin/data/Burkholderia/assembly/ \
    wangq@202.119.37.251:qyl/data/Burkholderia/ASSEMBLY
#本地运行

rsync -avP \
    /mnt/c/shengxin/data/Burkholderia/summary/ \
    wangq@202.119.37.251:qyl/data/Burkholderia/summary
```

## BioSample 生物样本
```shell
cd ~/qyl/data/Burkholderia

ulimit -n `ulimit -Hn`

cp ASSEMBLY/Burkholderia.assembly.tsv summary/
nwr template /mnt/c/shengxin/data/Burkholderia/assembly/Burkholderia.assembly.tsv \
    --bs

head BioSample/sample.tsv
#SAMD00000356    Parab_fer_NBRC_106233_GCF_000685035_1   Paraburkholderia_ferrariae
#SAMD00000357    Parab_mim_NBRC_106338_GCF_000739815_1   Paraburkholderia_mimosarum

bash BioSample/download.sh
#在本地下载

# Ignore rare attributes#忽略稀少属性
bash BioSample/collect.sh 10

datamash check < BioSample/biosample.tsv
#6472 lines, 107 fields

cp Biosample/attributes.lst summary/ 
cp Biosample/biosample.tsv summary/

#biosample.tsv 文件内容
#name   BioSample       sample name     observed biotic relationship    collection date environmental medium    geographic location     isolation and growth condition  latitude and longitude  locus_tag_prefix        number of replicons    project name     reference for biomaterial       strain  trophic level   source material identifiers     isolation sourchost     External Id     Submitter Id    culture collection      host health state       anonymized name common name    serovar  subject id      supplier_name   investigation type      sequencing method       assembly quality        assembly software       binning parameters      binning software        completeness score      completeness software   contamination score     metagenomic source      sample derived from     taxonomic identity marker       scientific_name local environmental context     broker name     isolate Alias   SRA accession   Title   alias   strain_synonym  substrain      anonymized_name  relationship to oxygen  collected by    passage history GOLD Stamp ID   Gram Staining   disease environment     Cell Shape      Motility        Temperature Range       Sporulation     Phenotypes      Isolation Site  Temperature Optimum     host disease    sample type     environmental package   biomaterial provider    depth   alternate strain name   organism modifier note  host age        note    host sex        host disease outcome    genotype        host tissue sampled     alternate_ID    elevation       identified by   host description        serotype        pathotype      subgroup subtype specimen_category       provider        collection method       uFORGE_Sample_ID        pH      altitudhost disease stage       host subject id isolate name alias      identification method   biovar  phylotype       sample size     metagenomic     derived from    host taxonomy ID        source type     Phylotype       temperature     sample storage location sample storage temperature      host infra specific name
Parab_fer_NBRC_106233_GCF_000685035_1   SAMD00000356    NBRC 106233     free living             Iron Mine       Brazil:Minas Gerais     17012573                BFE01S          Burkholderia ferrariae NBRC 106233 genome sequencing project   NBRC 106233      heterotroph     NBRC 106233

rsync -avP \
    /mnt/c/shengxin/data/Burkholderia/Biosample/ \
    wangq@202.119.37.251:qyl/data/Burkholderia/Biosample
```

## MinHash
```shell
cd /mnt/c/shengxin/data/Burkholderia

nwr template /mnt/c/shengxin/data/Burkholderia/summary/Burkholderia.assembly.tsv \
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
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/Biosample/ \
    /mnt/c/shengxin/data/Burkholderia/Biosample

rsync -avP \
    wangq@202.119.37.251:/share/home/wangq/qyl/data/Burkholderia/MinHash/ \
    /mnt/c/shengxin/data/Burkholderia/MinHash
    #本地运行   
```

### Condense branches in the minhash tree
```shell
mkdir -p cd /mnt/c/shengxin/data/Burkholderia/tree
cd /mnt/c/shengxin/data/Burkholderia/tree

nw_reroot ../MinHash/tree.nwk S_cer_S288C |
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
    rsvg-convert -o Burkholderia.minhash.png
#报错：超过了 32767 像素，这是 rsvg-convert 当前无法处理的最大尺寸限制。
nw_display -s -b 'visibility:hidden' -w 800 -v 15 minhash.species.newick | rsvg-convert -o Burkholderia.minhash.png
#还是过大

nw_display -s -b 'visibility:hidden' -w 600 -v 10 minhash.species.newick |
    rsvg-convert -o Burkholderia.minhash.png

```
## Count valid species and strains #计算有效物种数和菌株数

### For *genomic alignments* 用于*基因组比对*

```shell
cd /mnt/c/shengxin/data/Burkholderia
nwr template /mnt/c/shengxin/data/Burkholderia/assembly/Burkholderia.assembly.tsv \
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
| genus | #species | #strains |
|---|--:|--:|
| Burkholderia | 26 | 2273 |
| Caballeronia | 21 | 30 |
| Chitinasiproducens | 1 | 1 |
| Chitinimonas | 3 | 4 |
| Cupriavidus | 17 | 70 |
| Ephemeroptericola | 1 | 1 |
| Formosimonas | 1 | 1 |
| Lautropia | 2 | 5 |
| Limnobacter | 3 | 3 |
| Mycetohabitans | 2 | 2 |
| Mycoavidus | 1 | 1 |
| Pandoraea | 24 | 56 |
| Paraburkholderia | 70 | 164 |
| Pararobbsia | 2 | 2 |
| Polynucleobacter | 28 | 153 |
| Quisquiliibacterium | 1 | 1 |
| Ralstonia | 9 | 87 |
| Robbsia | 1 | 4 |
| Saccharomyces | 1 | 1 |
| Trinickia | 5 | 9 |
| Zeimonas | 2 | 2 |


| #family | genus | species | count |
| --- | --- | --- | ---:|
| Burkholderiaceae | Burkholderia | Burkholderia aenigmatica | 20 |
|  |  | Burkholderia ambifaria | 89 |
|  |  | Burkholderia anthina | 27 |
|  |  | Burkholderia arboris | 6 |
|  |  | Burkholderia catarinensis | 1 |
|  |  | Burkholderia cenocepacia | 395 |
|  |  | Burkholderia cepacia | 234 |
|  |  | Burkholderia contaminans | 86 |
|  |  | Burkholderia diffusa | 16 |
|  |  | Burkholderia dolosa | 20 |
|  |  | Burkholderia gladioli | 260 |
|  |  | Burkholderia glumae | 63 |
|  |  | Burkholderia guangdongensis | 1 |
|  |  | Burkholderia humptydooensis | 3 |
|  |  | Burkholderia lata | 16 |
|  |  | Burkholderia latens | 6 |
|  |  | Burkholderia mallei | 58 |
|  |  | Burkholderia mayonis | 2 |
|  |  | Burkholderia metallica | 7 |
|  |  | Burkholderia multivorans | 490 |
|  |  | Burkholderia oklahomensis | 9 |
|  |  | Burkholderia orbicola | 3 |
|  |  | Burkholderia paludis | 3 |
|  |  | Burkholderia perseverans | 1 |
|  |  | Burkholderia plantarii | 5 |
|  |  | Burkholderia pseudomallei | 1504 |
|  |  | Burkholderia pseudomultivorans | 10 |
|  |  | Burkholderia puraquae | 2 |
|  |  | Burkholderia pyrrocinia | 3 |
|  |  | Burkholderia reimsis | 1 |
|  |  | Burkholderia savannae | 4 |
|  |  | Burkholderia semiarida | 4 |
|  |  | Burkholderia seminalis | 15 |
|  |  | Burkholderia singularis | 1 |
|  |  | Burkholderia sola | 1 |
|  |  | Burkholderia stabilis | 2 |
|  |  | Burkholderia stagnalis | 101 |
|  |  | Burkholderia territorii | 37 |
|  |  | Burkholderia thailandensis | 27 |
|  |  | Burkholderia ubonensis | 297 |
|  |  | Burkholderia vietnamiensis | 139 |
|  | Caballeronia | Caballeronia arationis | 1 |
|  |  | Caballeronia arvi | 1 |
|  |  | Caballeronia calidae | 1 |
|  |  | Caballeronia catudaia | 1 |
|  |  | Caballeronia concitans | 1 |
|  |  | Caballeronia cordobensis | 2 |
|  |  | Caballeronia fortuita | 1 |
|  |  | Caballeronia glathei | 2 |
|  |  | Caballeronia glebae | 1 |
|  |  | Caballeronia grimmiae | 3 |
|  |  | Caballeronia hypogeia | 1 |
|  |  | Caballeronia insecticola | 1 |
|  |  | Caballeronia jiangsuensis | 1 |
|  |  | Caballeronia novacaledonica | 3 |
|  |  | Caballeronia pedi | 1 |
|  |  | Caballeronia peredens | 1 |
|  |  | Caballeronia ptereochthonis | 1 |
|  |  | Caballeronia sordidicola | 3 |
|  |  | Caballeronia telluris | 1 |
|  |  | Caballeronia temeraria | 1 |
|  |  | Caballeronia terrestris | 1 |
|  |  | Caballeronia turbans | 1 |
|  |  | Caballeronia zhejiangensis | 6 |
|  | Chitinasiproducens | Chitinasiproducens palmae | 1 |
|  | Chitinimonas | Chitinimonas arctica | 1 |
|  |  | Chitinimonas koreensis | 2 |
|  |  | Chitinimonas taiwanensis | 1 |
|  | Cupriavidus | Cupriavidus agavae | 1 |
|  |  | Cupriavidus alkaliphilus | 7 |
|  |  | Cupriavidus basilensis | 5 |
|  |  | Cupriavidus campinensis | 3 |
|  |  | Cupriavidus cauae | 2 |
|  |  | Cupriavidus gilardii | 17 |
|  |  | Cupriavidus laharis | 1 |
|  |  | Cupriavidus malaysiensis | 1 |
|  |  | Cupriavidus metallidurans | 15 |
|  |  | Cupriavidus nantongensis | 2 |
|  |  | Cupriavidus necator | 5 |
|  |  | Cupriavidus neocaledonicus | 3 |
|  |  | Cupriavidus numazuensis | 1 |
|  |  | Cupriavidus oxalaticus | 3 |
|  |  | Cupriavidus pampae | 1 |
|  |  | Cupriavidus pauculus | 2 |
|  |  | Cupriavidus pinatubonensis | 2 |
|  |  | Cupriavidus plantarum | 6 |
|  |  | Cupriavidus respiraculi | 2 |
|  |  | Cupriavidus taiwanensis | 37 |
|  |  | Cupriavidus yeoncheonensis | 1 |
|  | Ephemeroptericola | Ephemeroptericola cinctiostellae | 1 |
|  | Formosimonas | Formosimonas limnophila | 1 |
|  | Lautropia | Lautropia dentalis | 1 |
|  |  | Lautropia mirabilis | 4 |
|  | Limnobacter | Limnobacter alexandrii | 1 |
|  |  | Limnobacter humi | 1 |
|  |  | Limnobacter thiooxidans | 1 |
|  | Mycetohabitans | Mycetohabitans endofungorum | 1 |
|  |  | Mycetohabitans rhizoxinica | 1 |
|  | Mycoavidus | Mycoavidus cysteinexigens | 1 |
|  | Pandoraea | Pandoraea anapnoica | 1 |
|  |  | Pandoraea anhela | 1 |
|  |  | Pandoraea apista | 18 |
|  |  | Pandoraea aquatica | 1 |
|  |  | Pandoraea bronchicola | 1 |
|  |  | Pandoraea capi | 1 |
|  |  | Pandoraea captiosa | 1 |
|  |  | Pandoraea cepalis | 2 |
|  |  | Pandoraea commovens | 2 |
|  |  | Pandoraea communis | 2 |
|  |  | Pandoraea eparura | 1 |
|  |  | Pandoraea faecigallinarum | 1 |
|  |  | Pandoraea fibrosis | 2 |
|  |  | Pandoraea horticolens | 1 |
|  |  | Pandoraea iniqua | 2 |
|  |  | Pandoraea morbifera | 1 |
|  |  | Pandoraea norimbergensis | 1 |
|  |  | Pandoraea nosoerga | 3 |
|  |  | Pandoraea oxalativorans | 1 |
|  |  | Pandoraea pneumonica | 1 |
|  |  | Pandoraea pnomenusa | 10 |
|  |  | Pandoraea pulmonicola | 2 |
|  |  | Pandoraea soli | 1 |
|  |  | Pandoraea sputorum | 2 |
|  |  | Pandoraea terrae | 1 |
|  |  | Pandoraea terrigena | 1 |
|  |  | Pandoraea thiooxydans | 2 |
|  |  | Pandoraea vervacti | 1 |
|  | Paraburkholderia | Paraburkholderia acidicola | 1 |
|  |  | Paraburkholderia acidipaludis | 1 |
|  |  | Paraburkholderia acidiphila | 1 |
|  |  | Paraburkholderia acidisoli | 1 |
|  |  | Paraburkholderia agricolaris | 1 |
|  |  | Paraburkholderia antibiotica | 1 |
|  |  | Paraburkholderia aspalathi | 7 |
|  |  | Paraburkholderia atlantica | 7 |
|  |  | Paraburkholderia azotifigens | 1 |
|  |  | Paraburkholderia bonniea | 1 |
|  |  | Paraburkholderia bryophila | 4 |
|  |  | Paraburkholderia caballeronis | 9 |
|  |  | Paraburkholderia caffeinilytica | 3 |
|  |  | Paraburkholderia caffeinitolerans | 1 |
|  |  | Paraburkholderia caledonica | 5 |
|  |  | Paraburkholderia caribensis | 9 |
|  |  | Paraburkholderia diazotrophica | 1 |
|  |  | Paraburkholderia dilworthii | 1 |
|  |  | Paraburkholderia dinghuensis | 1 |
|  |  | Paraburkholderia dioscoreae | 1 |
|  |  | Paraburkholderia dipogonis | 1 |
|  |  | Paraburkholderia dokdonensis | 1 |
|  |  | Paraburkholderia domus | 7 |
|  |  | Paraburkholderia eburnea | 2 |
|  |  | Paraburkholderia edwinii | 1 |
|  |  | Paraburkholderia ferrariae | 1 |
|  |  | Paraburkholderia flava | 1 |
|  |  | Paraburkholderia franconis | 1 |
|  |  | Paraburkholderia fungorum | 17 |
|  |  | Paraburkholderia fynbosensis | 1 |
|  |  | Paraburkholderia gardini | 2 |
|  |  | Paraburkholderia ginsengisoli | 2 |
|  |  | Paraburkholderia ginsengiterrae | 2 |
|  |  | Paraburkholderia graminis | 3 |
|  |  | Paraburkholderia guartelaensis | 1 |
|  |  | Paraburkholderia haematera | 1 |
|  |  | Paraburkholderia hayleyella | 1 |
|  |  | Paraburkholderia heleia | 1 |
|  |  | Paraburkholderia hiiakae | 1 |
|  |  | Paraburkholderia hospita | 7 |
|  |  | Paraburkholderia humisilvae | 1 |
|  |  | Paraburkholderia kirstenboschensis | 1 |
|  |  | Paraburkholderia kururiensis | 4 |
|  |  | Paraburkholderia lacunae | 1 |
|  |  | Paraburkholderia lycopersici | 1 |
|  |  | Paraburkholderia madseniana | 3 |
|  |  | Paraburkholderia megapolitana | 3 |
|  |  | Paraburkholderia metrosideri | 1 |
|  |  | Paraburkholderia mimosarum | 4 |
|  |  | Paraburkholderia monticola | 1 |
|  |  | Paraburkholderia nemoris | 6 |
|  |  | Paraburkholderia pallida | 1 |
|  |  | Paraburkholderia panacisoli | 1 |
|  |  | Paraburkholderia phenazinium | 3 |
|  |  | Paraburkholderia phenoliruptrix | 7 |
|  |  | Paraburkholderia phosphatilytica | 1 |
|  |  | Paraburkholderia phymatum | 2 |
|  |  | Paraburkholderia piptadeniae | 1 |
|  |  | Paraburkholderia podalyriae | 1 |
|  |  | Paraburkholderia polaris | 1 |
|  |  | Paraburkholderia rhizosphaerae | 1 |
|  |  | Paraburkholderia rhynchosiae | 2 |
|  |  | Paraburkholderia ribeironis | 1 |
|  |  | Paraburkholderia sabiae | 2 |
|  |  | Paraburkholderia sacchari | 3 |
|  |  | Paraburkholderia saeva | 3 |
|  |  | Paraburkholderia sartisoli | 1 |
|  |  | Paraburkholderia silvatlantica | 5 |
|  |  | Paraburkholderia silviterrae | 1 |
|  |  | Paraburkholderia solisilvae | 1 |
|  |  | Paraburkholderia sprentiae | 2 |
|  |  | Paraburkholderia steynii | 1 |
|  |  | Paraburkholderia strydomiana | 2 |
|  |  | Paraburkholderia susongensis | 1 |
|  |  | Paraburkholderia tagetis | 1 |
|  |  | Paraburkholderia terrae | 4 |
|  |  | Paraburkholderia terricola | 4 |
|  |  | Paraburkholderia tropica | 14 |
|  |  | Paraburkholderia tuberum | 3 |
|  |  | Paraburkholderia ultramafica | 1 |
|  |  | Paraburkholderia unamae | 3 |
|  |  | Paraburkholderia xenovorans | 2 |
|  |  | Paraburkholderia youngii | 5 |
|  | Pararobbsia | Pararobbsia alpina | 1 |
|  |  | Pararobbsia silviterrae | 1 |
|  | Polynucleobacter | Polynucleobacter acidiphobus | 1 |
|  |  | Polynucleobacter aenigmaticus | 1 |
|  |  | Polynucleobacter alcilacus | 1 |
|  |  | Polynucleobacter antarcticus | 1 |
|  |  | Polynucleobacter arcticus | 1 |
|  |  | Polynucleobacter asymbioticus | 4 |
|  |  | Polynucleobacter bastaniensis | 1 |
|  |  | Polynucleobacter brandtiae | 1 |
|  |  | Polynucleobacter campilacus | 1 |
|  |  | Polynucleobacter corsicus | 1 |
|  |  | Polynucleobacter cosmopolitanus | 1 |
|  |  | Polynucleobacter difficilis | 1 |
|  |  | Polynucleobacter duraquae | 1 |
|  |  | Polynucleobacter finlandensis | 1 |
|  |  | Polynucleobacter hallstattensis | 1 |
|  |  | Polynucleobacter hirudinilacicola | 1 |
|  |  | Polynucleobacter ibericus | 1 |
|  |  | Polynucleobacter kasalickyi | 1 |
|  |  | Polynucleobacter nymphae | 1 |
|  |  | Polynucleobacter paludilacus | 1 |
|  |  | Polynucleobacter paneuropaeus | 113 |
|  |  | Polynucleobacter parvulilacunae | 1 |
|  |  | Polynucleobacter rarus | 1 |
|  |  | Polynucleobacter sinensis | 1 |
|  |  | Polynucleobacter sphagniphilus | 11 |
|  |  | Polynucleobacter tropicus | 1 |
|  |  | Polynucleobacter victoriensis | 1 |
|  |  | Polynucleobacter wuianus | 2 |
|  |  | Polynucleobacter yangtzensis | 3 |
|  | Quisquiliibacterium | Quisquiliibacterium transsilvanicum | 1 |
|  | Ralstonia | Ralstonia chuxiongensis | 1 |
|  |  | Ralstonia insidiosa | 23 |
|  |  | Ralstonia mannitolilytica | 14 |
|  |  | Ralstonia mojiangensis | 5 |
|  |  | Ralstonia nicotianae | 1 |
|  |  | Ralstonia pickettii | 3 |
|  |  | Ralstonia pseudosolanacearum | 33 |
|  |  | Ralstonia solanacearum | 121 |
|  |  | Ralstonia soli | 1 |
|  |  | Ralstonia syzygii | 7 |
|  |  | Ralstonia wenshanensis | 2 |
|  | Robbsia | Robbsia andropogonis | 4 |
|  | Trinickia | Trinickia caryophylli | 5 |
|  |  | Trinickia dabaoshanensis | 1 |
|  |  | Trinickia diaoshuihuensis | 1 |
|  |  | Trinickia dinghuensis | 1 |
|  |  | Trinickia fusca | 1 |
|  |  | Trinickia soli | 2 |
|  |  | Trinickia symbiotica | 3 |
|  | Zeimonas | Zeimonas arvi | 1 |
|  |  | Zeimonas sediminis | 1 |
| Saccharomycetaceae | Saccharomyces | Saccharomyces cerevisiae | 1 |

### For *protein families* #用于*蛋白质家族*

```shell
cd /mnt/c/shengxin/data/Burkholderia/

nwr template /mnt/c/shengxin/data/Burkholderia/assembly/Burkholderia.assembly.tsv \
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
| genus | #species | #strains |
|---|--:|--:|
| Burkholderia | 41 | 3962 |
| Caballeronia | 23 | 36 |
| Chitinasiproducens | 1 | 1 |
| Chitinimonas | 3 | 4 |
| Cupriavidus | 21 | 117 |
| Ephemeroptericola | 1 | 1 |
| Formosimonas | 1 | 1 |
| Lautropia | 2 | 5 |
| Limnobacter | 3 | 3 |
| Mycetohabitans | 2 | 2 |
| Mycoavidus | 1 | 1 |
| Pandoraea | 28 | 63 |
| Paraburkholderia | 83 | 217 |
| Pararobbsia | 2 | 2 |
| Polynucleobacter | 29 | 157 |
| Quisquiliibacterium | 1 | 1 |
| Ralstonia | 11 | 207 |
| Robbsia | 1 | 4 |
| Saccharomyces | 1 | 1 |
| Trinickia | 7 | 14 |
| Zeimonas | 2 | 2 |
## Collect proteins

```shell
cd /mnt/c/shengxin/data/Burkholderia/

nwr template /mnt/c/shengxin/data/Burkholderia/assembly/Burkholderia.assembly.tsv \
    --pro \
    --in ASSEMBLY/pass.lst \
    --not-in MinHash/abnormal.lst \
    --not-in ASSEMBLY/omit.lst

# * --pro: Protein/
#     * One TSV file
#         * species.tsv
#     * collect.sh

# collect proteins
bash Protein/collect.sh

cat Protein/counts.tsv |
    mlr --itsv --omd cat
```
| #item | count |
| --- | --- |
| Proteins | 29,719,225 |
| Unique headers and annotations | 6,971,502 |
| Unique proteins | 6,888,836 |
| all.replace.fa | 29,719,225 |
| all.annotation.tsv | 29,719,226 |
| all.info.tsv | 29,719,226 |
## Phylogenetics with bac120 #系统发育

### Find corresponding proteins by `hmmsearch` #通过`hmmsearch`查找相应的蛋白质

```shell
cd /mnt/c/shengxin/data/Burkholderia/
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
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw HMM/DddA-like.HMM - |
                grep '>>' |
                perl -nl -e ' m(>>\s+(\S+)) and printf qq(%s\t%s\n), \$1, {1}; '
        " \
        > Protein/replace.tsv

```
# 建立蛋白树
```shell
awk -F'\t' '{print $2 "_" $1}' Protein/replace.tsv > Protein/3.tsv

cd /mnt/c/shengxin/data/Burkholderia/
faops some Protein/all.replace.fa.gz <(tsv-select -f 1 Protein/3.tsv) Protein/DddA-like.fa

muscle -in Protein/DddA-like.fa -out Protein/DddA-like.aln.fa

FastTree Protein/DddA-like.aln.fa > Protein/DddA-like.aln.newick

nw_reroot Protein/DddA-like.aln.newick $(nw_labels Protein/DddA-like.aln.newick | grep -E "S_cer_S288C") |
    nw_order -c n - \
    > Protein/DddA-like.reoot.newick

nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/DddA-like.reoot.newick | rsvg-convert -o tree/DddA-like.reoot.png
```
## Phylogenetics with bac120 #系统发育

### Find corresponding proteins by `hmmsearch` #通过`hmmsearch`查找相应的蛋白质
```shell
#TIGRFAM
mkdir -p mnt/shengxin/data/Burkholderia/HMM/TIGRFAM
cd mnt/shengxin/data/Burkholderia/HMM/TIGRFAM
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

mkdir -p /mnt/c/shengxin/data/Burkholderia/HMM/bac120
cd /mnt/c/shengxin/data/Burkholderia/HMM/bac120


#下载bac120.tsv文件
mnt/shengxin/data/Burkholderia/HMM/bac120/

mkdir -p HMM

cat /mnt/c/shengxin/data/Burkholderia/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        tar --directory HMM -xzvf ../TIGRFAM/TIGRFAMs_14.0_HMM.tar.gz {}.HMM
    '

cat /mnt/c/shengxin/data/Burkholderia/HMM/bac120/bac120.tsv |
    sed '1d' |
    tsv-select -f 1 |
    grep -v '^TIGR' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        curl -L http://pfam.xfam.org/family/{}/hmm > HMM/{}.HMM
    '
#bac120进行系统发育
E_VALUE=1e-20

# Find all genes
cd /mnt/c/shengxin/data/Burkholderia/
cat Protein/species.tsv |cut -f 1 >Protein/1.tsv
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
cd /mnt/c/shengxin/data/Burkholderia/

cat HMM/bac120/bac120.tsv | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    tsv-summarize --quantile 1:0.25,0.5,0.75

cat HMM/bac120/bac120.tsv | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {}
        cat Protein/{}/replace.tsv |
            wc -l
    ' |
    paste - - |
    tsv-filter --invert --ge 2:1800 --le 2:2600 |
    cut -f 1 \
    > Protein/bac120.omit.lst

# Extract sequences
# Multiple copies slow down the alignment process
cat HMM/marker.lst |
    grep -v -Fx -f Protein/marker.omit.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> marker [{}]"

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
cat HMM/marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 "==> marker [{}]"
        if [ ! -s Protein/{}/{}.pro.fa ]; then
            exit
        fi
        if [ -s Protein/{}/{}.aln.fa ]; then
            exit
        fi

        muscle -quiet -in Protein/{}/{}.pro.fa -out Protein/{}/{}.aln.fa
    '

for marker in $(cat HMM/marker.lst); do
    echo >&2 "==> marker [${marker}]"
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
for marker in $(cat HMM/marker.lst); do
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
    > Protein/fungi61.aln.fas

cat Protein/species.tsv |
    tsv-join -f ASSEMBLY/pass.lst -k 1 |
    tsv-join -e -f MinHash/abnormal.lst -k 1 |
    tsv-join -e -f ASSEMBLY/omit.lst -k 1 |
    cut -f 1 |
    fasops concat Protein/fungi61.aln.fas stdin -o Protein/fungi61.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in Protein/fungi61.aln.fa -out Protein/fungi61.trim.fa -automated1

faops size Protein/fungi61.*.fa |
    tsv-uniq -f 2 |
    cut -f 2
#28706
#20432

# To make it faster
FastTree -fastest -noml Protein/fungi61.trim.fa > Protein/fungi61.trim.newick

nw_reroot Protein/fungi61.trim.newick S_cer_S288C |
    nw_order -c n - \
    > Protein/fungi61.reroot.newick

# png
nw_display -s -b 'visibility:hidden' -w 1200 -v 20 Protein/fungi61.reroot.newick |
    rsvg-convert -o tree/Burkholderia.marker.png
```
