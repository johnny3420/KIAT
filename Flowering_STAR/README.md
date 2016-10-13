## Comparing 50bp reads to 100bp reads

### Setup work space

```
cd
mkdir KIAT
cd KIAT/
mkdir -p Flowering_STAR/Alignments Flowering_STAR/Gene_Counts Flowering_STAR/50bp_Trimmed_Reads Flowering_STAR/100bp_Trimmed_Reads
```

### Creating adapter trimmed 50bp reads

```
./../50bp_flowering_trim.sh
```

### Creating adapter trimmed 100bp reads

```
./../100bp_flowering_trim.sh
```

### Running STAR on 50bp trimmed reads

```
./short_star_run.sh
```

### Running STAR on 100bp single end trimmed reads

```
./long_star_run.sh
```

### Collecting Gene Counts

```
./gene_counts.sh
```
### R work see R_Scripts folder

![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/Combined.DE.barplot.png "Combined DE barplot")
![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/DE.barplot.png "DE barplot")
![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/Mapping.barplot.4.png "Mapping Rates")

#### Up regulated genes

![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/Up_Regulated_Venn.png "Up Regulated Genes")

#### Down regulated genes

![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/Down_Regulated_Venn.png "Down Regulated Genes")

#### Same regulated genes

![alt text](https://github.com/johnny3420/KIAT/blob/master/Graphs/Same_Regulated_Venn.png "Same Regulated Genes")
