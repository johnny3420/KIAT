## Comparing 50bp reads to 100bp reads

### Setup work space

```
cd
mkdir KIAT
cd KIAT/
mkdir -p Flowering_STAR/Alignments Flowering_STAR/Gene_Counts Flowering_STAR/Trimmed_Reads
```

### Creating adapter trimmed reads 50bp reads

```
cd KIAT/Flowering_STAR/Trimmed_Reads
./../flowering_trim.sh
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
