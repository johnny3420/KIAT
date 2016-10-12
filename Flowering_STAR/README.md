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
### Running STAR on trimmed reads

```
./short_star_run.sh
```
