keepProp=$1
javac hydroplane/*.java
java hydroplane.Sieve -v  puc=$keepProp localDebug=false fn=/work-zfs/mschatz1/mkirsche/reads/ERR2173373.fastq ofn=uncontained_$keepProp.txt urf=uncontainednames_$keepProp.txt  nt=2

