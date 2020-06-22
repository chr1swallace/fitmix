fit a mixture to eqtl summary stats

Code adapted from original by Paul Kirk

## 1d clustering with stan

```{sh}
for tt in skin ebv blood ; do
    qR.rb -c 4 -j stn-${tt} -t "12:00:00" -r ./stan_1d.R --args t1=${tt}
done
```

## Run 1d clustering

```{sh}
rm ~/scratch/asemix1dnosort*.RData
k=2
id=$RANDOM
for i in 1 2 3; do
    for tt in skin ebv blood ; do
	for nits in 20000; do
        id=$(( $id + 1))
	    qR.rb  -j ${tt}-${id}-${nits} -t "12:00:00" -r ./fitmix_1d_run.R --args k=$k id=$id nits=$nits thin=10 t1=${tt}
	done
    done
done
```

## Run 2d clustering

```{sh}
rm ~/scratch/asemix2dv6*.RData

id=$RANDOM
for i in 1 2 3; do
    for nits in 10000 200000 500000; do
        id=$(( $id + 1))
	qR.rb  -j blood-ebv-${id}-${nits} -t "12:00:00" -r ./fitmix_2d_run_v6.R --args k=4 id=$id nits=$nits thin=10 t1="blood" t2="ebv"
        id=$(( $id + 1))
	qR.rb  -j blood-skin-${id}-${nits} -t "12:00:00" -r ./fitmix_2d_run_v6.R --args k=4 id=$id nits=$nits thin=10 t1="blood" t2="skin"
        id=$(( $id + 1))
	qR.rb  -j skin-ebv-${id}-${nits} -t "12:00:00" -r ./fitmix_2d_run_v6.R --args k=4 id=$id nits=$nits thin=10 t1="skin" t2="ebv"
    done
done

```

