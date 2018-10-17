for i in *.txt; do cat $i | sed '1s/^/\"RowID\"\t/' > temp ; mv $i $i.old ; mv temp $i ; done`
`rm *.old