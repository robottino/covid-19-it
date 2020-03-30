#!/bin/bash

mv "dpc-covid19-ita-andamento-nazionale.csv" "dpc-covid19-ita-andamento-nazionale-$(date +"%m-%d-%Y").csv"
curl "https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv" | cut -d\, -f3-12 | tail -n +2 > dpc-covid19-ita-andamento-nazionale.csv
