#!/bin/bash

echo "Downloading..."
curl "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv" > covid-19-data.csv
echo "...done."

cat covid-19-data.csv | grep Italy | cut -d\, -f5- | sed 's/,/\n/g' | (i=0; while read line; do echo $((i++))','$line; done) > covid-19-data-it.txt


