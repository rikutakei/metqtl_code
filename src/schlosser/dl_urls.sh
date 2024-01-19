#! /bin/bash

mkdir -p data/schlosser_metqtl/{plasma,urine}

parallel -j25 wget -P data/schlosser_metqtl/plasma/ {} ::: $(cat data/schlosser_metqtl/plasma_urls.txt)
parallel -j25 wget -P data/schlosser_metqtl/urine/ {} ::: $(cat data/schlosser_metqtl/urine_urls.txt)

