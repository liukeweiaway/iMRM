###iMRM

####Description:

iMRMF is able to simultaneously identify m6A, m5C, m1A, ψ and A-to-I modifications in Homo sapiens, Mus musculus and Saccharomyces cerevisiae.

####Installation:
- <span  style="color: #5bdaed; font-weight: bold">python3.6</span>
- pandas==0.23.3
- joblib==0.13.2
- xgboost==0.90

``` 
pip install pandas==0.23.3 joblib==0.13.2 xgboost==0.90
``` 
``` 
pip install -r yours/requirement.txt
``` 
####Optional arguments:
```
  -h, --help            show this help message and exit
  --addresses ADDRESSES
                        Liukeweiaway@hotmail.com
  -i INPUTFILE, --inputFile INPUTFILE
                        -i input.txt (The input file is a complete Fasta
                        format sequence.)
  -o OUTPUTFILE, --outputFile OUTPUTFILE
                        -o output.html (Results are saved under results
                        folder.)
  -s SPECIES, --species SPECIES
                        -m Human/Mouse/Yeast (Choose one from three species to
                        use.)
  -m MODIFICATION, --modification MODIFICATION
                        -m m6A/m5C/A-toI/pseudouridine/m1A/all (It should be
                        noted that A-to-I modification is limited by the
                        amount of data and can only be predicted in human.)
  -t THRESHOLD, --threshold THRESHOLD
                        -m low/normal/high (We offer 3 options based on the
                        difference in specificity, which are low, normal and
                        high.)
```
####Example:
```
python iMRM.py -i sequence.txt -o ccc.html -s Human -m all -t normal
```
***
Version number：V0.0.1 <br>
Updated date：2020-02-26 <br>
Web server: http://www.bioml.cn/XG_iRNA/home <br>
Email: Liukeweiaway@hotmai.com 
***