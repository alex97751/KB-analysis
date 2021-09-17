# KB-analysis
## Code For the Bachelor Thesis: Disambiguation and Normalization of Genomic Variants

Code by Alexander Harrisson. This was written for the Bachelor thesis


## Installation

Download the dbSNP DB as described in [SETH](https://rockt.github.io/SETH/). 
Download the Ensembl Database as described [here](https://m.ensembl.org/info/docs/webcode/mirror/install/ensembl-data.html).
Additionally the Ensembl DB for homo_sapiens_core_75_37 is being used as well.
After downloading and building these DBs, install the following Python packages:

```sh
pip install hgvs
pip install pandas
pip install requests
pip install mutalyzer_hgvs_parser
pip install mysql-connector-python
pip install psycopg2
pip install matplotlib
pip install matplotlib_venn
```
## Running the Code
Now you can download the KBs from their respective websites as seen in the thesis. Change the path to their respective locations on your local machine.
Put the function of the KB you want to analyse in the main function with either _API or _DB as suffix so we know what approach should be used.

Running
```sh
python3 parse.py
```
should now generate a .tsv file with found RSIDs and other datapoints for the specified KB. If you have a file from both the API and the Manual Search they can be merged by running:
```sh
python3 analytics.py
```
But first you have to comment out line 57. 
If you delete th comment sign in l.58 and 59 you will get information regarding the overlap of the two approaches.

You can now use these files and overlap different KBs with one another. Put the names of the files in the listofKBs and lenOfKBs list. After commenting out l.57-19 you can run:
```sh
python3 analytics.py
```
and it will overlap the KBs and give you the results. For depiction a Venn diagramm will be used if two KBs are in the list.
