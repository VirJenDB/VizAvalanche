# VizAvalanche
VizAvalanche – Hacking some nice plots!

VirJenDB (virjendb.org) is a large-scale virus genome sequence database containing approximately 15 million sequences, 910 thousand vOTUs and 90 metadata fields.
Currently statistics of these data can be found at https://virjendb.org/stats.
We aim to extend this page with informative and visually appealing summaries of the vOTU dataset.
The goal of this hackathon project is to create and add plots to dev.virjendb.org/stats.

The dataset underlying the existing statistics can be obtained from https://cloud.uni-jena.de/s/9BWwdkBCKYsi6yB for the duration of this hackathon.

The file `api.py` is the current code that creates most of these plots.

If you want to contribute to the plots or improve code:

* clone this repo
* eat some cookies
* implement your ideas
* share a pull request with us!

## wold map with source data

So one idea is to create a world map, that visualizes information like:

* where came the samples from
* where are the institutions that submitted data to our DB

Steps to create such a map could be:

1. find existing implementations for such maps (plotly, javascript+svg, …)
2. find out how to use these templates for our aims
3. if any data is required that is not in the provided CSV file: specify which data we need to provide
4. alter the templates to show our data

## general flow

* get a GitHub account
* clone this repository
* download the provided data
* IDE setup (+python environment)
    * we can help you with IntelliJ Idea
    * load the repository in you IDE
* start hacking!
* if you created something you consider useful: push it to your github account and create a pull request
  
