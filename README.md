# VizAvalanche
VizAvalanche – Hacking some nice plots!

## Short introduction video

https://bbb00.mirz.uni-jena.de/playback/presentation/2.3/9b2ee95a5c54d7d9a07a7a6b637f0f22d7266968-1778505675697

## Overiview

VirJenDB (virjendb.org) is a large-scale virus genome sequence database containing approximately 15 million sequences, 910 thousand vOTUs and 90 metadata fields.
Currently statistics of these data can be found at https://virjendb.org/stats.
We aim to extend this page with informative and visually appealing summaries of the vOTU dataset.
The goal of this hackathon project is to create and add plots to dev.virjendb.org/stats.

The dataset underlying the existing statistics can be obtained from https://cloud.uni-jena.de/s/9BWwdkBCKYsi6yB for the duration of this hackathon.

The file `api.py` is the current code that creates most of these plots – (only for inspiration) no need to modify that. 

If you want to contribute to the plots or improve code:

* clone this repo
* eat some cookies
* implement your ideas
* share a pull request with us!

## general flow

* get a GitHub account
* clone this repository
* download the provided data
* IDE setup (+python environment)
    * we can help you with IntelliJ Idea
    * load the repository in you IDE
* start hacking!
* if you created something you consider useful: push it to your github account and create a pull request
  

## Ideas

### sequence length distribution

create a plot showing how abundant different sequence lengths are in our database

### submissions over time

create a plot showing a time scale and the number of sequences added at different time points

### number of clusters by GC content

create a plot of GC content bins and how many clusters or single sequences are in that bins

### logarithmic plots

There already are some plots on the page, that need a log scale.
We want to figure out, how to apply this.

### wold map with source data

So one idea is to create a world map, that visualizes information like:

* where came the samples from
* where are the institutions that submitted data to our DB

Steps to create such a map could be:

1. find existing implementations for such maps (plotly, javascript+svg, …)
   → there probably are libraries out in the wild that are prepared to create such maps
2. find out how to use these templates for our aims
3. figure out, how to read data from the CSV file or – if required – from out API
4. if any data is required that is not in the provided CSV file: specify which data we need to provide

## Getting data out of our API:

Our web page has an API, that could be used to gather additional data.
For you, the up-to-date API is available at https://api2.virjendb.org/v2

### Important endpoints:

* https://api2.virjendb.org/v2#/Metadata/_public_v2_metadata_public_get – inspect fields available to the public
* https://api2.virjendb.org/v2#/Download/_process_download_v2_download_post – download requird data to a file
* https://api2.virjendb.org/v2#/Search/_do_search_v2_search_post – search the database
