Renthal lab tutorial on scRNA-seq
================

**Authors**: Shams Bhuiyan and Parth Bhatia

## Description

This repository goes over the fundamentals for performing single-cell
and single nuclei RNA-seq analysis of the trigeminal and dorsal root
ganglia using the Seurat. It borrows heavily from the [Seurat
vignette](https://satijalab.org/seurat/articles/get_started.html) and
the [Harvard Bioinformatics Core
workshop](https://hbctraining.github.io/scRNA-seq_online/). After
completion of these tutorials, you should be able to see how I integrate
scRNA-seq of the trigeminal ganglia, the reasoning behind each
step/function, and the arguments that I think about for each
step/function. Some disclaimers:

-   This is really a skeleton of how a Seurat analysis can be done, and
    how you perform your analysis needs to make sense in light of your
    biological question, your tissue, the sequencing platform, whether
    you have cells or nuclei, number of cells/nuclei, species and
    additional variables that has not applied to me.
-   Seurat is just one of at least 16 tools that can be used to analyze
    the scRNA-seq data. It would be in your own interests to learn about
    the advantages and disadvantages of the different tools.
-   This tutorial is not meant to teach you computational and
    statistical concepts behind each step/function, simply how its done
    and a high level explanation of why its done. I will link to
    resources where I can for what I read and watch when I want to learn
    more technical details, but to go over each statistical and
    computational technique in depth would be challenging and too time
    consuming. And in some cases, I’m not the expert on the matter!
-   If you just want to see the code and not go through the tutorial, my scripts for anchoring and integration can be found in the Scripts folder

## Prerequisite skills

-   Basic understanding of R, R syntax, Rstudio, R data structures data wrangling
    with dplyr/tidyverse, and data visualization with ggplot2. The
    Harvard Bioinformatics Core has a great workshop
    [here](https://hbctraining.github.io/Intro-to-R-flipped/) and I also
    recommend the [Code Academy lessons](https://www.codecademy.com/).
-   Basic understanding of bash/shell scripting is required to learn how
    to store/access data properly, and submit jobs to a job scheduler.
-   Its not necessary, but it would be really cool if you learned how to
    use git for the sake of code documentation and reproducability.

## Installation requirements



## Computational resources

I typically deal with data objects that are larger than 5 Gb. As of
2023, my 2021 Macbook Pro with the M1 chip crashes when I use it to
analyze any scRNA-seq datasets. Additionally, storage of scRNA-seq
datasets become difficult on local machines. Fortunately we have access
to 2 computational clusters for computational **POWER** and DropBox for
storage

-   O2: O2 is the server from HMS and the cluster that I primarily use.
    O2 has only 15 GB of storage space on the home directory. it is not
    recommended to keep large data files on your home directory. O2 has
    two options for larger datasets. O2 also has a scratch directory,
    which is never meant for permanent storage. It actually WILL CLEAR
    any files that are unmodified for 30 days. The second option on O2
    is standby, for which can be accessed via a separate cluster called
    the transfer cluster. Individual labs can purchase 10 TB of storage
    on the standby directory. Standby is not necessarily meant for
    permanent backups and long term storage for previous projects, but
    rather as back ups for ongoing projects. To make an account, go
    [here](https://harvardmed.atlassian.net/wiki/spaces/O2/overview)

-   ERIS: Your home directory on Eris has a soft quota of 160 GB (they
    will probably issue a warning if it exceeds 160 GB) and a limit of
    400 GB. It is not generally recommended to fill up your storage
    space in the home directory but technically, you can store files
    there. The other option is the scratch directory but it is only used
    for temporary storage. Personal scratch directories can have 600 GB
    of data and the data WILL GET DELETED after 30 days of being
    un-modified, so I’d be very careful with storing data in the
    personal scratch directory. Data will not get deleted in your home
    directory. Get access to Eris
    [here](https://rc.partners.org/about/who-we-are-risc/enterprise-research-infrastructure-services)

-   Dropbox: What we use for permanent storage is a lab Dropbox account
    since BWH provides labs with unlimited storage if you request for a
    team folder. Individual users who are a part of MGB also get 1 TB of
    individual storage each. Fortunately, we have a team folder for the
    Renthal lab. See Will Renthal for access to it.

## Applications
Before you start this tutorial, you should

- Know how to [start Rstudio on O2](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/2233335809/HMS+-+RStudio+on+O2)
- Install Seurat. Use this command in Rstudio `install.packages("Seurat")`
