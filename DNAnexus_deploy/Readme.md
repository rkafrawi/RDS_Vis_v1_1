<!-- dx-header -->
# rds_visualizer (DNAnexus Platform App)

RDS visualizing webapp. Accepts file uploads in .rds format and renders a series of user-query guided visualizations including Feature Plots, Dim Plots, Violin plots.

This is the source code for an app that runs on the DNAnexus Platform.
This Readme.md file is based on the DNAnexus documentation pertaining to running RStudio Shiny Server & Apps. For the official documentation on running RStudio Shiny Server & Apps, see
[https://documentation.dnanexus.com/](https://documentation.dnanexus.com/getting-started/developer-tutorials/web-app-let-tutorials/running-rstudio-shiny-server-and-apps).
<!-- /dx-header -->

## Instructions for Building the Web App and Deploying on DNAnexus Cloud Platform

Here I describe the step-by-step approach I took in order to build and deploy the RShiny web app to DNAnexus.

### 1. Logging in and running dx-app-wizard

Log into the DNAnexus platform using a token. Tokens can be generated under the user profile > API Tokens. Upon successful login, you will be prompted to select an existing project. When creating a token, make sure the token scope is set to all projects. Otherwise, you will only have VIEW level permissions when selecting a project. 

```
dx login --token <token>
# Enter your username, token
# Select the project where you want to deploy your app.
```

### 2. Start the app wizard to create a new app

The wizard will ask you questions to finalize the app's configuration.

```
dx-app-wizard

# Example answers:

App Name: rds_visualizer
Title []: RDS Visualizer
Summary []: RDS Visualizing RShiny app on DNAnexus
Version [0.0.1]:    # click <ENTER>
1st input name (<ENTER> to finish):    # click <ENTER>
1st output name (<ENTER> to finish):    # click <ENTER>
Timeout policy [48h]:    # click <ENTER>
Programming language: bash
Will this app need access to the Internet? [y/N]: y
Will this app need access to the parent project? [y/N]: y
Choose an instance type for your app [mem1_ssd1_v2_x4]:  # click <ENTER>

```

This results in a directory with the following layout:

rds_visualizer <br>
├── dxapp.json  <br>
├── src  <br>
&nbsp;│&emsp;└── rds_visualizer.sh  <br>
├── Readme.developer.md  <br>
├── Readme.md  <br>
├── resources  <br>
&nbsp;&emsp;└── test

Note that the 'App Name:' parameter is used for the name of the top folder in the directory as well as the name of the shell script in the src subfolder of the directory. 

### 3. Convert to web app

The [dxapp.json](./DNAnexus_deploy/dxapp.json) file contains the web app metadata. Within this file, find the following line:

```
    "version": "0.0.1",
    "inputSpec": [],
```

Modify it to the following to convert the app to a web app:

```
    "version": "0.0.1",
    "httpsApp": {
          "ports": [443],
          "shared_access": "VIEW"
      },
    "inputSpec": [],
```

### 4. Containerization

To build the web app, a docker image containing the required packages needs to be created using a rocker/shiny base. To do so, you can follow the logic of the following [Dockerfile](./DNAnexus_deploy/Dockerfile):

```
# base image: rocker/verse (with a specific version of R)
#   has R, RStudio, tidyverse, devtools, tex, and publishing-related packages
FROM rocker/shiny

MAINTAINER Rizky Kafrawi <rkafrawi@dfci.harvard.edu>

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgsl-dev \
    libxml2-dev \
    libglpk-dev \
    libcurl4-openssl-dev

RUN R -e "install.packages(c('shiny','shinydashboard','shinydashboardPlus', 'purrr', 'stringr', 'DT', 'shinyalert', 'ggplot2'))" 

RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

```
Note that the Seurat install is not done in line with the other required packages. For whatever reason, including Seurat in the following manner: install.packages(c(...,...,'Seurat')) will cause all the packages besides for Seurat to be installed. Without Seurat, the web app will not run. For organizational purposes, this Docker file and all associated images/tarball files will be kept in the `rds_visualizer/resources` subfolder.

Using the Dockerfile, you can then run the following commands in the CLI:

```
docker build -t <image_name> . 
docker save <image_name> | gzip > <image_name>.tar.gz
```

After building the tarball, make a new directory within the project to contain your docker images. Navigate to the new subdirectory and upload your tarball.

```
dx mkdir <Name>
dx cd <Name>
dx upload <image_name>.tar.gz
```

### 5. Create app code

The file that will be used to create the web app's code will be called [rds_visualizer.sh](./DNAnexus_deploy/rds_visualizer.sh). This script makes a directory in the selected DNAnexus project and pulls the source code from this repository for the build.

```
#!/bin/bash
# rds_visualizer 0.0.1
set -eux
main() {
 # get rds_visualizer app code
 mkdir rds_vizualizer
 url=https://raw.githubusercontent.com/rkafrawi/RDS_Vis_v1_1/main/Webapp_source
 wget -P rds_visualizer/ $url/DESCRIPTION $url/server.R $url/ui.R
 # pull and run Shiny Server docker image
 dx download $DX_PROJECT_CONTEXT_ID:/rk_shiny/rds_vis_maria.tar.gz
 docker load -i rds_vis_maria.tar.gz
 # attach our K-means app's folder as a volume
 docker run --rm -p 443:3838 -v $PWD/rds_visualizer:/srv/shiny-server/ rds_vis_maria
}
```

This shell script deviates from the one present in the official [documentation](https://documentation.dnanexus.com/getting-started/developer-tutorials/web-app-let-tutorials/running-rstudio-shiny-server-and-apps) by accessing the tarball the user uploaded to their tarball-containing-subdirectory to configure the web app environment instead of using the rocker/shiny base.

Note also that the `DESCRIPTION` file is required since the `DisplayMode` parameter is contained here and if you want to change the display to 'Showcase' mode, you would have to modify it within the  `DESCRIPTION` file.

### 6. Build and deploy

To build and deploy the app, navigate to the folder above 'rds_visualizer' and execute the following command:

```
dx build -f rds_visualizer
```
A successful build and deploy to the platform will result in the resulting applet ID being printed into the CLI:

```
{"id": "applet-GYjF85j07Xbqx1vYf1K59PqV"}
```
### 7. Start the app

Click on the app under the selected project to run the web app and click 'next.' Then, click on 'Start Analysis' at the top right-hand side of the DNAnexus website and click 'Launch Analysis' on the following popup.

## Parting Developer Notes


