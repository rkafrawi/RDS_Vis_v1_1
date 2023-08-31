<!-- dx-header -->
# rds_visualizer (DNAnexus Platform App)

RDS visualizing webapp. Accepts file uploads in .rds format and renders a series of user-query guided visualizations including Feature Plots, Dim Plots, Violin plots.

This is the source code for an app that runs on the DNAnexus Platform.
This Readme.md file is based on the DNAnexus documentation pertaining to running RStudio Shiny Server & Apps. For the official documentation on running RStudio Shiny Server & Apps, see
[https://documentation.dnanexus.com/](https://documentation.dnanexus.com/getting-started/developer-tutorials/web-app-let-tutorials/running-rstudio-shiny-server-and-apps).
<!-- /dx-header -->

<!-- Insert a description of your app here -->
The `DNAnexus_deploy/Dockerfile` was used to generate the rshiny image used to build the Docker image used in the web app. The `DNAnexus_deploy/dxapp.json` file contains the web app metadata. Lastly, the `DNAnexus_deploy/rds_visualizer.sh` script was used to build and deploy the web app.

## Instructions for Web App Build/Deploy
Here I describe the step-by-step approach I took in order to build and deploy the RShiny webapp to DNAnexus.

### 1.Logging in and running dx-app-wizard

Log into the DNAnexus platform using a token. Tokens can be generated under the user profile > API Tokens. The user will be prompted to select an existing project. When creating a token, make sure the token scope is set to all projects. Otherwise, the user will only have VIEW level permissions when selecting a project. 

```
dx login --token <token>
# Enter your username, token
# Select the project where you want to deploy your app.
```

### 2. Start the app wizard to create a new app:

The wizard will ask the user questions to finalize the app's configuration.

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

### 3. Create web app

The dxapp.json file contains the web app metadata. Within the dxapp.json file, find the following line:

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

### 4. Create Docker Image

To run the web app, a docker image containing the required packages needs to be created using a rocker/shiny base. To do so, users can follow the logic of the following [Dockerfile](./DNAnexus_deploy/Dockerfile):

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

Using the Dockerfile, users can then run the following commands in a terminal session:

```
docker build -t <image_name> . 
docker save <image_name> | gzip > <image_name>.tar.gz
```

### 5. Create app code
The file that will be used to create the web app's code will be called rds_visualizer.sh. Ths script makes a directory in the selected DNAnexus project.

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
### 6. Build and deploy

```
dx build -f rds_visualizer
```

### 7. Start the app
Click on the app under the selected project to run the web app and click 'next.' Then, click on 'Start Analysis' at the top right-hand side of the DNAnexus website and click 'Launch Analysis' on the following popup.

## Parting Developer Notes


