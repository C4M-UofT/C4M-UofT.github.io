---
layout: page
title: Software Requirements
menu: Software Requirements
weight: 10
id: software
---

This course will be using teach you the basics of programming in `Python 3`. You can run `Python` code in many different ways, but we will be specifically using `Jupyter` notebooks. `Jupyter` notebooks can be run directly in your browser, allow you to mix written text with code, and have other benefits that we will discuss throughout the year.
The installation instructions vary slightly depending on your computer's operating system. Please click on the relevant link to jump the the correct instructions: [MacOS](#installation-macos--linux), [Linux](#installation-macos--linux), [Windows](#installation-windows).

## Installation (MacOS / Linux)

### Install Python 3
1. Open your `terminal`. On MacOS, you can do this by going to `Applications/Utilities/Terminal` or, using  Spotlight by pressing the `command + space` keys, and searching for "Terminal").
2. Copy & paste the following into your terminal and press the `return` key:
MacOS
```shell
$ curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output ~/miniconda.sh
```
Linux
```shell
$ curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --output ~/miniconda.sh
```
3. Once this has finished, copy & paste the following and hit the `return` key:
```shell
$ bash ~/miniconda.sh
```
4. Agree to the terms of service. Press `return`, `space`, and then `return` keys.
5. Hit `return` when you are asked to accept the default install location.
6. Enter "yes" when prompted by _"Do you wish the installer to prepend the Miniconda3 install location to PATH in your ~/.bashrc?"_

### Create a virtual environment
1. Run all of these commands from your **terminal**.
2. First, let's create an environment called `C4M`
```shell
conda create -n C4M python=3.6
```
3. When conda asks you to proceed, type `y`
```
proceed ([y]/n)? y
```

### Install Jupyter
1. Activate the virtual environment that you just created with
```shell
source activate C4M
```
2. Copy & paste the following to install `jupyter`
```shell
conda install -c conda-forge jupyterlab
```
Answer `y` if prompted. This might take a little bit.
3. Finally, run
```shell
conda install nb_conda
```
Again, answer `y` if prompted.

### Run Hello World
Every time you wish to open a notebook, you need to do the following:

1. Activate your environment
```shell
$ source activate C4M
(C4M) $ # you should notice your command prompt change when the environment is active!
```
This will open a page in your browser. 
2. Find `hello_world.ipynb` on your computer. You can download the `hello_world.ipynb` notebook from <a href="notebooks/hello_world.ipynb">here</a>.
3. Follow the instructions in the notebook to make sure you installed everything correctly.

## Installation (Windows)

### Install Python
1. Download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe). If you are running a 32-bit version of Windows, download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe). If you don't know, you can check which version you are running [here](https://support.microsoft.com/en-us/help/15056/windows-7-32-64-bit-faq).
2. Double click the `.exe` file to start the installation. Accept the terms of service, and leave all the default values when installing.

### Create a Virtual Environment Jupyter
1. Type **Anaconda** in the search box, choose **Anaconda Prompt** from the list. Run all these command from there.
2. First, let's create an environment called `C4M`
```shell
conda create -n C4M python=3.6
```
3. 
When conda asks you to proceed, type `y`
```
proceed ([y]/n)? y
```

### Install Jupyter
1. Activate this environment that you just created with
```shell
source activate C4M
```
2. Copy & paste the following to install `jupyter`
```shell
conda install -c conda-forge jupyterlab
```
Answer `y` if prompted. This might take a little bit.
3. Finally, run
```shell
conda install nb_conda
```
Again, answer `y` if prompted.

### Test Hello World
Every time you wish to open a notebook, you need to do the following:
1. First, activate your environment
```shell
$ conda activate C4M
(C4M) $
```
2. Then, run `jupyter`
```shell
(C4M) $ jupyter lab
```
This will open a page in your browser. 
3. Find `hello_world.ipynb` on your computer. You can download the `hello_world.ipynb` notebook from <a href="notebooks/hello_world.ipynb">here</a>.
4. Follow the instructions in the notebook to make sure you installed everything correctly.