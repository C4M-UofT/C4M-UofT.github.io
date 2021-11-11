---
layout: page
title: Software
menu: Software
weight: 10
id: software
---
If you would like to skip to the installation instructions, please click on the link that corresponds to your computer's operating system: [MacOS](#installation-macos--linux), [Linux](#installation-macos--linux), [Windows](#installation-windows).

## Python
Software can be written using a variety of programming languages, each with their own advantages and disadvantages. 
This course will teach you the basics of programming in [Python](https://www.python.org/){:target="_blank"}. 
Python is a popoular programming language for a variety of reasons; here are just a few:
* Beginner-friendly (e.g., simple syntax, object-oriented)
* Popular for data science
* Hundreds of libraries and frameworks for special functionality

Even though this program only covers one programming language, the concepts you learn will make it easier for you to pick up other languages you may encounter.

## Jupyter
When people think of programs, most people think of files with lines of code.
[Jupyter](https://jupyter.org/){:target="_blank"} notebooks are designed to make programming more interactive, allowing you to mix English text, Python code, and visualizations all in the same document.]

## PyCharm
In theory, you can write all of your programs in a basic text editor (e.g., Notepad) and then run them in your computer's console.
However, just like how Microsoft Word provides enhanced capabilities to document editing (e.g., spelling and grammar checking, buttons for formatting), there are many integrated development environments (IDE) specifically designed to help you write programs.
Some of the features provided by most IDEs include a debugging tool, code auto-complete, syntax checking, and package management.
During the later stages of the course, we will use [PyCharm](https://www.jetbrains.com/pycharm/){:target="_blank"} for Python development.

## Google Colaboratory (Google Colab)
[Google Colab](https://colab.research.google.com/notebooks/intro.ipynb){:target="_blank"} is a browser-based environment that allows you to work on Jupyter notebooks without any setup. 
We will start off the course using Google Colab so that you build confidence in your skills.

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