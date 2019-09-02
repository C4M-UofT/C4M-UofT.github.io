---
layout: page
title: Software Requirements
menu: Software Requirements
weight: 10
id: software
---

# Installation

Installation involves installing `Python 3`, setting up your local environment and testing that it works. The following steps walk you through this.

## Install Python

### MacOS / Linux

First, open your `terminal`.

On MacOS, you can do this by going to `Applications/Utilities/Terminal` or, using  Spotlight by pressing the `command + space` keys, and searching for "Terminal").

Copy & paste the following into your terminal and press the `return` key:

__MacOS__

```shell
$ curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output ~/miniconda.sh
```

__Linux__

```shell
$ curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --output ~/miniconda.sh
```

Once this has finished, copy & paste the following and hit the `return` key:

```shell
$ bash ~/miniconda.sh
```

Agree to the terms of service.

Press `return`, `space`, and then `return` keys.

Hit `return` when you are asked to accept the default install location.

Enter "yes" when prompted by _"Do you wish the installer to prepend the Miniconda3 install location to PATH in your ~/.bashrc?"_


### Windows

Download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe).

If you are running a 32-bit version of Windows, download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe). If you don't know, you can check which version you are running [here](https://support.microsoft.com/en-us/help/15056/windows-7-32-64-bit-faq).

Double click the `.exe` file to start the installation. Accept the terms of service, and leave all the default values when installing.

## Setup your environment

We need to **create** a virtual environment and install `jupyter`. Note, you only have to do this ONCE!

The instructions here are almost identical for MacOS / Linux and Windows.

- MacOS / Linux, you run all these commands from your **terminal**.
- On Windows, type **Anaconda** in the search box, choose **Anaconda Prompt** from the list. Run all these command from there.

First, lets create an environment called `C4M`

```shell
conda create -n C4M python=3.6
```

When conda asks you to proceed, type `y`

```
proceed ([y]/n)? y
```

Activate this environment that you just created with

```shell
source activate C4M
```

Then copy & paste the following to install `jupyter`

```shell
conda install -c conda-forge jupyterlab
```

> Answer `y` if prompted. This might take a little bit.

And finally, run

```shell
conda install nb_conda
```

Again, answer `y` if prompted.

## Using the jupyter notebooks

Everything in this class will happen through the `jupyter` notebooks. Notebooks are somewhere we can mix code and english, and run the code right in our browsers. Every time you wish to open a notebook, you need to

First, activate your environment

__MacOS / Linux__

```shell
$ source activate C4M
(C4M) $ # you should notice your command prompt change when the environment is active!
```

__Windows__

```shell
$ conda activate C4M
(C4M) $
```

Then, run `jupyter`

```shell
(C4M) $ jupyter lab
```

This will open a page in your browser. Use it to find `hello_world.ipynb` on your computer.

You can download the `hello_world.ipynb` notebook from <a href="notebooks/hello_world.ipynb">here</a>.

Follow the instructions in the notebook to make sure you installed everything correctly.
