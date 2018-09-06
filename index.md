# Overview

Welcome to the webpage for the 2018 Computing for Medicine workshops.

## Course Staff

- Marzyeh Ghassemi
- John Giorgi
- Daniyal Liaqat

## Table of contents

1. [Installation](#installation)
2. [Schedule](#schedule)

# Installation

Installation involves installing `Python 3`, setting up your local environment and testing that it works. The following steps walk you through this.

## Install Python

### MacOS / Linux

First, open your `terminal`

> On MacOS, you can do this by going to `Applications/Utilities/Terminal` or, using  Spotlight by pressing the `command + space` keys, and searching for "Terminal").

Copy & paste the following into your terminal and press the `return` key:

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
```

Once this has finished, copy & paste the following and hit the `return` key:

```
bash ~/miniconda.sh
```

Agree to the terms of service.

> Press `return`, `space`, and then `return` keys.

Hit `return` when you are asked to accept the default install location.

Enter "yes" when prompted by _"Do you wish the installer to prepend the Miniconda3 install location to PATH in your ~/.bashrc?"_

### Windows

Download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe).

> If you are running a 32-bit version of Windows, download the installer from [here](https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86.exe). If you don't know, you can check which version you are running [here](https://support.microsoft.com/en-us/help/15056/windows-7-32-64-bit-faq).

Double click the `.exe` file to start the installation. Accept the terms of service, and leave all the default values when installing.

## Setup your environment

We need to **create** a virtual environment and install `jupyter`. Note, you only have to do this ONCE!

The instructions here are almost identical for MacOS / Linux and Windows.

- MacOS / Linux, you run all these commands from your **terminal**.
- On Windows, type **Anaconda** in the search box, choose **Anaconda Prompt** from the list. Run all these command from there.

First, create an environment called `C4M`

```bash
conda create -n C4M python=3.6
```

> When conda asks you to proceed, type `y`

Activate the environment that you just created with

```bash
source activate C4M
```

Then copy & paste the following to install `jupyter`

```bash
conda install -c conda-forge jupyterlab
```

> Answer `y` if prompted. This might take a little bit.

And finally, run

```bash
conda install nb_conda
```

> Again, answer `y` if prompted.

## Using the jupyter notebooks

Everything in this class will happen through the `jupyter` notebooks. Notebooks are somewhere we can mix code and english, and run the code right in our browsers. Every time you wish to open a notebook, you need to

**First**, activate your environment

```bash
source activate C4M
(C4M) # you should notice your command prompt change when the environment is active!
```

**Then**, run `jupyter`

```bash
(C4M) jupyter lab
```

This will open a page in your browser. Use it to find `hello_world.ipynb` on your computer.

> You will have to download the `hello_world` notebook from the [course website](https://c4m-uoft.github.io).

Follow the instructions in the notebook to make sure you installed everything correctly.

> Note: Once the notebook is open, you need to make sure you see `Python [conda env:C4M]` in the top right-hand corner. If not, click this area and select it from the dropdown.

# Schedule

## Phase 1 (Fall 2018)

### Level 1

Session | Date | Time | Location | Resources
------- | ---- | ---- | -------- | ---------
Session 1 | Wednesday September 19, 2018 | 9:30 – 12:30 pm | MS 3281 | <ul>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/software.shtml">Software installation</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/operators_builtins.html">Arithmetic operators, built-in functions</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/variables.html">Variables</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/strings.html">Intro to Strings</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/functions.html">Functions</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/converting_between_types.html">Converting Between Types</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/programs.html">Writing Programs</a></li>          </ul>
Session 2 | Wednesday October 3, 2018 | 9:30 – 12:30 pm | MS 3287 | <ul><li>Recap: writing programs</li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/booleans.html">Booleans</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/if.html"><code>if</code> statements</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/nesting_reusing_functions.html">Nesting, reusing functions</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/more_str_ops_methods.html">More String Operations, Methods</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/docstrings_help_dir.html">Docstrings, <code>help</code>, <code>dir</code></a></li>              </ul>
Session 3 | Wednesday October 17, 2018 | 9:30 – 12:30 pm | DSC Innovation Lab, Gerstein Library | <ul>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/for_over_str.html"><code>for</code> loops over <code>str</code></a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/lists.html">Lists, loops, and <code>range</code></a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/list_mutability.html">List methods and mutability</a></li>                     </ul>
Session 4 | Wednesday November 21, 2018 | 9:30 – 12:30 pm | MS 3281 | <ul><li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/parallel.html">Parallel lists and strings</a></li><li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/nested_lists_loops.html">Nested lists and loops</a></li><li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/while.html"><code>while</code> loops</a></li>            </ul>

### Level 2

Session | Date | Time | Location | Resources
------- | ---- | ---- | -------- | ---------
Session 1 | Wednesday October 17, 2018 | 1:00 – 4:00pm | DSC Innovation Lab, Gerstein Library | <ul>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/software.shtml">Software installation</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/operators_builtins.html">Arithmetic operators, built-in functions</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/variables.html">Variables</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/strings.html">Intro to Strings</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/more_str_ops_methods.html">More String Operations, Methods</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/functions.html">Functions</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/docstrings_help_dir.html">Docstrings, <code>help</code>, <code>dir</code></a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/converting_between_types.html">Converting Between Types</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session1/programs.html">Writing Programs</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/booleans.html">Booleans</a></li>            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session2/if.html"><code>if</code> statements</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/for_over_str.html"><code>for</code> loops over <code>str</code></a></li>          </ul>
Session 2 | Wednesday November 21, 2018 | 1:00 – 4:00pm | MS 3281 | <ul>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/lists.html">Lists, loops, and <code>range</code></a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session3/list_mutability.html">List methods and mutability</a></li>                            <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/parallel.html">Parallel lists and strings</a></li>            </ul>
Session 3 | Wednesday December 19, 2018 | 9:30 – 12:30 pm | MS 3281 | <ul><!--              <li>Practice building and modifying <code>list</code>s</li>-->              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/nested_lists_loops.html">Nested lists and loops</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level1/session4/while.html"><code>while</code> loops</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase1/level2/memory_model.html">Python Memory Model</a></li>               </ul>

## Phase 2 (Winter 2019)

Session | Date | Time | Location | Resources
------- | ---- | ---- | -------- | ---------
Session 1 | Wednesday, January 23, 2019 | 1:00 – 4:00 pm | MS 3281 |             <ul>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/dictionaries.html">Dictionaries</a></li>              <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/files.html">Files</a>                 (<a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/story.txt">story.txt</a>, <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/january06.txt">january06.txt</a>)              </li>              <li>Bonus material: <br/><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/testing_debugging.html">Testing and Debugging</a> (<a href="http://c4m.cdf.toronto.edu/cohort3/phase2/session1/broken_is_teenager.pdf">broken_is_teenager.pdf</a>)</li>   </ul>        <ul>            <li>Project 1 preparation exercises                  <ul>                    <li>Exercise Set 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                    <li>Exercise Set 2:                       <ul>                      <li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                      <li>Part 2 <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/prep_exercises/project1_exercise2_partb.pdf">handout</a> (submit on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a>)</li>                     </ul>                    <li>Exercise Set 3:                      <ul>                      <li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li><li>Part 2 <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/prep_exercises/project1_exercise3_partb.pdf">handout</a>, <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/prep_exercises/tester.py">tester.py</a>, and <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/prep_exercises/sym_data1.txt">sym_data1.txt</a> (submit on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a>)</li>                      </ul>                  </ul>              </li>            </ul> </ul>
Session 2 | Wednesday February 27, 2019 | 1:00 – 4:00 pm | MS 3281 | <ul>              <li>Project 1: Medical Document Retrieval</li>              <ul>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/project1_worksheet.pdf">Worksheet</a></li>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/C4MPhaseIIProject1.pdf">Handout</a></li>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/project1.py">project1.py</a> (starter code)</li>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/wikipages.zip">wikipages.zip</a> (data)</li>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project1/tester.py">tester.py</a> (basic tests)</li>              </ul>            </ul>
Session 3 | Wednesday March 20, 2019 | 1:00 – 4:00 pm | MS 3281 |
Session 4 | Wednesday March 27, 2019 | 1:00 – 4:00 pm | MS 3281 | <ul>              <li>Project 2: Human Mobility and Epidemic Modelling</li>              <ul>                <li><a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project2/C4MPhaseIIProject2.pdf">Handout</a> (sample data <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project2/cities.txt">cities.txt</a>)</li>                <li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a> </li>                <li>Part 2 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a>; on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a> (starter code: <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project2/dijkstra.py">dijkstra.py</a>, tester: <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project2/dijkstra_tester.py">dijkstra_tester.py</a>)                <li>Part 3 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                <li>Part 4 on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a> (starter code: <a href="http://c4m.cdf.toronto.edu/cohort3/phase2/project2/simulation.py">simulation.py</a>)              </li>            </ul>            <li>Video tutorials:</li>              <ul>                <li><a href = "https://www.youtube.com/watch?v=NvzktOyLhdM">Representing Travel Distances in Graphs</a> &mdash; useful for Part 1</li>                <li><a href = "https://www.youtube.com/watch?v=O2EKA8yIw0E">Dijkstra's Algorithm</a> &mdash; explains algorithm in Part 2</li>                <li><a href = "https://www.youtube.com/watch?v=qugrJ8t3Wzg">Simulation Intro</a>; <a href = "https://www.youtube.com/watch?v=z8pO3wexTP0">Implementing a Simulation with Python</a> (<a href = "http://c4m.cdf.toronto.edu/summerphase/hw/simulation.py">code</a>) &mdash; useful for 4</li>              </ul>            </ul>
Session 5 | Wednesday April 3, 2019 | 1:00 – 4:00 pm | MS 3281 |
Session 6 | Wednesday April 17, 2019 | 1:00 – 4:00 pm | DSC Innovation Lab, Gerstein Library |

## Phase 3 (Fall 2018 - Winter 2019)

Session | Speaker | Date | Time | Location | Resources
------- | ------- | ---- | ---- | -------- | ---------
Session 1 | Fanny Chevalier | Tuesday, October 9, 2018 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library | <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar6/C4M_seminar6_part1.pdf">Part 1 slides</a> <br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar6/C4M_seminar6_part2.pdf">Part 2 slides</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar6/C4MPhase3Project6.pdf">Project Handout</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar6/foodfacts.json">foodfacts.json</a>
Session 2 | Jared Simpson | Tuesday, October 16, 2018 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library |  <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar4/C4M_seminar4_part1.pdf">Part 1 slides</a><br>          <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar4/C4M_seminar4_part2.pdf">Part 2 slides</a><br>          <a href="https://play.library.utoronto.ca/BNGVy0H2eeqd">Video</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar4/Seminar4Project.pdf">Project Handout</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar4/C4MProject4.zip">Project Starter code and data</a>
Session 3 | Frank Rudzick | Tuesday, November 20, 2018 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library | <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar1/C4M_seminar1_part1.pdf">Part 1 slides</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar1/C4M_seminar1_part2.pdf">Part 2 slides</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar1/Seminar1Project.pdf"> Project Handout</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar1/Seminar1Project.zip">Project Starter code and data</a>
Session 4 | Chris J. McIntosh | Tuesday, February 12, 2019 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library | <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar2/C4M_seminar2_part1.pdf">Part 1 slides</a><br>          <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar2/C4M_seminar2_part2.pdf">Part 2 slides</a><br>          <a href="https://play.library.utoronto.ca/9lVvPztBVHO_">Video For Part 2 Only</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar2/Seminar2Project.pdf">Project Handout</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar2/C4MProject2.zip">Project Starter code and data</a>
Session 5 | Michael  Brudno | Tuesday, March 26, 2019 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library | <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/C4M_seminar5_part1.pdf">Part 1 slides</a> <br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/C4M_seminar5_part2.pdf">Part 2 slides</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/SeminarProject5.pdf">Project Handout</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/C4MProject5.zip">Project Starter code and data</a>             Demo code:<br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/demo/default_parameters.py">default_parameters.py</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/demo/search_recursive.py">search_recursive.py</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/demo/reverse_recursive.py">reverse_recursive.py</a><br>            <a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar5/demo/tree_recursive.py">tree_recursive.py</a><br>            <a href="https://play.library.utoronto.ca/uYvlpi3QRPwZ">Video</a><br>
Session 6 | Marzyeh Ghassemi | Tuesday, April 9, 2019 | 4:00 – 6:00 pm | DSC Innovation Lab, Gerstein Library
