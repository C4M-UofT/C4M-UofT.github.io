# Overview

Welcome to the webpage for the 2018 Computing for Medicine workshops.

# Course Staff

- Marzyeh Ghassemi
- John Giorgi
- Daniyal Liaqat

# Installation

Installation involves installing `Python 3`, setting up your local environment and testing that it works. The following steps walk you through this.

## Install Python

### MacOS / Linux

First, open your `terminal`.

> On MacOS, you can do this by going to `Applications/Utilities/Terminal` or, using  Spotlight by pressing the `command + space` keys, and searching for "Terminal").

Copy & paste the following into your terminal and press the `return` key:

__MacOS__

```bash
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh --output ~/miniconda.sh
```

__Linux__

```bash
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --output ~/miniconda.sh
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

First, lets create an environment called `C4M`

```bash
conda create -n C4M python=3.6
```

When conda asks you to proceed, type `y`

```
proceed ([y]/n)? y
```

Activate this environment that you just created with

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

First, activate your environment

__MacOS / Linux__

```bash
source activate C4M
(C4M) # you should notice your command prompt change when the environment is active!
```

__Windows__

```bash
conda activate C4M
(C4M)
```

Then, run `jupyter`

```bash
(C4M) jupyter lab
```

This will open a page in your browser. Use it to find `hello_world.ipynb` on your computer.

> You can download the `hello_world.ipynb` notebook from <a href="notebooks/hello_world.ipynb">here</a>.

Follow the instructions in the notebook to make sure you installed everything correctly.

# Schedule

## Phase 1 (Fall 2018)

### Level 1
<table>
	<thead>
    <tr>
      <th>Sessions</th>
      <th>Date</th>
      <th>Time</th>
      <th>Location</th>
      <th>Resources</th>
      <th>Homework</th>
    </tr>
	</thead>
    <tbody>
    <tr>
      <td>Session 1 </td>
      <td>Wednesday September 19, 2018</td>
      <td>9:30 – 12:30 pm</td>
      <td>MS 3281</td>
      <td>
      	<ul>
      	    <li>Arithmetic operators, built-in functions [<a href="examples/operators_builtins.html">html</a>, <a href="notebooks/operators_builtins.ipynb">notebook</a>]</li>
      	    <li>Variables [<a href="examples/variables.html">html</a>, <a href="notebooks/variables.ipynb">notebook</a>]</li>
      	    <li>Intro to Strings [<a href="examples/strings.html">html</a>, <a href="notebooks/strings.ipynb">notebook</a>]</li>
      	    <li>Functions [<a href="examples/functions.html">html</a>, <a href="notebooks/functions.ipynb">notebook</a>]</li>
      	    <li>Converting Between types [<a href="examples/converting_between_types.html">html</a>, <a href="notebooks/converting_between_types.ipynb">notebook</a>]</li>
      	    <li>Writing Programs [<a href="examples/programs.html">html</a>, <a href="notebooks/programs.ipynb">notebook</a>]</li>
      	</ul>
      </td>
      <td></td>
    </tr>
    <tr>
    	<td>Session 2</td>
    	<td>Wednesday October 3, 2018</td>
    	<td>9:30 – 12:30 pm</td>
    	<td>MS 3287</td>
    	<td>
    		<ul>
    			<li>Recap: writing programs</li>
    			<li>Booleans [<a href="examples/booleans.html">html</a>, <a href="notebooks/booleans.ipynb">notebook</a>]</li>
    			<li><code>if</code> statements [<a href="examples/if.html">html</a>, <a href="notebooks/if.ipynb">notebook</a>]</li>
    			<li>Nesting, reusing functions [<a href="examples/nesting_reusing_functions.html">html</a>, <a href="notebooks/nesting_reusing_functions.ipynb">notebook</a>]</li>
    			<li>More String Operations, Methods [<a href="examples/more_str_ops_methods.html">html</a>, <a href="notebooks/more_str_ops_methods.ipynb">notebook</a>]</li>
    			<li>Docstrings, <code>help</code>, <code>dir</code> [<a href="examples/docstrings_help_dir.html">html</a>, <a href="notebooks/docstrings_help_dir.ipynb">notebook</a>]</li>
    		</ul>
    	</td>
    	<td>
    		Phase I: Level 1 Session 2 Homework on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Oct. 18, 2018, 11:59 p.m</b>
    	</td>
    </tr>
    <tr>
		<td>Session 3</td>
		<td>Wednesday October 17, 2018</td>
		<td>1:00 – 4:00pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><code>for</code> loops over <code>str</code> [<a href="examples/for_over_str.html">html</a>, <a href="notebooks/for_over_str.ipynb">notebook</a></li>
				<li>Lists, loops, and <code>range</code> [<a href="examples/lists.html">html</a>, <a href="notebooks/lists.ipynb">notebook</a></li>
				<li>List methods and mutability [<a href="examples/list_mutability.html">html</a>, <a href="notebooks/list_mutability.ipynb">notebook</a>]</li>
			</ul>
		</td>
		<td>
    		Phase I: Level 1 Session 3 Homework on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Nov. 22, 2018, 11:59 p.m</b>
    	</td>
	</tr>
	<tr>
		<td>Session 4</td>
		<td>Wednesday November 21, 2018</td>
		<td>1:00 – 4:00pm</td>
		<td>MS 3281</td>
		<td>
			<ul>
				<li>Parallel lists and strings [<a href="examples/parallel.html">html</a>, <a href="notebooks/parallel.ipynb">notebook</a>]</li>
				<li>Nested lists and loops [<a href="examples/nested_lists_loops.html">html</a>, <a href="notebooks/nested_lists_loops.ipynb">notebook</a>]</li>
				<li><code>while</code> loops [<a href="examples/while.html">html</a>, <a href="notebooks/while.ipynb">notebook</a>]</li>
			</ul>
		</td>
		<td>
    		Phase I: Level 1 Session 4 In Class Assignment on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Nov. 22, 2018, 11:59 p.m</b>
    	</td>
    </tr>
  </tbody>
</table>

### Level 2
<table>
  <thead>
    <tr>
      <th>Sessions</th>
      <th>Date</th>
      <th>Time</th>
      <th>Location</th>
      <th>Resources</th>
      <th>Homework</th>
    </tr>
  </thead>
  <tbody>
  	<tr>
  		<td>Session 1</td>
  		<td>Wednesday October 17, 2018</td>
  		<td>1:00 – 4:00pm |</td>
  		<td>DSC Innovation Lab, Gerstein Library</td>
  		<td>
			PRE-WORK FOR LEVEL 2 STUDENTS:
  		 	<ul>
			    <li>Arithmetic operators, built-in functions [<a href="examples/operators_builtins.html">html</a>, <a href="notebooks/operators_builtins.ipynb">notebook</a>]</li>
			    <li>Variables [<a href="examples/variables.html">html</a>, <a href="notebooks/variables.ipynb">notebook</a>]</li>
			    <li>Intro to Strings [<a href="examples/strings.html">html</a>, <a href="notebooks/strings.ipynb">notebook</a>]</li>
			    <li>More String Operations, Methods [<a href="examples/more_str_ops_methods.html">html</a>, <a href="notebooks/more_str_ops_methods.ipynb">notebook</a>]</li>
			    <li>Functions [<a href="examples/functions.html">html</a>, <a href="notebooks/functions.ipynb">notebook</a>]</li>
			    <li>Docstrings, <code>help</code>, <code>dir</code> [<a href="examples/docstrings_help_dir.html">html</a>, <a href="notebooks/docstrings_help_dir.ipynb">notebook</a>]</li>
			    <li>Converting Between Types [<a href="examples/converting_between_types.html">html</a>, <a href="notebooks/converting_between_types.ipynb">notebook</a>]</li>
			    <li>Writing Programs [<a href="examples/programs.html">html</a>, <a href="notebooks/programs.ipynb">notebook</a>]</li>
			    <li>Booleans [<a href="examples/booleans.html">html</a>, <a href="notebooks/booleans.ipynb">notebook</a>]</li>
			    <li><code>if</code> statements [<a href="examples/if.html">html</a>, <a href="notebooks/if.ipynb">notebook</a>]</li>
			</ul>
					MATERIAL COVERED:
			<ul>
	  		    <li><code>for</code> loops over <code>str</code> [<a href="examples/for_over_str.html">html</a>, <a href="notebooks/for_over_str.ipynb">notebook</a>]</li>
			    <li>Lists, loops, and <code>range</code> [<a href="examples/lists.html">html</a>, <a href="notebooks/lists.ipynb">notebook</a>]</li>
			    <li>List methods and mutability [<a href="examples/list_mutability.html">html</a>, <a href="notebooks/list_mutability.ipynb">notebook</a>]</li>				
  		   	</ul>
  		</td>
  		<td>
    		Phase I: Level 2 Session 1 Homework on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Nov. 22, 2018, 11:59 p.m</b>
    	</td>
  	</tr>
	<tr>
		<td>Session 2</td>
		<td>Wednesday November 21, 2018</td>
		<td>1:00 – 4:00pm |</td>
		<td>MS 3281</td>
		<td>
		 	<ul>
			    <li>Parallel lists and strings [<a href="examples/parallel.html">html</a>, <a href="notebooks/parallel.ipynb">notebook</a>]</li>
			 	<li>Nested lists and loops [<a href="examples/nested_lists_loops.html">html</a>, <a href="notebooks/nested_lists_loops.ipynb">notebook</a>]</li>
			 	<li><code>while</code> loops [<a href="examples/while.html">html</a>, <a href="notebooks/while.ipynb">notebook</a>]</li>
		    </ul>
		</td>
		<td>
    		Phase I: Level 2 Session 2 Homework on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Dec. 20, 2018, 11:59 p.m</b>
    	</td>
	</tr>
	<tr>
		<td>Session 3</td>
		<td>Wednesday December 19, 2018</td>
		<td>9:30 – 12:30 pm</td>
		<td>MS 3281</td>
		<td>
			<ul>
			 	<li>Python Memory Model [<a href="examples/memory_model.html">html</a>]</li>
		 	</ul>
		</td>
		<td>
    		Phase I: Level 2 Session 3 Homework on <a href="https://pcrs.teach.cs.toronto.edu/C4M-2018-09">PCRS</a>
    		<br>Due date: <b>Dec. 27, 2018, 11:59 p.m</b>
    	</td>
	</tr>
  </tbody>
</table>

## Phase 2 (Winter 2019)
<table>
  <thead>
    <tr>
      <th>Sessions</th>
      <th>Date</th>
      <th>Time</th>
      <th>Location</th>
      <th>Resources</th>
    </tr>
  </thead>
  <tbody>
	<tr>
		<td>Session 1</td>
		<td>Wednesday January 23, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>MS 3281</td>
		<td>
			<ul>              
			<li><a href="examples/dictionaries.html">Dictionaries</a></li>              
			<li><a href="examples/files.html">Files</a> (<a href="examples/story.txt">story.txt</a>, <a href="examples/january06.txt">january06.txt</a>)</li>              
			<li>Bonus material:</li>
				<ul>
					<li><a href="examples/testing_debugging.html">Testing and Debugging</a></li>
					<li><a href="examples/broken_is_teenager.pdf">broken_is_teenager.pdf</a></li>
				</ul>
				<li>Project 1 preparation exercises</li>
					<ul>                    
						<li>Exercise Set 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                    
						<li>Exercise Set 2:</li>
							<ul>                      
								<li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                      
								<li>Part 2 <a href="projects/project1/project1_exercise2_partb.pdf">handout</a> (submit on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a>)</li>
							</ul>                   
					</ul>
				<li>Exercise Set 3:</li>
					<ul>                      
						<li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>
						<li>Part 2 <a href="projects/project1/project1_exercise3_partb.pdf">handout</a>, <a href="projects/project1/ex3_tester.py">tester.py</a>, and <a href="projects/project1/sym_data1.txt">sym_data1.txt</a> (submit on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a>)</li>
					</ul>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 2</td>
		<td>Wednesday February 27, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>MS 3281</td>
		<td>
			<ul>              
				<li>Project 1: Medical Document Retrieval</li>
					<ul>                
						<li><a href="projects/project1/project1_worksheet.pdf">Worksheet</a></li>                
						<li><a href="projects/project1/C4MPhaseIIProject1.pdf">Handout</a></li>                
						<li><a href="projects/project1/project1.py">project1.py</a> (starter code)</li>                
						<li><a href="projects/project1/wikipages.zip">wikipages.zip</a> (data)</li>                
						<li><a href="projects/project1/tester.py">tester.py</a> (basic tests)</li>
					</ul>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 3</td>
		<td>Wednesday March 20, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>MS 3281</td>
		<td> Work on projects. </td>
	</tr>
	<tr>
		<td>Session 4</td>
		<td>Wednesday March 27, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>MS 3281</td>
		<td>
			<ul>              
				<li>Project 2: Human Mobility and Epidemic Modelling</li>
				<ul>                
					<li><a href="projects/project2/C4MPhaseIIProject2.pdf">Handout</a> (sample data <a href="projects/project2/cities.txt">cities.txt</a>)</li>                
					<li>Part 1 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a> </li>                
					<li>Part 2 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a>; on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a> (starter code: <a href="projects/project2/dijkstra.py">dijkstra.py</a>, tester: <a href="projects/project2/dijkstra_tester.py">dijkstra_tester.py</a>)</li>
					<li>Part 3 on <a href="https://pcrs.teach.cs.toronto.edu/C4M17">PCRS</a></li>                
					<li>Part 4 on <a href="https://markus.teach.cs.toronto.edu/c4m-2017-09">MarkUs</a> (starter code: <a href="projects/project2/simulation.py">simulation.py</a>) </li>
				</ul>
			<li>Video tutorials:</li>
				<ul>                
					<li><a href = "https://www.youtube.com/watch?v=NvzktOyLhdM">Representing Travel Distances in Graphs</a> &mdash; useful for Part 1</li>                
					<li><a href = "https://www.youtube.com/watch?v=O2EKA8yIw0E">Dijkstra's Algorithm</a> &mdash; explains algorithm in Part 2</li>                
					<li><a href = "https://www.youtube.com/watch?v=qugrJ8t3Wzg">Simulation Intro</a>; <a href = "https://www.youtube.com/watch?v=z8pO3wexTP0">Implementing a Simulation with Python</a> (<a href = "http://c4m.cdf.toronto.edu/summerphase/hw/simulation.py">code</a>) &mdash; useful for 4</li>
				</ul>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 5</td>
		<td>Wednesday April 3, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>MS 3281</td>
		<td> Work on projects. </td>
	</tr>
	<tr>
		<td>Session 6</td>
		<td>Wednesday April 17, 2019</td>
		<td>1:00 – 4:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td> Work on projects. </td>
	</tr>
  </tbody>
</table>


## Phase 3 (Fall 2018 - Winter 2019)
<table>
  <thead>
    <tr>
      <th>Sessions</th>
      <th>Speaker</th>
      <th>Date</th>
      <th>Time</th>
      <th>Location</th>
      <th>Resources</th>
    </tr>
  </thead>
  <tbody>
	<tr>
		<td>Session 1</td>
		<td>Fanny Chevalier</td>
		<td>Tuesday, October 9, 2018</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><a href="seminars/C4M_seminar6_part1.pdf">Part 1 slides</a></li>
				<li><a href="seminars/C4M_seminar6_part2.pdf">Part 2 slides</a></li>
				<li><a href="seminars/C4MPhase3Project6.pdf">Project Handout</a></li>
				<li><a href="seminars/foodfacts.json">foodfacts.json</a></li>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 2</td>
		<td>Jared Simpson</td>
		<td>Tuesday, October 16, 2018</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><a href="seminars/C4M_seminar4_part1.pdf">Part 1 slides</a></li>
				<li><a href="seminars/C4M_seminar4_part2.pdf">Part 2 slides</a></li>
				<li><a href="https://play.library.utoronto.ca/BNGVy0H2eeqd">Video</a></li>
				<li><a href="seminars/Seminar4Project.pdf">Project Handout</a></li>
				<li><a href="seminars/C4MProject4.zip">Project Starter code and data</a></li>
			</ul>
		</td>
	</tr>
	<tr>
	<td>Session 3</td>
		<td>Frank Rudzick</td>
		<td>Tuesday, November 20, 2018</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><a href="seminars/C4M_seminar1_part1.pdf">Part 1 slides</a></li>
				<li><a href="seminars/C4M_seminar1_part2.pdf">Part 2 slides</a></li>
				<li><a href="seminars/Seminar1Project.pdf"> Project Handout</a></li>
				<li><a href="seminars/Seminar1Project.zip">Project Starter code and data</a></li>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 4</td>
		<td>Chris J. McIntosh</td>
		<td>Tuesday, February 12, 2019</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><a href="seminars/C4M_seminar2_part1.pdf">Part 1 slides</a></li>
				<li><a href="seminars/C4M_seminar2_part2.pdf">Part 2 slides</a></li>
				<li><a href="https://play.library.utoronto.ca/9lVvPztBVHO_">Video For Part 2 Only</a></li>
				<li><a href="seminars/Seminar2Project.pdf">Project Handout</a></li>
				<li><a href="http://c4m.cdf.toronto.edu/cohort2/phase3/seminar2/C4MProject2.zip">Project Starter code and data</a> </li>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 5</td>
		<td>Michael  Brudno</td>
		<td>Tuesday, March 26, 2019</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td>
			<ul>
				<li><a href="seminars/C4M_seminar5_part1.pdf">Part 1 slides</a> </li>
				<li><a href="seminars/C4M_seminar5_part2.pdf">Part 2 slides</a></li>
				<li><a href="seminars/SeminarProject5.pdf">Project Handout</a></li>
				<li><a href="seminars/C4MProject5.zip">Project Starter code and data</a> </li>
				<li>Demo code:</li>
					<ul>
						<li><a href="seminars/default_parameters.py">default_parameters.py</a></li>
						<li><a href="seminars/search_recursive.py">search_recursive.py</a></li>
						<li><a href="seminars/reverse_recursive.py">reverse_recursive.py</a></li>
						<li><a href="seminars/tree_recursive.py">tree_recursive.py</a></li>
					</ul>
				<li><a href="https://play.library.utoronto.ca/uYvlpi3QRPwZ">Video</a><br> </td></li>
			</ul>
		</td>
	</tr>
	<tr>
		<td>Session 6</td>
		<td>Marzyeh Ghassemi</td>
		<td>Tuesday, April 9, 2019</td>
		<td>4:00 – 6:00 pm</td>
		<td>DSC Innovation Lab, Gerstein Library</td>
		<td></td>
	</tr>
  </tbody>
</table>
