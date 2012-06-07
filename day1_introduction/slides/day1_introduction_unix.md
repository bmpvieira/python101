Introduction
============
Python 101 team

---

# Course overview

- **Day 1**: Installing and configuring Python; Introduction to Python programming
- **Day 2**: Data structures
	- Strings
	- Lists
	- Tuples
	- Dictionaries
	- Sets
- **Day 3**: Control flow 
	- Conditional statements (if, elif)
	- Loops (while, for)
- **Day 4**: Functions and Input/Output
- **Day 5**: Modules (*sys*, *os*, *re*, *Bio*) and Do It Yourself day

---

# The UNIX shell

- The UNIX shell acts as an interface between the **user** and the **kernel**, which is the hub of the operating system

![unix shell](https://github.com/bmpvieira/python101/raw/master/day1_introduction/slides/assets/unix_shell.png)

- Through the shell, the user has a wide range of programs at its disposal that allow file manipulations, navigation through the directory structure, etc. Here we will provide a short tour through the shell using some of the most useful programs:
	- **ls**
	- **cd**
	- **touch**
	- **cp**
	- **mv**
	- **rm**
	- **mkdir**

---

# A tour through the shell

## Asking for help

<ul>
<li> Linux provides documentation and help for all of its programs:
<ul>
<li> using the **man** command before the program name

```python
man ls # Manual for the ls command
 #forceoutput
```
<li> using the **whatis** command before the program name

```python
whatis ls # Brief description of the ls command
 #forceoutput
```

<li> using the **--help**/**-h** option after the program name
</ul>

```python
ls --help # Provides information on the usage and options of the ls command
 #forceoutput
```
<ul/>

---

# A tour through the shell

## Navigating through the directory structure

- We will start in the home directory (~), which is familiar to most unix users. To view the contents of the home or any other directory, we can use the **ls** command

```python
ls # Shows the contents inside the current dir
ls ~/ # Shows contents inside the home dir
ls path/to/dir # Shows contents inside the specified dir
 #forceoutput
```

- The **ls** command has some useful options:
	- **-a**: Shows hidden files
	- **-d**: Lists directories instead of contents
	- **-l**: Long listing format, with aditional details

---

# A tour through the shell

## Navigating through the directory structure

- To change our current directory from home, we can use the **cd** command

```python
cd path/to/dir
cd # Returns you to the home directory
cd ~/ # Equivalent to "cd"
cd ./ # Stays on the current directory
cd ../ # Goes to the parent directory
 #forceoutput
```

- If you ever get lost and wish to know your location you can use the **pwd** command

```python
pwd # Returns the current directory
 #forceoutput
```

---

# A tour through the shell

## Creating new directories

- We can create a new directory for the material of this course using the **mkdir** command

```python
mkdir ~/Python101 # Creates the Python101 directory in the home dir
mkdir Python101 # Creates the Python101 directory in the current dir
 #forceoutput
```

## Manipulating files

- To create new files in which python scripts can be written, use the **touch** command

```python
touch my_script.py # Creates a new file name my_script.py in the current dir
 #forceoutput
```

- To copy or move/rename files or directories use the **cp** and **mv** commands, respectively

```python
cp my_script.py /usr/local/bin # copies the my_script.py file in the current dir to the /usr/local/bin/ dir
cp -r my_folder/ ~/ # Copies my_folder/ dir recursively to the HOME dir
mv my_script.py ~/ # Moves the my_script.py file to the home dir
mv my_script.py my_program.py # Renames my_script.py to my_program.py
 #forceoutput
```

---

# A tour through the shell

## Manipulating files

- To remove files or directories, use the **rm** command

```python
rm my_script.py # Removes my_script.py in the current dir
rm ~/my_script.py # Removes my_script.py from the HOME dir from any dir
rm -r my_folder/ # Removes my_folder/ in the current dir recursively
 #forceoutput
```

- To change the file owner and group, use the **chown** command

```python
chown USER FILENAME # Changes the user owner of the FILENAME
chown USER:GROUP FILENAME # Changes the user and group ownership of the FILENAME
chown :GROUP FILENAME # Changes the group owner of the FILENAME
 #forceoutput
```

---

# A tour through the shell

## Manipulating files

<ul>
<li> To change file access permissions, use the **chmod** command

```python
chmod MODE FILENAME
\# Symbolic mode
chmod -rwx FILENAME # Removes ("-" symbol) read ("r"), write ("w") and executable ("x") permissions of the FILENAME for the owner, group and others
chmod +rwx FILENAME # Gives ("+" symbol) read, write and executable permissions to FILENAME for the owner, group and others
\# Numerical mode
chmod 755 FILENAME # Changes FILENAME permisions to read, write and executable for the owner ("7") and read and executable for the group and others
 #forceoutput
```
<li> **Numerical mode rational**: File permissions are represented by a three-digit octal number. These numbers are added accordingly, using these values:
<ul>
	<li> 4 = read (r)
	<li> 2 = write (w)
	<li> 1 = execute (x)
	<li> 0 = no permissions (-)
</ul>
</ul>


---

# A tour through the shell

## Using the terminal screen
<ul>
<li> To clear the terminal window from previous commands so that the output of subsequent commands can be more clearly visualized, use the **clear** command

```python
clear
 #forceoutput
```

<li> To display the contents of a file on the screen, use the **cat** command

```python
cat FILENAME
 #forceoutput
```
<li> To display the contents of a file on the screen one page at a time, use the **less** command

```python
less FILENAME # Use the keywords "space" and "b" to move forward and backwards, respectively
 #forceoutput
```
</ul>

---

# A tour through the shell

## Using the terminal screen

- The initial or final contents of a file can also be visualized on screen using the **head** and **tail** commands

```python
head FILENAME # prints the first 10 lines of the FILENAME
head -n 50 # The -n option changes the lines read to 50
less FILENAME # prints the last 10 lines of the FILENAME
less -n 30 # Same as the -n options for the head command
 #forceoutput
```

## Filtering the screen output

- The output that is printed on the terminal screen can be filtered/searched using the **grep** command

```python
grep [options] PATERN FILENAME # grep searches the FILENAME for a given PATTERN and prints lines with that pattern on the screen
ls | grep PATERN # grep can be used with the pipe (|) to filter the output of other shell commands
 #forceoutput
```

---

# A tour through the shell

## Redirection and piping

- To redirect output the output of any shell command, use the **>** or **>>** symbols

```python
cat FILENAME > new_file.txt # The output of the cat command is redirected to a new file, instead of the terminal screen. The ">" symbol overwrites any contents in the original new_file.txt
grep PATERN FILENAME >> new_file.txt # The output of the grep command is redirected to a new file. The ">>" symbol appends the output to the end of the new_file.txt
ls | grep PATERN > new_file.txt 
 #forceoutput
```

- Pipes ("**|**") allow separate processes to communicate implicitly, enabling narrow tools to be combined in complex ways

```python
ls | grep PATERN # allows communication between "ls" and "grep". Searches for PATERN using grep, from the output of ls
cat FILENAME | wc # The FILENAME output is fed to the "wc" command, which counts the no. of words in a file
 #forceoutput
```

---
# What is python?

Python is an **interpreted**, **object-oriented**, **high-level** programming language, that emphasizes in code readability:

- **Interpreted**: Programs can be run as soon as they are written, no need to compile;
- **Object-oriented**: The programming paradigm, where "objects" can be somewhat separated from the rest of the program. An object is an entity that has a: 
	- **Identity**: Unique identifier of the object
	- **Value**: 
	- **Type**: whether it is a string, number, etc. 
- **High-level**: Abstracts from the computer it is being run on;

---

# Two ways of running python code

## The interactive mode

- Ideal for **training** and **testing blocks of code** on the fly. During an interactive session, instructions are directly executed. However, the code cannot be saved between different sessions. Some interactive interpreters include:
	- [Python's built-in interpreter](http://docs.python.org/tutorial/interpreter.html)
	- [Python IDLE](http://docs.python.org/library/idle.html)
	- [Dreampie](http://dreampie.sourceforge.net/)

## Script mode

- Used to construct Python programs following this simple procedure:
	- Create and edit a text file with Python code using a text editor ([geany](http://www.geany.org/), [vim](http://www.vim.org/)).
	- Save the file with the filename extension "**.py**"
	- Run the script with the Python interpreter

---

# The interactive mode

In its easiest form, the iteractive mode can be initialized by simply typing "**python**" in the command line prompt.

![Python interpreter](https://github.com/bmpvieira/python101/raw/master/day1_introduction/slides/assets/python_interpreter.png)

- The symbols "**>>>**" are the prompt of the interactive mode. To execute an instruction, write on the prompt line and press enter.
- The symbols "**...**" are the secondary prompt. They represent continuation lines, i.e., when the instruction is not yet complete and ready to be executed.
- To exit the interactive mode, you can use the keywords "**Ctrl + D**" or type "**quit()**"

---

# The script mode

Writing scripts will make the majority of your coding experience. A python script that contains a set of instructions can be executed in more than one way:

- **Shebang way**:
	1. The first line of the script file must contain a "shebang" (#!) followed by the path of the Python executable (e.g. #!/usr/bin/python). This will tell the program loader which interperter should be used to execute the file.
	1. Make the script file executable (usually with the "chomd +x filename" command)
	1. Run the script by typing "./my_script.py" in the terminal prompt.
- **Script argument way**:
	1. Run your python interpreter with the script as an argument (e.g. python my_script.py)

---

# Introducing Python programming

## Data structures (Day 2)

- Basic operators  ("+","-","/","\*")
- Variable assignment 
- Standard data structures (Strings, numbers, lists and dictionary)

## Control flow (Day 3)

- Control flow is used in Python whenever you want to:
	- **Make choices**, using *conditional statements* when conditions are met
	- **Perform repeated actions**, using a set of statements executed *n* times in a *loop* until a condition is met

---

# Introducing Python programming

## Functions (Day 4)

- Blocks of code that perform specific tasks (procedures) repeatedly
- How to create new functions
- Using arguments to bring flexibility to those procedures


## Input/Output (Day 4)

- The three input/output channels: standard input, standard output and standard error
- Handling input files (opening and parsing) and output files (writing)
- Interactive scripts that request user input
- Printing messages to the terminal

---

# Introducing Python programming

## Modules (Day 5)

- What are modules, how do they work and how to use them
- Overview of some of the most useful modules:
	- *re*, regular expressions
	- *sys*, system tools
	- *os*, operative system tools and commands
	- *biopython*, the bioinformatic's module of choice


---

# Online resources

Some of the most usefull resources that will help you solving problems and getting new ideas during your coding sessions:

- [Official Python tutorial](http://docs.python.org/tutorial/): From the official Python documentation, this resource covers all basic and advance built-in features of Python.
- [Python wiki](http://wiki.python.org/moin/): Includes several oficial Python resources, such as the [beginners guide](http://wiki.python.org/moin/BeginnersGuide), [beginners errors](http://wiki.python.org/moin/BeginnerErrorsWithPythonProgramming) and [Python books](http://wiki.python.org/moin/PythonBooks). 
- [Stack Overflow](http://stackoverflow.com/): An extremely usefull language-independent collaboratively edited question and answer site for programmers. Most of the questions and problems that you may face when coding have been probably already answered in this forum.

## [Google](google.com) is your friend

- Giving the current flood of digital data, the importance of google's search engine to sort over all information cannot be overstated. The trick, however, is knowing how to put your problem in order to obtain the most usefull answers. 

---

