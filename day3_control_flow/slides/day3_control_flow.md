Control flow
=================================
Francisco Pina-Martins

---

Day Overview:
=============

Today we will go around the basics of control flow:
----------------------------------------------
- What types are there?
- What does it do?
- How do I work with it?

---
Types:
======

There are essentially 2 types of flow control:
----------------------------------------------
- Conditionals
- Loops

---

Conditionals:
=============

"Boolean Operators"
-------------------
Let's take a short trip back to the land of _High School Mathematics_:

    <  » Less than;
    <= » Less or equal than;
    >  » Greater than;
    >= » Greater or equal than;
    != » Not equal to;
    <> » Not equal to (alternative);
    == » Equal to (Note the double "=");

"Boolean expressions"
---------------------

    x is y
    x is not y
    x in y
    x not in y

---

Conditional Statements:
=======================
When are they used?
---------------------
When we want our program to do different things [if](http://docs.python.org/reference/compound_stmts.html#if "if statement documentation") a determined condition is met.

How do they work?
-----------------
Let's look at some pseudo-code:

    if «Condition(s) to be met»:
        «Do something»
        «Do something else»
    elif «Another condition»:
        «Do something different»
    else:
        «Do something else enteirely»
    «Do something every time,»
    «Since this part is not indented»

Take special care with:
-----------------------

- You can have as many "elif"s as you wish;
- You can only have one "else" and it has to be after the last "elif";

---
Real code example:
==================
Let's look at a real code example:
----------------------------------

```python
sequence = "ATG"
if sequence == "ATG":
    print "We have a start codon!"
elif sequence in ["TGA", "TAG", "TAA"]:
    print  "We have a stop codon!"
else:
    print "Our sequence is neither a start nor a stop codon."
```

- Try changing the value of _sequence_ and see the different results.

Likewise, we can use other boolean operators:

```python
sequence = "ATG"
if len(sequence) == 3:
    print "This can be a codon."
elif len(sequence) > 3:
    print "This is too long for a codon."
else:
    print "This is too short for a codon."
```

- Notice the use of the [len()](http://docs.python.org/library/functions.html#len "len\(\) function documentation") function. It is used to return the length of an object, in this case, the length of the string _sequence_.

---

The _for_ loop:
===============
When are they used?
-------------------
When we want our program to do the same thing to a lot of things.
The [for loop](http://docs.python.org/reference/compound_stmts.html#for "for loop documentation") will do something __for__ every value in an [iterable](http://docs.python.org/glossary.html#term-iterable "iterable definition").

How do they work?
-----------------
Let's look at another pseudo-code example:

    for «item» in «iterable»:
        «Do somthing with item»
        «Do something else with item»
    «Do something after the loop is done»

Take special care with:
-----------------------

- An _iterable_ can be any iterable object, such as:
    - A string, a tuple a list or a dictionary;
        - Characters in a string;
        - Elements of lists and tuples;
        - Keys and values of dictionaries;
    - Integers and floats are __not__ iterable;
    - A _list_ of integers, however is iterable;

---

Real code example:
==================
Once again, let's look at a real code example:
----------------------------------------------

```python
for numbers in range(4):
    print(numbers)
```

- Running this code will print the numbers from 0 to 3 (remeber python starts to count from 0), each followed by a newline character.

- Also make note of the [range()](http://docs.python.org/library/functions.html#range "range\(\) function documentation")  function. It is used in this case to create a list of integers from 0 to 5 on the fly. It is a very versatile function, you can read more about it in the documentation.

- Another example could be:

```python
sequences = ["ATGCTAGCTGATC", "ATGCCCTGATTAT"]
for i in sequences:
    print(i)
```

Now that was easy, wasn't it? Let's make it a bit more difficult...

---

Nested loops:
=============

Sometimes we have some code that we want to run __x__ times and some code within that code that we want to run __y__ times.

- In this example we want to find which sequences are common to both lists:

```python
sequences1=["ATGTCTA", "TCGATCGA", "GCCCTAGT"]
sequences2=["ATCGCTA", "GCTATATT", "TCGATCGA"]
for i in sequences1:
    for j in sequences2:
        if i == j:
            print "Sequence %s is common to both lists" %(j)
```

Take special care with:
-----------------------
 - Nested loops can look like a good idea at first, but they usually have a great impact on performance. If you are working with large datasets, you are advised to avoid them.

---

The _while_ loop:
=================
When are they used?
-------------------
The [while](http://docs.python.org/reference/compound_stmts.html#while "while loop documentation") loop is used when we want to combine the functions of the _if_ statement and the _for_ loop (sort of).

How do they work?
-----------------
Here is some more pseudo-code as an example:

    while «Condition is true»:
        «Do something»
        «Do something else»
    «Do something after Condition is not true»

Take special care with:
-----------------------
- Make sure the contents of your _while_ loop alter the condition being verified, otherwise you may get caught in an "infinite loop".

---

Real code example:
==================
Let's look at another real code example:
----------------------------------------

```python
number=0
while number <= 3:
    print number
    number += 1
```

Running this code will yield the same result as our first _for_ loop, but it's done in a diffrent way.

As you can see, the _while_ loop will test against a condition and run the code in it while the condition is true.

Here's another example (a bit more bio and a bit less abstract). Let's call it an ORF generator:

```python
import random
ORF = "ATG"
bases = ["A","T","G","C"]
stops = ("TGA","TAG","TAA")
while ORF.endswith(stops) == False:
    ORF += random.choice(bases)
print ORF
```

Wow, wait a minuite, what is this? Let's look at it in parts. (Next slide please!)

---

The Mighty ORF generator:
=========================

    import random

This will import the functions from the _random_ module. Don't worry about it for now. We will have more fun with modules later.

Then, we declare our variables: _ORF_, _bases_ and _stops_, so far so good.

Finally the loop:

    while ORF.endswith(stops) == False:

What this means - "While the variable ORF does __not__ end with any of the content of _stops_ do this:"

    ORF += random.choice(bases)

What this means - "__Add__ a random character from _bases_ to ORF."

Here is the documentation for the used functions: [endswith()](http://docs.python.org/library/stdtypes.html#str.endswith "endswith\(\) documentation"), [random.choice()](http://docs.python.org/library/random.html#random.choice "random.choice\(\) documentation").

Can you see something wrong with this?

---

Deeper into control flow:
=========================

[Break](http://docs.python.org/reference/simple_stmts.html#break "break statement documentation"), [continue](http://docs.python.org/reference/simple_stmts.html#continue "continue statement documentation"), _else_ on loops and [pass](http://docs.python.org/reference/simple_stmts.html#pass "pass statement documentation"):
----------------------------------------------------------------------------------------------------------------

- Break
    - will immediately stop any _for_ or _while_ loop;
- Continue
    - will immediately continue with the next iteration of the loop;
- Else on loops
    - will do something _after and only_ the loop is finished. Breaking the loop will not run this code;
- Pass
    - will do absolutely nothing;

---

Real code examples:
===================
(We don't really need pseudo-code for this)
-------------------------------------------

```python
breakpoint = 4
skippoint = 2
for i in range(1,6):
    if i == skippoint:
        continue
    elif i == breakpoint:
        print("loop broke at " + str(breakpoint))
        break
    print i
else:
    print "loop never reached %s and never broke" %(breakpoint)
```

Take special care with:
-----------------------

- The way the [range()](http://docs.python.org/library/functions.html#range "range\(\) function documentation") function was used; In this case we also defined the _start_ of the count;
- The [str()](http://docs.python.org/library/functions.html#str "str\(\) function documentation") function - it will convert any object (in this case an _integer_) into a _string_. This is required to concatenate the variables in the [print()](http://docs.python.org/library/functions.html#print "print\(\) function documentation") function;
- Try to change the **breakpoint** and **skippoint** variables for different results;

---

Special type of iteration - dictionaries:
=========================================
- When "looping" through a dictionary, we can use a special function - [items()](http://docs.python.org/library/stdtypes.html#dict.items "dict.items() documentation")

```python
d = {"one":"1", "two":"2", "three":"3"}
for key,value in d.items():
    print key + " - " + value
```

What's so special about this?

Note that we are _iterating_ two variables at the same time. This can be tricky to master at first, but it is a very useful function once you've gotten the hang of it.

Take special care with:
-----------------------

- Dictionaries will not preserve the order that the _key:value_ pairs are stored in;
    - This means that when you iterate through a dictionary, your _key:value_ pairs can turn up in any order;
- You can do something similar with two (or more) lists by using the [zip()](http://docs.python.org/library/functions.html#zip "zip\(\) function documentation") function;

---

Biological examples:
====================

Let's suppose we have a dictionary of 3 lists with several species each and we wish to know in which of these lists (if at all) we can find our species - _Homo sapiens_

```python
listset = {"reptiles":["Lacerta lepida", "Psammodromus algirus",
"Aspidoscelis ironata"], "plants":["Arabidopsis thaliana", "Quercus suber",
"Vitis vinifera", "Ricinus comunis"], "mammals":["Mus musculus",
"Canis lupus", "Homo sapiens"]}
species = "Homo sapiens"
for lists in listset:
    if species in listset[lists]:
        print(species + " can be found in the following list: " + lists)
        break
else:
    print(species + " could not be found in any of the lists.")
```

Take special care with:
-----------------------

- Notice that when defining _listset_, the code is split along several lines; you can read more about this [here](http://www.python.org/dev/peps/pep-0008/#maximum-line-length "PEP8 style guide");
- In line 7, we are calling the values in the dictionary, not the keys;
- Try changing the variable _species_ and see the results;

---

Biological examples (part II):
==============================
In this example we have a string with 3 "columns" divided by tabs ("\t") in python. Let's suppose that we wish to extract the Fst value for each column into a list.

```python
datastring = """# Locus ID	Overall Pi	Fst
2	0.4	0.1666666667
3	0.5	0.0000000000
4	0.1	0.1095890411
5	0.2	0.2068965517"""
datalist = datastring.splitlines()
fsts = []
for lines in datalist:
    if lines.startswith("#"):
        pass
    else:
        values = lines.split("\t")
        fst = values[2]
        fsts.append(fst)
print(fsts)
['0.1666666667', '0.0000000000', '0.1095890411', '0.2068965517'] #forceoutput
```

Take special care with:
-----------------------

- The [splitlines()](docs.python.org/library/stdtypes.html#str.endswith "splitlines\(\) type documentation") type; this built-in will split a string into a list where each element is a line of the string;
- The [startswith()](http://docs.python.org/library/stdtypes.html#str.startswith "startswith\(\) function documentation") function; it is pretty much self explanatory;
- The [split()](http://docs.python.org/library/stdtypes.html#str.split "split\(\) function documentation") function; it will split a string into a list of words eliminating the separator.
- You have to test this in IDLE or equivalent.

---


