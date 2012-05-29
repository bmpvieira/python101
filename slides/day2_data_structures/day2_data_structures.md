Data Structures
===============

Bruno Vieira

Overview
========

* Introduction to operators
* What are variables and data structures?
* What are the types available?
* What are the methods availabe for each type?
* What purpose each has?
* Application examples?

Arithmetic operators
====================

* ``+`` Addition
* ``-`` Subtraction
* ``*`` Multiplication
* ``**`` Exponent
* ``/`` Division
* ``//`` Floor Division
* ``%`` Modulus

Arithmetic examples: + - * /
============================

```python
print 2 + 2
print 4 - 2
print 2 * 4
print 8 / 2
print 2 * 3 + 3 * 2
print (2 * 3) + (3 * 2)
print 2 * (3 + 3) * 2
```

Arithmetic examples: / //
=========================

```python
print 8 / 3
print 8.0 / 3
print 8.0 // 3
2<br />2.6666666666666665<br />2 #forceoutput
```

 Arithmetic examples: ** and operators style
============================================

```python
print 2 ** 4
print 2**4
print 2 **4
print 2** 4
```
[PEP 8 -- Style Guide for Python Code](http://www.python.org/dev/peps/pep-0008)

Arithmetic examples
===================

```python
print 8.0 / 3
print 8. / 3
```

Arithmetic examples
===================

```python
print 5.0 / 3
print 5.0 // 3
print 5.0 % 3
```

Variables
=========

"...a variable is a storage location and an associated symbolic name..." in [Wikipedia](http://en.wikipedia.org/wiki/Variable_(computer_science\))

Assignment statement
====================

```python
a = 1
print a
```

This is different from the equal operator ``==``

```python
a = 1
print a == 1
print a == 2
```

Variables usefulness
=====================

```python
print "Interesting genes number: "
print 15
print "Uninteresting genes number"
print 30
print "Total genes: "
print 15 + 30
```

Variables usefulness
=====================

```python
igenes = 15 # Interesting genes
ugenes = 30 # Uninteresting genes
print "Interesting genes number: "
print igenes
print "Uninteresting genes number: "
print ugenes
print "Total genes: "
print igenes + ugenes
```

Variables forbidden names
=========================
The following identifiers are used as reserved words, or keywords of the language, and cannot be used as ordinary identifiers.

    and       del       from      not       while
    as        elif      global    or        with
    assert    else      if        pass      yield
    break     except    import    print     None
    class     exec      in        raise     _*
    continue  finally   is        return    __*__
    def       for       lambda    try       __*

More information
----------------
[Python Documentation: Identifiers](http://docs.python.org/reference/lexical_analysis.html#identifiers)

Assignment operators
====================

* ``=`` Simple assignment operator
* ``+=`` Add AND assignment operator
* ``-=`` Subtract AND assignment operator
* ``*=`` Multiply AND assignment operator
* ``/=`` Divide AND assignment operator
* ``//=`` Floor Dividion and assigns a value
* ``%=`` Modulus AND assignment operator
* ``**=`` Exponent AND assignment operator

Assignment examples
===================

```python
a = 1
a = a + 3
print a
```

```python
a = 1
a += 3
print a
```

Assignment examples
===================

```python
counter = 0
\# Do something and increment counter
counter += 1
\# Repeat something in a loop and increment counter each time
counter += 1
\# Check expected counter number reached and stop working
print counter
```

Assignment examples
===================

```python
a = 2
a *= 20; print a
a /= 3; print a
a //= 2; print a
a %= 4; print a
a **= 6; print a
```

Data types
==========

* **Numbers (integers and floats)**
* **Strings**
* Lists
* Tuples
* Sets
* Dictionaries

More information
----------------
[Python Documentation](http://docs.python.org/library/stdtypes.html#numeric-types-int-float-long-complex)

Integers and floats
===================

* Numbers are created by numeric literals or as the result of built-in functions and operators
* Numeric literals containing a decimal point or an exponent sign yield floating point numbers
* Python fully supports mixed arithmetic

Integers and floats assignment
==============================

```python
a = 1
b = 2.0
print a + b
```

Strings
=======

* String literals are written in single or double quotes
* In triple-quoted strings, unescaped newlines and quotes are allowed (and are retained)
* Strings are immutable sequence types: such objects cannot be modified once created

Strings assignment
==================

```python
a = "hello"
b = "world"
print a + b
```

Strings assignment
==================
Please try the following lines in the editor below
--------------------------------------------------

    print "hello world"
    print 'hello world'
    print 'hello" "world'
    print "hello' 'world"
    print 'hello"
    print "world'

```python
print "hello world"
```

Strings assignment multiline
============================
    """ or '''

```python
seq = """ABCD
EFGH
IJKLMNOPQ"""
print seq
```

Integers + Strings
==================

```python
a = 1
b = "hello"
print a + b
```

Integers * Strings
==================

```python
a = 5
b = "hello "
print a * b
```

Strings slices
==============

    variable = "string"
    variable[start:stop:step]

```python
seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
print seq[1:5]
```

Strings slices
==============

```python
seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
print seq[:]
print seq[0]
print seq[0:5]
print seq[:5]
print seq[5:]
```

Strings slices with negative indices
====================================

```python
seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
print seq[-1]
print seq[-5:-1]
print seq[-0]
```

Strings slices with steps
=========================

```python
seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
print seq[0:10:2]
print seq[-1:-11:-2]
print seq[:-11:-1]
print seq[::-1] # This one is very useful
print seq[::-2]
```

Strings slices out of range
=========================

```python
seq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
print seq[5:100]
```

&nbsp;
===============

<h1>Strings methods</h1>



Strings methods
===============

capitalize center count decode encode endswith expandtabs find format index isalnum isalpha isdigit islower isspace istitle isupper join ljust lower lstrip partition replace rfind rindex rjust rpartition rsplit rstrip split splitlines startswith strip swapcase title translate upper zfill

Total: 38
---------

![CHALLENGE FAILED](https://raw.github.com/bmpvieira/python101/master/slides/day2_data_structures/assets/sad-challenge-failed.png "CHALLENGE FAILED")

Strings methods
===============

count endswith find islower isupper join lower lstrip replace rfind rsplit rstrip split splitlines startswith strip swapcase translate upper

Total: 19
---------

![CHALLENGE CONSIDERED](https://raw.github.com/bmpvieira/python101/master/slides/day2_data_structures/assets/determined-challenge-considered.png "CHALLENGE CONSIDERED")

Strings methods
===============

* count
* find
* replace
* startswith and endswith
* islower and isupper
* lower, upper and swapcase
* join
* strip
* translate
* split and splitlines &larr; later with *Lists*

Total: 10
---------

<div style="position: absolute; bottom: 10%; right: 20%">
![CHALLENGE ACCEPTED](https://raw.github.com/bmpvieira/python101/master/slides/day2_data_structures/assets/determined-challenge-accepted.png "CHALLENGE ACCEPTED")
</div>

Strings methods
===============
More information
----------------
[Python Documentation](http://docs.python.org/library/stdtypes.html#string-methods)

String count
============

    str.count(sub[, start[, end]])

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.count("A")
```

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.count("A", 0 , 10)
2 #forceoutput
```

String replace
==============

    str.replace(old, new[, count])

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.replace("T", "U")
```

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.replace("T", "U", 3)
TCCTGGUGGUGUATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA #forceoutput
```

String translate
================

    str.translate(table[, deletechars])
    string.maketrans(from, to)

```python
seq = "TCCTGGAGGAGAATGGAGGTCGAGANNNNNGGTGTGGTGGGGGTGANNNNN"
print seq.translate(None, "AT")
CCGGGGGGGGGCGGGCCGCGGGGGGGGGGGGGGGG #forceoutput
```

```python
import string # We'll learn about this line later
seq = "TCCTGGAGGAGAATGGAGGTCGAGANNNNNGGTGTGGTGGGGGTGANNNNN"
table = string.maketrans("AGTC", "1234")
print seq.translate(table, "N")
34432212212113221223421212232322322222321 #forceoutput
```

String startswith and endswith
==============================

    str.startswith(suffix[, start[, end]])

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.startswith("TCCT")
print seq.startswith("TGCT")
print seq.startswith("ATG", 12, 20)
True<br />False<br />True #forceoutput
```

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.endswith("TGA")
True #forceoutput
```

String find and rfind
=====================

    str.find(sub[, start[, end]])

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.find("GTC")
print seq.rfind("GTC")
```

```python
seq = "TCCTGGAGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
print seq.find("ATG", 0, 10)
-1 #forceoutput
```

String islower and isupper
==========================

    str.islower()

```python
a = "T"
cod1 = "ATG"
cod2 = "AtG"
print a.isupper()
print cod1.isupper()
print cod2.isupper()
True<br />True<br />False<br /> #forceoutput
```

String lower, upper and swapcase
==========================

    str.lower()

```python
cod = "AtG"
print cod.lower()
print cod.upper()
```

```python
cod = "AtG"
print cod.swapcase()
aTg #forceoutput
```

String join
===========

    str.join(iterable)

```python
str = "TCCTGG"
print ":".join(str)
```

String strip, lstrip and rstrip
===============================

    str.strip([chars])

```python
seq = "NNNNGCGCGCTGGAGGAGGTGAGAAGTTTAGGGTAAAAAAAAANNNN"
seq = seq.strip("N")
print seq
print seq.lstrip("GC")
print seq.rstrip("A")
print seq.strip("ATG")
GCGCGCTGGAGGAGGTGAGAAGTTTAGGGTAAAAAAAAA<br />TGGAGGAGGTGAGAAGTTTAGGGTAAAAAAAAA<br />GCGCGCTGGAGGAGGTGAGAAGTTTAGGGT<br />CGCGC #forceoutput
```

Data Structures
===============

* Numbers
* Strings
* **Lists**
* Tuples
* Sets
* Dicts

Lists
=====

* List of objects (Duh!)
* [Mutable sequence type (allow in-place modification of the object)](http://docs.python.org/library/stdtypes.html#mutable-sequence-types)


Lists assignment
================

    [a, b, c]

```python
a = []
b = [1,2,3,4,5]
c = ["a", "b", "c", "d", "f"]
print a
print b[0]
print c[-1]
```

Lists assignment
================

```python
a = [0] * 5
b = [1] * 5
c = ['NA'] * 5
d = range(5) # range([start], stop[, step]) [1]
print a
print b
print c
print d
```

#####[1] [Python Built-in Functions](http://docs.python.org/library/functions.html#built-in-functions)

Lists assignment
================

    str.split([sep[, maxsplit]])

```python
raw_data = "TCCTGGAGGAG;GTCAAGGGTCCAGCT;GGAGAAGTTTAGGG;TGTGGTG;GGGGTGA"
sequences = raw_data.split(";")
print sequences
```
```python
raw_data = "TCCTGGAGGAG;GTCAAGGGTCCAGCT;GGAGAAGTTTAGGG;TGTGGTG;GGGGTGA"
sequences = raw_data.split(";", 3)
['TCCTGGAGGAG', 'GTCAAGGGTCCAGCT', 'GGAGAAGTTTAGGG;TGTGGTG;GGGGTGA'] #forceoutput
```

Lists assignment
================

    str.splitlines([keepends])

```python
raw_data = """TCCTGGAGGAG
GTCAAGGGTCCAGCT
GGAGAAGTTTAGGG"""
sequences = raw_data.splitlines()
['TCCTGGAGGAG', 'GTCAAGGGTCCAGCT', 'GGAGAAGTTTAGGG'] #forceoutput
```

Lists assignment
================

    str.splitlines([keepends])

```python
raw_data = """TCCTGGAGGAG
GTCAAGGGTCCAGCT
GGAGAAGTTTAGGG"""
sequences = raw_data.splitlines(True)
['TCCTGGAGGAG\n', 'GTCAAGGGTCCAGCT\n', 'GGAGAAGTTTAGGG'] #forceoutput
```


Lists are a mutable sequence type
=================================

```python
mylist = ["A", "T", "G", "C"]
print mylist[0]
mylist[0] = "T"
print mylist
```

```python
myseq = "ATGC"
print myseq[0]
```

Multi-dimensional lists
=======================

```python
mylist = [["A", "T", "G", "C"], [1, 2, 3, 4]]
print mylist
print mylist[0]
print mylist[1][2]
```

Lists methods
=============
* count
* index
* append
* insert
* remove
* pop
* reverse
* sort

More information
----------------
[Python Documentation](http://docs.python.org/library/stdtypes.html#mutable-sequence-types)

Lists count, index, append and insert
=====================================

```python
mylist = ["A", "T", "G", "C", "A", "T"]
print mylist.count("A")
print mylist.index("C")
mylist.append("TAIL")
mylist.insert(3, "MIDDLE")
print mylist
```

Lists remove
============

```python
mylist = ["A", "T", "G", "C", "A", "T"]
mylist.remove("A")
print mylist
mylist.remove("A")
print mylist
```

Lists pop
===========

```python
mylist = ["A", "T", "G", "C", "A", "T"]
mylist.pop(2)
print mylist
print mylist.pop(1)
print mylist
```

Lists reverse and sort
======================

```python
mylist = [1, 1, 3, 5, 2, 4]
mylist.sort()
print mylist
mylist.reverse()
print mylist
```

Lists with Python Built-in Functions
====================================

```python
mylist=range(5)
print mylist
print len(mylist)
print min(mylist)
print max(mylist)
```

More information
----------------
[Python Documentation](http://docs.python.org/library/functions.html#built-in-functions)


Lists with join
===============

```python
mylist = ["A", "T", "G", "C", "A", "T"]
print "".join(mylist)
```

Data Structures
===============

* Numbers
* Strings
* Lists
* **Tuples**
* Sets
* Dicts

Tuples
======

* Lists that are immutable (like strings)
* Useful for storing heterogeneous data in which order has semantic value (like coordinates)
* Fast!!! 

More information
----------------
[Python Docs: Tuples and Sequences](http://docs.python.org/tutorial/datastructures.html#tuples-and-sequences)

Tuples assignment
=================

```python
coord1 = 12, 35 # pair of coordinates
coord2 = (32, 12)
coordinates = [coord1, coord2]
print coordinates 
```

Data Structures
===============

* Numbers
* Strings
* Lists
* Tuples
* **Sets**
* Dicts

Sets
====

* Unordered collection with no duplicate elements
* Uses include membership testing and eliminating duplicate entries
* Also support mathematical operations like union, intersection, difference, and symmetric difference.

More information
----------------
[Python Docs: Sets](http://docs.python.org/tutorial/datastructures.html#sets)

Sets
====

```python
palette = set(['blue', 'red', 'green', 'red'])
print "blue" in palette
print "magenta" in palette
```

Sets
====

```python
p1 = set(['blue', 'red', 'green', 'red'])
p2 = set(['yellow', 'green', 'blue', 'yellow', 'blue'])
print p1 - p2 # colors in p1 but not in p2
print p1 | p2 # colors in either p1 or p2
print p1 & p2 # colors in both p1 and p2
print p1 ^ p2 # colors in p1 or p2 but not both
set(['red'])<br />set(['blue', 'green', 'yellow', 'red'])<br />set(['blue', 'green'])<br />set(['red', 'yellow']) #forceoutput
```

Data Structures
===============

* Numbers
* Strings
* Lists
* Tuples
* Sets
* **Dicts**

Dictionaries
============

* Unordered set of "key: value" pairs, with the requirement that the keys are unique
* Known in other languages as "associative memories" or "associative arrays"
* Indexed by keys (unlike sequences, which are indexed by a range of numbers)
* Indices can be any immutable type (strings, numbers, or tuples of immutable objects)
* Usefull for storing, extracting or deleting values with a key

Dictionaries assignment
=======================

    {'a': 1, 'b': 2, 'c': 3}
    {'key1': "Value", 'key2': "Value", 1: "Another Value", 'd': 42}

```python
sequences = {'s1': "AGTAGCGT", 's2': "ATGAC", 'primer': "AGCTGCTAG"}
print sequences['primer']
del sequences['s2']
print sequences
```

Dictionaries assignment
=======================
```python
sequences = {'s1': "AGTAGCGT", 's2': "ATGAC", 'primer': "AGCTGCTAG"}
sequences['s1'] = "AAAAAAAA"
print sequences
print sequences.items()
print sequences.keys()
print sequences.values()
print 'primer' in sequences
```

Wrap up
=======

* Arithmetic operators ``+ - * / // % **``
* Assignment operators ``+= -= *= /= //= %= **=``
* Numbers and Strings ``a = 1; b = "Hello World"``
* String methods ``count, find, join, translate, split, etc``
* Lists and methods ``a = [1, 2]; append, pop, reverse, sort, etc``
* Some Built-in functions ``range, len, min, max``
* Tuples ``a = 1, 2, 3; b = (1, 2, 3)``
* Sets ``a = set([1, 2, 3])``
* Dictionaries ``a = {a: 1, b: 2, c: 3}``