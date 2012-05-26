Data Structures
===============

Bruno Vieira

Overview
========

* Introduction to operators
* What are variables and data structures?
* What are the types available?
* What purpose each has?
* Application examples?

Arithmetic operators:
=====================

* ``+`` Addition
* ``-`` Subtraction
* ``*`` Multiplication
* ``**`` Exponent
* ``/`` Division
* ``//`` Floor Division
* ``%`` Modulus

Arithmetic examples
===================

```python
print 2 + 2
print 4 - 2
print 2 * 4
print 8 / 2
print 2 * 3 + 3 * 2
print (2 * 3) + (3 * 2)
print 2 * (3 + 3) * 2
```

 Arithmetic examples 
====================

```python
print 2 ** 4
print 2**4
print 2 **4
print 2** 4
```
[PEP 8 -- Style Guide for Python Code](http://www.python.org/dev/peps/pep-0008)

Arithmetic examples (Python 2 vs Python 3)
==========================================

Python 2
--------

```python
print 8 / 3
print 8.0 / 3
print 8.0 // 3
2<br />2.6666666666666665<br />2 #forceoutput
```
Python 3
--------

```python
print(8 / 3)
print(8.0 / 3)
print(8.0 // 3)
```

Python 2 vs Python 3
====================

More information:

[Python 2 or Python 3](http://wiki.python.org/moin/Python2orPython3)
[What’s New in Python](http://docs.python.org/py3k/whatsnew/index.html)

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

Assignment operators
====================

* ``=`` Simple assignment operator
* ``+=`` Add AND assignment operator
* ``-=`` Subtract AND assignment operator
* ``*=`` Multiply AND assignment operator
* ``/=`` Divide AND assignment operator
* ``%=`` Modulus AND assignment operator
* ``**=`` Exponent AND assignment operator
* ``//=`` Floor Dividion and assigns a value

Assignment examples
===================
```python
a = 1
a = a + 3
print a
a = 1
a += 3
print a
```

Assignment examples
===================

```python
a = 1
a += 3
print a
a -= 1
a -= 1
print a
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

Integers and floats
===================

```python
a = 1
b = 2
print a + b
```

Strings
=======
```python
a = "hello"
b = "world"
print a + b
```

Integers + Strings
==================

```python
a = 1
b = "world"
print a + b
```

Strings
=======

```python
seq = "TCCTGGAGGAGAATGGAGGTCAAGGGTCCAGCTGGAGAAGTTTAGGGTGTGGTGGGGGTGA"
```


Lists
=====

```python
a = [1,2,3,4,5]
b = ["a", "b", "c", "d", "f"]
print a[0]
print b[-1]
```

Tuples

# Dictionaries

# Sets

# Wrap up