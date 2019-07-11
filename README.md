Basic Polynomial Algebra Subprograms (BPAS)
===========================================

<p align="center"> 
<img src="http://www.bpaslib.org/bpas.png" alt="BPAS Logo">
</p>

The BPAS library provides support for arithmetic operations with polynomials
on modern computer architectures. Typical operations are
addition, multiplication, division, evaluation and interpolation.
The BPAS library also supports polynoial system solving
and real root isolation. The code is mainly written in C++
for ease of use, with underlying C for performance. The CilkPlus
extension is used in places for parallelism targeting multi-core processors.

Class Overview
--------------

Many different ring classes are supplied:
 - Integer
 - RationalNumber
 - ComplexRationalNumber
 - Fraction
 - SmartFraction
 - SmallPrimeField
 - BigPrimeField
 - GeneralizedFermatPrimeField

Many different polynomal classes are supplied:
 - DenseUnivariateIntegerPolynomial (DUZP)
 - DenseUnivariateRationalPolynomial (DUQP)
 - SparseMultivariateIntegerPolynomial (SMZP)
 - SparseMultivariateRationalPolynomial (SMQP)
 - SparseUnivariatePolynomial&lt;Ring&gt; (SUP)
 - DenseUnivariatePolynomial&lt;Field&gt; (DUP)
 - SmallPrimeFieldDistributedDenseMultivariateModularPolynomial (SDMP)
 - DistributedDenseMultivariateModularPolynomial (DDMP) 

For polynomial system solving, the TriangularSet and RegularChain classes are both templated
by a polynomial over a field :
 - TriangularSet&lt;Field, RecursivePoly&gt;
 - RegularChain&lt;Field, RecursivePoly&gt;



Library Structure
-----------------

The BPAS library makes use of abstract classes to define the inteface of many
common types. Classes which begin with "BPAS" are abstract. Concrete classes
then make use of multiple inheritance to satisfy many interfaces. 
For example, SMQP is both a BPASMultivariatePolynomial and a BPASGCDDomain. 
Indeed, these both further inherit from BPASRing, leading to [diamond inheritance](https://en.wikipedia.org/wiki/Multiple_inheritance#The_diamond_problem)
in SMQP.

It is possible that two BPASRing subclasses are *incompatible*. Consider
a RationalNumber and a DenseUnivariateIntegerPolynomial. Operators between these two rings 
is not well-defined. To restrict BPASRing subclasses to operate
polymorphically while also maintaining mathematical compatibility we make use 
of the ["Curiously Recurring Template Pattern"](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) (CTRP) and heavy use of templating. 

Moreover, it is advantageous to implement *conditional exporting* of particular functions.
That is, depending on the run-time characteritics of a particular instance, it should
provide certain funcionality. For example, consider SparseUnivariatePolynomial<Integer>
versus SparseUnivaraitePolynomial<RationalNumber>. The first is a BPASGCDDomain but the 
second is a BPASEuclideanDomain. To capture that the characteristics of a polynomial
change with its coefficeint ring, we implement *conditional exporting* using 
[introspection](https://en.cppreference.com/w/cpp/types/is_base_of) and 
[conditional inheritance](https://en.cppreference.com/w/cpp/types/conditional). 

With so much interesting inheritance and templating occuring within the classes
of BPAS, each class provides two versions of its class inheritance diagram. 
The first shows a simplified *semantic* inheritance diagram, free from excessive templating
(in particular, removing the CRTP) and providing a semantic over-view of the class
and its interactions with 
other classes, both contrete and abstract. The second inheritance diagram
provides the full picture of inheritance as it is actually programmed, providing
a full view of the templating and inheritance of a class.



Questions
---------
Questions, comments, bug reports: create an "issue" and we'll get back to you.
