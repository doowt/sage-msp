# Table of contents
1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
    1. [Complete Access Structures](#cas)
		1. [Standard CAS](#cas_standard)
		2. [Threshold CAS](#cas_threshold)
	2. [Monotone Span Programs](#msp)
		1. [Standard MSP](#msp_standard)
		2. [Replicated Secret-Sharing](#msp_replicated)
		3. [Shamir's Secret-Sharing](#msp_shamir)
		4. [Multiplicative MSP](#msp_mult)

<a name="introduction"></a>
# Sage-MSP 
This repository provides SageMath modules to allow creation and
manipulation of monotone span programs (MSPs) and complete access
structures (CASs).  The module is intended to be purely educational.

This package can:
* Create a CAS from an under-specified access structure (i.e. where
  not every set is declared as qualified or unqualified)
* Create an MSP representing replicated secret-sharing from any CAS
* Create an instance of Shamir's secret-sharing from a threshold
  access structure
* Compute the recombination vectors for each subset of parties for a
  given MSP (i.e. determine how to combine the shares of any qualified
  subset of parties to retrieve the secret)
* Compute an access structure from a given MSP
* Compute the error-checking matrix for `Q_2` access structures for
  use in Multi-Party Computation protocols such as those by Smart and
  Wood [[1](#ref_1)]
* Convert any MSP computing a `Q_2` access structure into a
  multiplicative MSP by at most doubling the total number of shares,
  i.e. where parties can compute an additive share of the product of
  two secrets by local computation, using the technique described by Cramer et al.
  [[2](#ref_2)]

<a name="prerequisites"></a>
## Prerequisites 

This package has been tested with SageMath 9.5 and Python 3.10.6.

<a name="installation"></a>
## Installation 

Clone this repository, start `sage-ipython` and run `load("msp.sage")` to load
the MSP and CAS classes, or `load("cas.sage")` to load just the CAS
class.  The `examples.sage` file provides some basic example usage.

<a name="usage"></a>
## Usage 

<a name="cas"></a>
### Complete Access Structures 

A complete access structure is an access structure where every subset
is either qualified or unqualified.  A complete access structure can
be derived given a set of maximally unqualified sets, or minimally
qualified sets.  See [[3](#ref_3)] for more detail.

The file `cas.sage` introduces a sage class `CAS` for complete access
structures.  An access structure is defined from a matrix `M` encoding
the sets, and a bit `b`.

```
my_cas = CAS(M, b)
```

If `b == 1`, then `M` represents the MUSs (maximally unqualified
sets).  Any set that is not a subset of any MUS specified is
considered to be qualified.

If `b == 0` then `M` represents the MQSs (minimally qualified sets).
Any set that is not a superset of any MQS specified is considered to
be unqualified.

Each set is represented by one row of the matrix, with a 1 in position
`(i,j)` iff the `i`'th MUS or MQS contains `P_j` (with party indexing
starting at 0).  For example, for the access structure defined as

```
delta_plus = {{0, 1}, {0, 2}, {0, 3}, {1, 2, 3}}
```

the access structure is declared as
```
my_cas = CAS(matrix([[1,1,0,0], [1,0,1,0], [1,0,0,1], [0,1,1,1]]), 1)
```

A subclass for threshold access structures is also defined:
```
my_threshold_cas = Threshold_CAS(n, t)
```
where `n` and `t` are integers, the number of parties and the threshold,
respectively.  (The threshold `t` for this class is defined as the smallest
number of parties in a set to be qualified; note that some authors choose this to be `t+1`.)

#### Objects, Attributes, and Methods
The Python help function can be used to determine the methods
associated with an object (e.g. run `help(CAS)`).

<a name="cas_standard"></a>
##### Standard CAS 
To initialise a CAS, run the following:

```
my_cas = CAS(M, b)
```

| Initialisation argument | Description                                                                                                      |
|-------------------------|------------------------------------------------------------------------------------------------------------------|
| `M`                     | A matrix `M = (m_ij)` whose rows represent sets of parties, where`m_ij = 1 ` iff party `P_j` is in the `i`th set |
| `b`                     | if `b == 1` then `M` represents unqualified sets, and if `b == 0` then `M` represents qualified sets             |

###### Attributes

| Variable Name | Description                                                                                   |
|---------------|-----------------------------------------------------------------------------------------------|
| `n`           | Number of parties, which is the number of columns of `M`                                      |
| `delta`       | Unqualified sets as a matrix `A = (a_ij)`: `a_ij == 1` iff the `i`th set contains party `P_j` |
| `delta_plus`  | Maximal elements of `delta`                                                                   |
| `gamma`       | Qualified sets as a matrix `A = (a_ij)`: `a_ij == 1` iff the `i`th set contains party `P_j`   |
| `gamma_minus` | Minimal elements of `gamma`                                                                   |

###### Methods

| Method Name   | Description                                                                    |
|---------------|--------------------------------------------------------------------------------|
| `print_all()` | Prints all information about an object                                         |
| `is_q(l)`     | Returns true iff the access structure is `Q_l` where `l` is a positive integer (as defined by Beaver and Wool [[4](#ref_4)]) |

<a name="cas_threshold"></a>
##### Threshold Complete Access Structure 
To initialise a Threshold CAS, run the following:

```
my_threshold_cas = Threshold_CAS(n, t)
```

| Initialisation argument | Description                                                     |
|-------------------------|-----------------------------------------------------------------|
| `n`                     | Number of parties                                               |
| `t`                     | Threshold, where any `t` parties or more can recover the secret |


<a name="msp"></a>
### Monotone Span Programs
The file `msp.sage` introduces a sage class for monotone span programs.

There are three subclasses:
* `Shamir_MSP`: The secret-sharing scheme corresponding to this MSP is due to Shamir
[[6](#ref_6)].
* `Replicated_MSP`: The secret-sharing scheme corresponding to this
  MSP is due to Ito et al. [[5](#ref_5)].
* `Mult_MSP`: This class is defined from an input MSP, and creates a
new multiplicative MSP from the input if the access structure is
`Q_2`.  The MSP returned is the same as the input MSP if the MSP is
multiplicative, and otherwise converts the input MSP to a
multiplicative MSP via the method due to Cramer et al. [[2](#ref_2)], which at most
doubles the size (i.e. total number of shares).

#### Objects, Attributes, and Methods
The Python help function can be used to determine the methods
associated with an object (e.g. run `help(MSP)`).

<a name="msp_standard"></a>
#### Standard MSP 
To initialise a standard MSP, run the following:

```
my_msp = MSP(AS, F)
```

| Initialisation Argument | Description                             |
|-------------------------|-----------------------------------------|
| `AS`                    | Access structure of type `CAS`          |
| `F`                     | A field (e.g. `QQ`) |


##### Attributes
| Variable Name              | Description                                                                                                                                                                                                                                                                                      |
|----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `M`                        | MSP matrix                                                                                                                                                                                                                                                                                       |
| `N`                        | Error-detection matrix; the cokernel of `M`.  Used for Multi-Party Computation with `Q_2` access structures as described by [[1](#ref_1)]                                                                                                                                                        |
| `target`                   | MSP target vector                                                                                                                                                                                                                                                                                |
| `rowmap`                   | MSP map of rows of M to parties                                                                                                                                                                                                                                                                  |
| `comms_indicator`          | Encodes which parties will be received from in recombination of the secret; used to determine the communication channels required in recombination                                                                                                                                               |
| `list_of_reconst_matrices` | For each `i`, provides a matrix `A` such that `A` times the subvector of shares held by `P_i` after shares are transferred in a round of communication is the fully reconstructed vector of shares.  Used for Multi-Party Computation with `Q_2` access structures as described by [[1](#ref_1)] |
| `recomb_channels`          | Matrix of communication, as determined by `comms_indicator` and as described by [[1](#ref_1)], in order to reconstruct the secret                                                                                                                                                                |
| `recomb_vectors`           | For each `i`, provides a (public) vector `lambda` so that after party `P_i` receives shares from other parties when a secret is opened to obtain a vector `v`, the secret can be computed by `P_i` as `lambda * v`                                                                               |
		  
##### Methods

| Method Name                  | Description                                                                                                                                                                                                                                                                                 |
|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `compute_cas()`              | Computes the access structure from an input MSP.  Initialise any access structure on the correct number of parties, the MSP matrix, the field, the rowmap and the target vectors, and the access structure will be recomputed, along with all other MSP properties.  (See `examples.sage`.) |
| `compute_all_msp_properties` | Initialises all MSP properties and sanitises the MSP when changes are made or MSP is set up manually by defining the matrix |
| `print_all()`                | Prints all information about an object                                                                                                                                                                                                                                                      |


<a name="msp_replicated"></a>
#### Replicated Secret-Sharing 
To initialise replicated secret-sharing, run the following:
```
my_replicated_msp = Replicated_MSP(AS, F)
```

| Initialisation Argument | Description                             |
|-------------------------|-----------------------------------------|
| `AS`                    | An access structure of type `CAS`       |
| `F`                     | A field (e.g. `QQ`) |


##### Additional Attributes

| Attribute Name    | Description                                                                                                                        |
|-------------------|------------------------------------------------------------------------------------------------------------------------------------|
| `parties_to_sets` | Vector encoding a map of parties to sets to specify which shares each party should open                                            |
| `sets_to_parties` | Vector encoding the map of sets to parties to specify which party opens each share (where shares correspond to MSP matrix columns) |

##### Additional Methods

| Method Name                 | Description                                                                                                                                                        |
|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `check_share_consistency()` | Replicated secret-sharing involves the same value being distributed to multiple parties, and this routine verifies that each share is the same wherever it is held |
| `check_share_correctness()` | Performs the error-detection step described by [[1](#ref_1)] |


<a name="msp_shamir"></a>
#### Shamir's Secret-Sharing 
To initialise Shamir's secret-sharing, run the following:
```
my_replicated_msp = Shamir_MSP(AS, F)
```

| Initialisation Argument | Description                                 |
|-------------------------|---------------------------------------------|
| `AS`                    | An access structure of type `Threshold_CAS` |
| `F`                     | A field (e.g. `QQ`)     |


No additional methods or attributes are provided.

<a name="msp_mult"></a>
#### Multiplicative MSP 
To initialise a multiplicative MSP, run the following:
```
my_multiplicative_msp = Mult_MSP(MSP)
```

| Initialisation Argument | Description                               |
|-------------------------|-------------------------------------------|
| `MSP`                   | An MSP computing a `Q_2` access structure |

##### Additional Attributes

| Attribute Name          | Description                                                                                                                                                                                          |
|-------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `list_of_mult_matrices` | For each `i`, provides a matrix `A` for which if the share of two secrets `a` and `b` owned by `P_i` are `v` and `w`, then `transpose(v) * A * w` is the additive summand of `a * b` owned by `P_i`. |
| `mult_AS` | Access Structure of the product of two secrets, i.e. which subsets of parties hold enough shares to compute the product of two secrets |

## References
[1] <a name="ref_1">[Error-Detecting in Monotone Span Programs with Application to Communication Efficient Multi-Party Computation](https://ia.cr/2018/467)</a>

[2] <a name="ref_2">[On codes, matroids and secure multi-party computation from linear secret sharing schemes](http://ia.cr/2004/245)</a>

[3] <a name="ref_3">[New Monotone Span Programs from Old](https://ia.cr/2004/282)</a>

[4] <a name="ref_4">[Quorum-Based Secure Multi-party Computation](https://link.springer.com/content/pdf/10.1007/BFb0054140.pdf)

[5] <a name="ref_5">[How to Share a Secret](https://web.mit.edu/6.857/OldStuff/Fall03/ref/Shamir-HowToShareASecret.pdf)</a>

[6] <a name="ref_6">[Multiple assignment scheme for sharing secret](https://dl.acm.org/doi/10.1007/BF02620229)</a>

## Authors
Tim Wood

