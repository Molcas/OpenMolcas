   WORK IN PROGRESS

Utilities for handling sparse matrices

Sparse matrices are stored in a Modified Compressed Sparse Row format.
A matrix A (nxn) is stored in two one dimensional vectors of length (m): V and IJ
The vector V contains the non-zero elements of A, the vector IJ contains
some indices (integers). Such that:

- The first (n) elements of V are the diagonal elements of A, even if they are zeros.
- The element (n+1) of V is undefined (see below)
- The elements from (n+2) to (m) of V are the non-zero off-diagonal elements of A,
  in row order: first for row 1, then row 2, etc., and sorted by column

- The first (n) elements of IJ contain the indices for finding where in V are
  stored the elements for each row of A. Element (i) of IJ is the index of V where
  the first element of row (i) of A is found. If there are elements of A for this row,
  element (i) of IJ is the same as element (i+1). Element (1) of IJ is always n+2,
  element (n+1) of IJ is always m+1.
- Elements from (n+2) to (m) of IJ are the colums to which the corresponding element
  of V belongs (in A).

- Optionally, for symmetric A, only the lower triangle is stored, and then element
  (n+1) of V is set to 1.0. Otherwise (non-symmetric storage), it is set to 0.0.

Thus, m=n+k+1, where k is the number of non-zero off-diagonal elements of A.


EXAMPLE (n=6, k=9, m=16)

  i    1  2  3  4  5  6     j

A = [ 11 12 13 14  0  0 ]   1
    [  0 22 23  0  0  0 ]   2
    [  0  0 33 34 35 36 ]   3
    [  0  0  0 44 45  0 ]   4
    [  0  0 53  0  0  0 ]   5
    [  0  0  0  0  0 66 ]   6


   i    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16

 V = [ 11 22 33 44  0 66  * 12 13 14 23 34 35 36 45 53 ]
IJ = [  8 11 12 15 16 17 17  2  3  4  3  4  5  6  5  3 ]
(element * will have to be 0.0, because the matrix is not symmetric)


Storing the matrix per rows makes it easy to compute the product with a vector on the
right, but it's trickier to compute the transpose. Matrix products are more difficult
because the result does not have the same sparsity pattern.
