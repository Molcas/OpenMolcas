!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!  ********************************************************
!  ** Public-domain library routines used by casvb only. **
!  ********************************************************
!  **********************
!  ** EISPACK ROUTINES **
!  **********************

subroutine balanc(nm,n,a,low,igh,scale)
! this subroutine is a translation of the algol procedure balance,
! num. math. 13, 293-304(1969) by parlett and reinsch.
! handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
! this subroutine balances a real matrix and isolates
! eigenvalues whenever possible.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.
!
!    n is the order of the matrix.
!
!    a contains the input matrix to be balanced.
!
! on output
!
!    a contains the balanced matrix.
!
!    low and igh are two integers such that a(i,j)
!      is equal to zero if
!       (1) i is greater than j and
!       (2) j=1,...,low-1 or i=igh+1,...,n.
!
!    scale contains information determining the
!       permutations and scaling factors used.
!
! suppose that the principal submatrix in rows low through igh
! has been balanced, that p(j) denotes the index interchanged
! with j during the permutation step, and that the elements
! of the diagonal matrix used are denoted by d(i,j).  then
!    scale(j) = p(j),    for j = 1,...,low-1
!             = d(j,j),      j = low,...,igh
!             = p(j)         j = igh+1,...,n.
! the order in which the interchanges are made is n to igh+1,
! then 1 to low-1.
!
! note that 1 is returned for igh if igh is zero formally.
!
! the algol procedure exc contained in balance appears in
! balanc  in line.  (note that the algol roles of identifiers
! k,l have been reversed.)
!
! Questions and comments should be directed to Alan K. Cline,
! Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
! Electronic mail to cline@cs.utexas.edu.
!
! this version dated january 1989. (for the IBM 3090vf)
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

integer i, j, k, l, m, n, nm, igh, low
real*8 a(nm,n), scale(n)
real*8 c, f, g, r, s, b2, radix
logical first, noconv, skip

radix = 16.0d0

b2 = radix*radix
k = 1
l = n
! .......... search for rows isolating an eigenvalue and push them down ..........
first = .true.
do
  if (first) then
    first = .false.
  else
    if (l == 1) exit
    l = l-1
  end if
  ! .......... for j=l step -1 until 1 do -- ..........
  do j=l,1,-1

    skip = .false.
    do i=1,l
      if (i == j) cycle
      if (a(j,i) /= 0.0d0) then
        skip = .true.
        exit
      end if
    end do

    if (.not. skip) then
      m = l
      call rc_exchange()
      exit
    end if
  end do
end do

if (first .or. (l /= 1)) then

  ! .......... search for columns isolating an eigenvalue and push them left ..........
  first = .true.
  do
    if (first) then
      first = .false.
    else
      k = k+1
    end if

    do j=k,l

      skip = .false.
      do i=k,l
        if (i == j) cycle
        if (a(i,j) /= 0.0d0) then
          skip = .true.
          exit
        end if
      end do

      if (.not. skip) then
        m = k
        call rc_exchange()
        exit
      end if
    end do
  end do
  ! .......... now balance the submatrix in rows k to l ..........
  do i=k,l
    scale(i) = 1.0d0
  end do
  ! .......... iterative loop for norm reduction ..........
  do
    noconv = .false.

    do i=k,l
      c = 0.0d0
      r = 0.0d0

      do j=k,l
        if (j == i) cycle
        c = c+abs(a(j,i))
        r = r+abs(a(i,j))
      end do
      ! .......... guard against zero c or r due to underflow ..........
      if ((c == 0.0d0) .or. (r == 0.0d0)) cycle
      g = r/radix
      f = 1.0d0
      s = c+r
      do
        if (c >= g) exit
        f = f*radix
        c = c*b2
      end do
      g = r*radix
      do
        if (c < g) exit
        f = f/radix
        c = c/b2
      end do
      ! .......... now balance ..........
      if ((c+r)/f >= 0.95d0*s) cycle
      g = 1.0d0/f
      scale(i) = scale(i)*f
      noconv = .true.

      do j=k,n
        a(i,j) = a(i,j)*g
      end do

      do j=1,l
        a(j,i) = a(j,i)*f
      end do

    end do

    if (.not. noconv) exit
  end do

end if

low = k
igh = l

return

contains

! .......... in-line procedure for row and column exchange ..........
subroutine rc_exchange()
  integer i
  scale(m) = j
  if (j == m) return
  do i=1,l
    f = a(i,j)
    a(i,j) = a(i,m)
    a(i,m) = f
  end do
  do i=k,n
    f = a(j,i)
    a(j,i) = a(m,i)
    a(m,i) = f
  end do
end subroutine rc_exchange

end subroutine balanc
