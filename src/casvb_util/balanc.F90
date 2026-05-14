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

subroutine balanc(nm,n,a,low,igh,scl)
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
!    scl contains information determining the
!       permutations and scaling factors used.
!
! suppose that the principal submatrix in rows low through igh
! has been balanced, that p(j) denotes the index interchanged
! with j during the permutation step, and that the elements
! of the diagonal matrix used are denoted by d(i,j).  then
!    scl(j) = p(j),    for j = 1,...,low-1
!           = d(j,j),      j = low,...,igh
!           = p(j)         j = igh+1,...,n.
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

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nm, n
real(kind=wp), intent(inout) :: a(nm,n)
integer(kind=iwp), intent(out) :: low, igh
real(kind=wp), intent(out) :: scl(n)
integer(kind=iwp) :: i, j, k, l, m
real(kind=wp) :: b2, c, f, g, r, s
logical(kind=iwp) :: first, noconv, skip
real(kind=wp), parameter :: rdx = 16.0_wp

b2 = rdx*rdx
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
      if (a(j,i) /= Zero) then
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
        if (a(i,j) /= Zero) then
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
  scl(k:l) = One
  ! .......... iterative loop for norm reduction ..........
  do
    noconv = .false.

    do i=k,l
      c = Zero
      r = Zero

      do j=k,l
        if (j == i) cycle
        c = c+abs(a(j,i))
        r = r+abs(a(i,j))
      end do
      ! .......... guard against zero c or r due to underflow ..........
      if ((c == Zero) .or. (r == Zero)) cycle
      g = r/rdx
      f = One
      s = c+r
      do while (c < g)
        f = f*rdx
        c = c*b2
      end do
      g = r*rdx
      do while (c >= g)
        f = f/rdx
        c = c/b2
      end do
      ! .......... now balance ..........
      if ((c+r)/f >= 0.95_wp*s) cycle
      g = One/f
      scl(i) = scl(i)*f
      noconv = .true.

      a(i,k:) = a(i,k:)*g

      a(1:l,i) = a(1:l,i)*f

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
  integer(kind=iwp) :: i
  scl(m) = j
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
