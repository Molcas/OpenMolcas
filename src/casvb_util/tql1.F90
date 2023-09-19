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
!********************************************************
!** Public-domain library routines used by casvb only. **
!********************************************************
!**********************
!** EISPACK ROUTINES **
!**********************

subroutine tql1(n,d,e,ierr)
! this subroutine is a translation of the algol procedure tql1,
! num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
! wilkinson.
! handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
! this subroutine finds the eigenvalues of a symmetric
! tridiagonal matrix by the ql method.
!
! on input
!
!    n is the order of the matrix.
!
!    d contains the diagonal elements of the input matrix.
!
!    e contains the subdiagonal elements of the input matrix
!      in its last n-1 positions.  e(1) is arbitrary.
!
!  on output
!
!    d contains the eigenvalues in ascending order.  if an
!      error exit is made, the eigenvalues are correct and
!      ordered for indices 1,2,...ierr-1, but may not be
!      the smallest eigenvalues.
!
!    e has been destroyed.
!
!    ierr is set to
!      zero       for normal return,
!      j          if the j-th eigenvalue has not been
!                 determined after 30 iterations.
!
! calls pythag for  sqrt(a*a + b*b) .
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

integer i, j, l, m, n, ii, l1, l2, mml, ierr
real*8 d(n), e(n)
real*8 c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2, pythag
logical skip

ierr = 0
c3 = 0.0d0 ! dummy initialize
s2 = 0.0d0 ! dummy initialize
if (n == 1) return

do i=2,n
  e(i-1) = e(i)
end do

f = 0.0d0
tst1 = 0.0d0
e(n) = 0.0d0

do l=1,n
  j = 0
  h = abs(d(l))+abs(e(l))
  if (tst1 < h) tst1 = h
  ! .......... look for small sub-diagonal element ..........
  do m=l,n
    tst2 = tst1+abs(e(m))
    if (tst2 == tst1) exit
    ! .......... e(n) is always zero, so there is no exit through the bottom of the loop ..........
  end do

  if (m /= l) then
    do
      if (j == 30) then
        ! .......... set error -- no convergence to an eigenvalue after 30 iterations ..........
        ierr = l
        return
      end if
      j = j+1
      ! .......... form shift ..........
      l1 = l+1
      l2 = l1+1
      g = d(l)
      p = (d(l1)-g)/(2.0d0*e(l))
      r = pythag(p,1.0d0)
      d(l) = e(l)/(p+sign(r,p))
      d(l1) = e(l)*(p+sign(r,p))
      dl1 = d(l1)
      h = g-d(l)

      do i=l2,n
        d(i) = d(i)-h
      end do

      f = f+h
      ! .......... ql transformation ..........
      p = d(m)
      c = 1.0d0
      c2 = c
      el1 = e(l1)
      s = 0.0d0
      mml = m-l
      ! .......... for i=m-1 step -1 until l do -- ..........
      do ii=1,mml
        c3 = c2
        c2 = c
        s2 = s
        i = m-ii
        g = c*e(i)
        h = c*p
        r = pythag(p,e(i))
        e(i+1) = s*r
        s = e(i)/r
        c = p/r
        p = c*d(i)-s*g
        d(i+1) = h+s*(c*g+s*d(i))
      end do

      p = -s*s2*c3*el1*e(l)/dl1
      e(l) = s*p
      d(l) = c*p
      tst2 = tst1+abs(e(l))
      if (tst2 <= tst1) exit
    end do
  end if
  p = d(l)+f
  ! .......... order eigenvalues ..........
  skip = .false.
  if (l /= 1) then
    ! .......... for i=l step -1 until 2 do -- ..........
    do ii=2,l
      i = l+2-ii
      if (p >= d(i-1)) then
        skip = .true.
        exit
      end if
      d(i) = d(i-1)
    end do
  end if

  if (.not. skip) i = 1
  d(i) = p
end do

return

end subroutine tql1
