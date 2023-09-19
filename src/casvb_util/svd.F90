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

subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
! this subroutine is a translation of the algol procedure svd,
! num. math. 14, 403-420(1970) by golub and reinsch.
! handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
!
! this subroutine determines the singular value decomposition
!      t
! a=usv  of a real m by n rectangular matrix.  householder
! bidiagonalization and a variant of the qr algorithm are used.
!
! on input
!
!    nm must be set to the row dimension of two-dimensional
!      array parameters as declared in the calling program
!      dimension statement.  note that nm must be at least
!      as large as the maximum of m and n.
!
!    m is the number of rows of a (and u).
!
!    n is the number of columns of a (and u) and the order of v.
!
!    a contains the rectangular input matrix to be decomposed.
!
!    matu should be set to .true. if the u matrix in the
!      decomposition is desired, and to .false. otherwise.
!
!    matv should be set to .true. if the v matrix in the
!      decomposition is desired, and to .false. otherwise.
!
! on output
!
!    a is unaltered (unless overwritten by u or v).
!
!    w contains the n (non-negative) singular values of a (the
!      diagonal elements of s).  they are unordered.  if an
!      error exit is made, the singular values should be correct
!      for indices ierr+1,ierr+2,...,n.
!
!    u contains the matrix u (orthogonal column vectors) of the
!      decomposition if matu has been set to .true.  otherwise
!      u is used as a temporary array.  u may coincide with a.
!      if an error exit is made, the columns of u corresponding
!      to indices of correct singular values should be correct.
!
!    v contains the matrix v (orthogonal) of the decomposition if
!      matv has been set to .true.  otherwise v is not referenced.
!      v may also coincide with a if u is not needed.  if an error
!      exit is made, the columns of v corresponding to indices of
!      correct singular values should be correct.
!
!    ierr is set to
!      zero       for normal return,
!      k          if the k-th singular value has not been
!                 determined after 30 iterations.
!
!    rv1 is a temporary storage array.
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
!
!  TESTING FOR EQUALITY OR ZERO NOW DONE USING THRESHOLD, TT, 110100.

integer i, j, k, l, m, n, ii, i1, kk, k1, ll, l1, mn, nm, its, ierr
real*8 a(nm,n), w(n), u(nm,n), v(nm,n), rv1(n)
real*8 c, f, g, h, s, x, y, z, tst1, tst2, scale, pythag, thresh
logical matu, matv, skip
save thresh
data thresh/1d-16/

! The next two initializations are to appease a compiler
f = 0.0d0
h = 0.0d0

ierr = 0
l = 0 ! dummy initialize
l1 = 0 ! dummy initialize

do i=1,m
  do j=1,n
    u(i,j) = a(i,j)
  end do
end do
! .......... householder reduction to bidiagonal form ..........
g = 0.0d0
scale = 0.0d0
x = 0.0d0

do i=1,n
  l = i+1
  rv1(i) = scale*g
  g = 0.0d0
  s = 0.0d0
  scale = 0.0d0

  if (i <= m) then
    do k=i,m
      scale = scale+abs(u(k,i))
    end do

    if (abs(scale) >= thresh) then
      do k=i,m
        u(k,i) = u(k,i)/scale
        s = s+u(k,i)**2
      end do

      f = u(i,i)
      g = -sign(sqrt(s),f)
      h = f*g-s
      u(i,i) = f-g

      if (i /= n) then
        do j=l,n
          s = 0.0d0

          do k=i,m
            s = s+u(k,i)*u(k,j)
          end do

          f = s/h

          do k=i,m
            u(k,j) = u(k,j)+f*u(k,i)
          end do
        end do
      end if

      do k=i,m
        u(k,i) = scale*u(k,i)
      end do
    end if
  end if

  w(i) = scale*g
  g = 0.0d0
  s = 0.0d0
  scale = 0.0d0

  if ((i <= m) .and. (i /= n)) then
    do k=l,n
      scale = scale+abs(u(i,k))
    end do

    if (abs(scale) >= thresh) then
      do k=l,n
        u(i,k) = u(i,k)/scale
        s = s+u(i,k)**2
      end do

      f = u(i,l)
      g = -sign(sqrt(s),f)
      h = f*g-s
      u(i,l) = f-g

      do k=l,n
        rv1(k) = u(i,k)/h
      end do

      if (i /= m) then
        do j=l,m
          s = 0.0d0

          do k=l,n
            s = s+u(j,k)*u(i,k)
          end do

          do k=l,n
            u(j,k) = u(j,k)+s*rv1(k)
          end do
        end do
      end if

      do k=l,n
        u(i,k) = scale*u(i,k)
      end do
    end if
  end if

  x = max(x,abs(w(i))+abs(rv1(i)))
end do
! .......... accumulation of right-hand transformations ..........
if (matv) then
  ! .......... for i=n step -1 until 1 do -- ..........
  do ii=1,n
    i = n+1-ii

    if (i /= n) then
      if (abs(g) >= thresh) then
        do j=l,n
          ! .......... double division avoids possible underflow ..........
          v(j,i) = (u(i,j)/u(i,l))/g
        end do

        do j=l,n
          s = 0.0d0

          do k=l,n
            s = s+u(i,k)*v(k,j)
          end do

          do k=l,n
            v(k,j) = v(k,j)+s*v(k,i)
          end do
        end do
      end if

      do j=l,n
        v(i,j) = 0.0d0
        v(j,i) = 0.0d0
      end do
    end if

    v(i,i) = 1.0d0
    g = rv1(i)
    l = i
  end do
end if
! .......... accumulation of left-hand transformations ..........
if (matu) then
  ! ..........for i=min(m,n) step -1 until 1 do -- ..........
  mn = n
  if (m < n) mn = m

  do ii=1,mn
    i = mn+1-ii
    l = i+1
    g = w(i)

    if (i /= n) then
      do j=l,n
        u(i,j) = 0.0d0
      end do
    end if

    if (abs(g) < thresh) then
      do j=i,m
        u(j,i) = 0.0d0
      end do
    else
      if (i /= mn) then
        do j=l,n
          s = 0.0d0

          do k=l,m
            s = s+u(k,i)*u(k,j)
          end do
          ! .......... double division avoids possible underflow ..........
          f = (s/u(i,i))/g

          do k=i,m
            u(k,j) = u(k,j)+f*u(k,i)
          end do
        end do
      end if

      do j=i,m
        u(j,i) = u(j,i)/g
      end do
    end if

    u(i,i) = u(i,i)+1.0d0
  end do
end if
! .......... diagonalization of the bidiagonal form ..........
tst1 = x
! .......... for k=n step -1 until 1 do -- ..........
do kk=1,n
  k1 = n-kk
  k = k1+1
  its = 0

  do
    ! .......... test for splitting. for l=k step -1 until 1 do -- ..........
    skip = .false.
    do ll=1,k
      l1 = k-ll
      l = l1+1
      tst2 = tst1+abs(rv1(l))
      if (abs(tst2-tst1) < thresh) then
        skip = .true.
        exit
      end if
      ! .......... rv1(1) is always zero, so there is no exit through the bottom of the loop ..........
      tst2 = tst1+abs(w(l1))
      if (abs(tst2-tst1) < thresh) exit
    end do
    if (.not. skip) then
      ! .......... cancellation of rv1(l) if l greater than 1 ..........
      c = 0.0d0
      s = 1.0d0

      do i=l,k
        f = s*rv1(i)
        rv1(i) = c*rv1(i)
        tst2 = tst1+abs(f)
        if (abs(tst2-tst1) < thresh) exit
        g = w(i)
        h = pythag(f,g)
        w(i) = h
        c = g/h
        s = -f/h

        if (matu) then
          do j=1,m
            y = u(j,l1)
            z = u(j,i)
            u(j,l1) = y*c+z*s
            u(j,i) = -y*s+z*c
          end do
        end if

      end do
    end if
    ! .......... test for convergence ..........
    z = w(k)
    if (l == k) then
      ! .......... convergence ..........
      if (z >= 0.0d0) exit
      ! .......... w(k) is made non-negative ..........
      w(k) = -z

      if (matv) then
        do j=1,n
          v(j,k) = -v(j,k)
        end do
      end if
      exit
    else
      ! .......... shift from bottom 2 by 2 minor ..........
      if (its == 30) then
        ! .......... set error -- no convergence to a singular value after 30 iterations ..........
        ierr = k
        return
      end if
      its = its+1
      x = w(l)
      y = w(k1)
      g = rv1(k1)
      h = rv1(k)
      f = 0.5d0*(((g+z)/h)*((g-z)/y)+y/h-h/y)
      g = pythag(f,1.0d0)
      f = x-(z/x)*z+(h/x)*(y/(f+sign(g,f))-h)
      ! .......... next qr transformation ..........
      c = 1.0d0
      s = 1.0d0

      do i1=l,k1
        i = i1+1
        g = rv1(i)
        y = w(i)
        h = s*g
        g = c*g
        z = pythag(f,h)
        rv1(i1) = z
        c = f/z
        s = h/z
        f = x*c+g*s
        g = -x*s+g*c
        h = y*s
        y = y*c

        if (matv) then
          do j=1,n
            x = v(j,i1)
            z = v(j,i)
            v(j,i1) = x*c+z*s
            v(j,i) = -x*s+z*c
          end do
        end if

        z = pythag(f,h)
        w(i1) = z
        ! .......... rotation can be arbitrary if z is zero ..........
        if (abs(z) >= thresh) then
          c = f/z
          s = h/z
        end if
        f = c*g+s*y
        x = -s*g+c*y

        if (matu) then
          do j=1,m
            y = u(j,i1)
            z = u(j,i)
            u(j,i1) = y*c+z*s
            u(j,i) = -y*s+z*c
          end do
        end if

      end do

      rv1(l) = 0.0d0
      rv1(k) = f
      w(k) = x
    end if

  end do

end do

return

end subroutine svd
