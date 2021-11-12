!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1983, Per Ake Malmqvist                                *
!               1994,1995, Niclas Forsberg                             *
!***********************************************************************

!  Contains:
!    Cholesky    : Performs a Cholesky factorization of a positive
!                 definite matrix A.
!    Dool_MULA   : Solves the system of linear equations A*X = B.
!    SolveSecEq  : Solves the seqular equation A*C = S*C*D.
!    Polfit      : Fit a polynomial of requested dimension to the input
!                  data.
!    calcNorm    : Calculates norm of a vector.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

subroutine Dool_MULA(A,LA1,LA2,B,LB1,LB2,det)
!  Purpose:
!    Solve A*X = B
!
!  Input:
!    A      : Real*8 two dimensional array
!    B      : Real*8 two dimensional array
!
!  Output:
!    B      : Real*8 two dimensional array - contains
!             the solution X.
!  Written by:
!    Per-AAke Malmquist
!    Dept. of Theoretical Chemistry, Lund University, 1983.
!
!  Modified by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One

implicit real*8(a-h,o-z)
real*8 A(LA1,LA2)
real*8 B(LB1,LB2)
integer, allocatable :: iPiv(:), jPiv(:)
real*8, allocatable :: Buf(:)

! Initialize.
ip = -9999999
jp = -9999999
n = LA2
m = LB2
call mma_allocate(Buf,n,label='Buf')
call mma_allocate(iPiv,n,label='iPiv')
call mma_allocate(jPiv,n,label='jPiv')

do i=1,n
  iPiv(i) = i
  jPiv(i) = i
end do
det = One
do i=1,n
  ! Now find better pivot element.
  Amax = -One
  do k=i,n
    do l=i,n
      Am = abs(A(iPiv(k),jPiv(l)))
      if (Amax <= Am) then
        Amax = Am
        ip = k
        jp = l
      end if
    end do
  end do
  if (ip /= i) then
    det = -det
    iTemp = iPiv(i)
    iPiv(i) = iPiv(ip)
    iPiv(ip) = iTemp
  end if
  if (jp /= i) then
    det = -det
    jTemp = jPiv(i)
    jPiv(i) = jPiv(jp)
    jPiv(jp) = jTemp
  end if
  ip = iPiv(i)
  jp = jPiv(i)
  diag = A(ip,jp)
  !Buf(i) = diag
  Buf(i) = diag
  det = det*diag
  do k=i+1,n
    kp = iPiv(k)
    c = A(kp,jp)/diag
    A(kp,jp) = c
    do l=i+1,n
      lp = jPiv(l)
      A(kp,lp) = A(kp,lp)-c*A(ip,lp)
    end do
  end do
end do

! First resubstitution step.
do j=1,m
  do i=2,n
    ip = iPiv(i)
    sum = B(ip,j)
    do k=1,i-1
      sum = sum-A(ip,jPiv(k))*B(iPiv(k),j)
    end do
    B(ip,j) = sum
  end do
end do
!!
!!---- Second resubstitution step.
do j=1,m
  do i=n,1,-1
    ip = iPiv(i)
    sum = B(ip,j)
    do k=i+1,n
      sum = sum-A(ip,jPiv(k))*B(iPiv(k),j)
    end do
    B(ip,j) = sum/Buf(i)
  end do
end do

! Reorganization part.
do j=1,m
  do i=1,n
    Buf(i) = B(iPiv(i),j)
  end do
  do i=1,n
    B(jPiv(i),j) = Buf(i)
  end do
end do

call mma_deallocate(Buf)
call mma_deallocate(iPiv)
call mma_deallocate(jPiv)

end subroutine Dool_MULA
!####
subroutine SolveSecEq(A,n,C,S,D)
!  Purpose:
!    Solve the secular equation SAC = CD.
!
!  Input:
!    A        : Real*8 two dimensional array
!    S        : Real*8 two dimensional array
!
!  Output:
!    C        : Real*8 two dimensional array
!    D        : Real*8 array
!
!  Calls:
!    Jacob
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

implicit real*8(a-h,o-z)
integer n
real*8 A(n,n), S(n,n), C(n,n), D(n)
real*8, allocatable :: Asymm(:,:), Scr(:), T(:,:), Temp(:,:)

! Initialize.
!D write(u6,*) 'SolveSecEq test prints.'
!D write(u6,*) 'Matrix A:'
!D do i=1,n
!D   write(u6,'(1x,5F16.8)') (A(i,j),j=1,n)
!D end do

! get memory for temporary matrices.
call mma_allocate(Scr,n*(n+1)/2,label='Scr')
call mma_allocate(T,n,n,label='T')
call mma_allocate(Temp,n,n,label='Temp')
call mma_allocate(Asymm,n,n,label='Asymm')

! Transform S to lower packed storage in Scratch.
k = 1
do i=1,n
  do j=1,i
    Scr(k) = S(i,j)
    k = k+1
  end do
end do

! Turn T into a unit matrix.
T(:,:) = Zero
do i=1,n
  T(i,i) = One
end do

! Diagonalize Scratch and scale each column of T with the square
! root of the corresponding eigenvalue.
call Jacob(Scr,T,n,n)
do j=1,n
  jj = j*(j+1)/2
  T(:,j) = T(:,j)*sqrt(Scr(jj))
end do

! Make A symmetric and transform it to lower packed storage in
! Scratch.
call DGEMM_('N','N',n,n,n,One,A,n,T,n,Zero,Temp,n)
call DGEMM_('T','N',n,n,n,One,T,n,Temp,n,Zero,Asymm,n)
k = 1
do i=1,n
  do j=1,i
    Scr(k) = Asymm(i,j)
    k = k+1
  end do
end do

! Diagonalize Scratch.
call Jacob(Scr,T,n,n)
call JacOrd(Scr,T,n,n)

! Store the eigenvalues in array D.
do i=1,n
  ii = i*(i+1)/2
  D(i) = Scr(ii)
end do
C(:,:) = T

! Free memory of temporary matrices.
call mma_deallocate(Scr)
call mma_deallocate(T)
call mma_deallocate(Temp)
call mma_deallocate(Asymm)

end subroutine SolveSecEq
!####
subroutine PolFit(ipow,nvar,var,yin,ndata,coef,nterm,stand_dev,max_err,diff_vec,use_weight)
!  Purpose:
!    Fit a polynomial to yin.
!
!  Written by:
!    P-AA Malmquist
!
!  Modified by:
!    Niclas Forsberg

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
parameter(mxdeg=6)
integer ipow(nvar,nterm)
real*8 var(ndata,nvar)
real*8 yin(ndata)
real*8 coef(nterm,1)
real*8 stand_dev, max_err
real*8 yfit(ndata)
real*8 vpow(0:mxdeg,nvar)
real*8 equmat(nterm,nterm)
real*8 rhs(nterm), term(nterm)
real*8 diff_vec(ndata)
logical use_weight
real*8, allocatable :: weight(:)

! Initialize.
NrOfVar = nvar
nPolyTerm = nterm
myCoef2 = 1
rhs(:) = Zero
equmat(:,:) = Zero

! Set up weight vector.
call mma_allocate(weight,ndata,label='weight')
if (use_weight) then
  e_min = yin(1)
  e_max = yin(1)
  do idata=2,ndata
    if (yin(idata) < e_min) e_min = yin(idata)
    if (yin(idata) > e_max) e_max = yin(idata)
  end do
  e_range = e_max-e_min
  do idata=1,ndata
    weight(idata) = One/(One+1.0e3_wp*((yin(idata)-e_min)/e_range))
  end do
else
  !vv weight(:) = One
  weight(:) = 0.1_wp
end if

! Accumulate equation matrix and right-hand-side.
do idata=1,ndata
  ! Calculate powers of individual variable values.
  do ivar=1,NrOfVar
    pow = One
    vpow(0,ivar) = One
    do i=1,mxdeg
      pow = pow*var(idata,ivar)
      vpow(i,ivar) = pow
    end do
  end do
  ! Calculate value of each polynomial term at this point.
  do iterm=1,nPolyTerm
    ip = ipow(iterm,1)
    t = vpow(ip,1)
    do ivar=2,NrOfVar
      ip = ipow(iterm,ivar)
      t = t*vpow(ip,ivar)
    end do
    term(iterm) = t
  end do
  ! Accumulate equmat and rhs.
  do iterm=1,nPolyTerm
    rhs(iterm) = rhs(iterm)+yin(idata)*term(iterm)*weight(idata)
    do jterm=1,iterm
      equmat(iterm,jterm) = equmat(iterm,jterm)+term(iterm)*term(jterm)*weight(idata)
    end do
  end do
end do

! Set upper triangle of equmat by symmetry.
do iterm=1,nPolyTerm-1
  do jterm=iterm+1,nPolyTerm
    equmat(iterm,jterm) = equmat(jterm,iterm)
  end do
end do

! Solve the resulting equation system.
do i=1,nterm
  coef(i,1) = rhs(i)
end do
call Dool_MULA(equmat,nPolyTerm,nPolyTerm,coef,nPolyTerm,MyCoef2,det)
if (abs(det) == Zero) write(u6,*) 'WARNING!! Determinant=0 in PolFit'

! Calculate fitted result.
do idata=1,ndata
  ! Calculate powers of individual variable values.
  do ivar=1,NrOfVar
    pow = One
    vpow(0,ivar) = One
    do i=1,mxdeg
      pow = pow*var(idata,ivar)
      vpow(i,ivar) = pow
    end do
  end do
  ! Calculate value of each polynomial term. Add to yfit.
  pol = Zero
  do iterm=1,nPolyTerm
    ip = ipow(iterm,1)
    t = vpow(ip,1)
    do ivar=2,NrOfVar
      ip = ipow(iterm,ivar)
      t = t*vpow(ip,ivar)
    end do
    pol = pol+coef(iterm,1)*t
  end do
  yfit(idata) = pol
end do

! Calculate standard deviation and maximum error.
sum = Zero
do idata=1,ndata
  diff = abs(yin(idata)-yfit(idata))
  diff_vec(idata) = diff
  if (idata == 1) then
    max_err = diff
  else
    if (diff > max_err) max_err = diff
  end if
  sum = sum+diff**2
end do
stand_dev = sqrt(sum/ndata)

call mma_deallocate(weight)

end subroutine PolFit
!####
subroutine factor(exponent,nder,rfactor)
!  Purpose:
!    Calculate coefficient for n'th derivative.

use Definitions, only: wp

implicit none
integer exponent, nder
integer i, isum, ntot
real*8 rfactor

isum = 1
ntot = nder+exponent
if (nder > 0) then
  do i=ntot,ntot-nder+1,-1
    isum = isum*i
  end do
else if (nder < 0) then
  isum = 0
end if
rfactor = real(isum,kind=wp)

end subroutine factor
!####
subroutine Cholesky(A,L,nd)
!  Purpose:
!    Perform a Cholesky factorization of the matrix A:
!
!                       A = L~ * L
!
!    where L is a lower triangular matrix and L~ is L transposed.
!
!  Input:
!    A      : real*8 two dimensional array - positive
!             definite matrix.
!
!  Output:
!    L      : real*8 two dimensional array - lower
!             triangular matrix.
!
!  Calls:

use Constants, only: Zero, One
use Definitions, only: u6

use stdalloc, only: mma_allocate, mma_deallocate

!implicit none
!VV: all calls use nxn
real*8 A(nd,nd), L(nd,nd), dd
integer n, j
integer iRow, jRow
real*8, allocatable :: D(:)

! Initialize.
n = nd
call mma_allocate(D,n,label='D')
L(:,:) = A

! Take care of the n-1 last rows.
if (n > 1) then
  do iRow=n,2,-1
    D(iRow) = L(iRow,iRow)
    L(iRow,:) = L(iRow,:)/D(iRow)
    do jRow=iRow-1,1,-1
      do j=1,n
        L(jRow,j) = L(jRow,j)-L(jRow,iRow)*L(iRow,j)
      end do
    end do
  end do
end if

! Take care of the first row.
D(1) = L(1,1)
L(1,1) = One

! Multiply Llow with the square root of d.
do iRow=1,n
  if (D(iRow) < Zero) then
    write(u6,*) 'Error in Cholesky!!! Matrix not positive definite.'
    call Abend()
  end if
  dd = sqrt(D(iRow))
  L(iRow,:) = L(iRow,:)*dd
end do
call mma_deallocate(D)

!Llow(:,:) = L

end subroutine Cholesky
