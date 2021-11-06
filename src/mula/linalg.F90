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

implicit real*8(a-h,o-z)
real*8 A(LA1,LA2)
real*8 B(LB1,LB2)
#include "WrkSpc.fh"

! Initialize.
ip = -9999999
jp = -9999999
n = LA2
m = LB2
call GetMem('Buf','Allo','Real',ipBuf,n)
call GetMem('iPiv','Allo','Inte',ipiPiv,n)
call GetMem('jPiv','Allo','Inte',ipjPiv,n)

do i=1,n
  iWork(ipiPiv+i-1) = i
  iWork(ipjPiv+i-1) = i
end do
det = 1.0d0
do i=1,n
  ! Now find better pivot element.
  Amax = -1.0d0
  do k=i,n
    do l=i,n
      Am = abs(A(iWork(ipiPiv+k-1),iWork(ipjPiv+l-1)))
      if (Amax <= Am) then
        Amax = Am
        ip = k
        jp = l
      end if
    end do
  end do
  if (ip /= i) then
    det = -det
    iTemp = iWork(ipiPiv+i-1)
    iWork(ipiPiv+i-1) = iWork(ipiPiv+ip-1)
    iWork(ipiPiv+ip-1) = iTemp
  end if
  if (jp /= i) then
    det = -det
    jTemp = iWork(ipjPiv+i-1)
    iWork(ipjPiv+i-1) = iWork(ipjPiv+jp-1)
    iWork(ipjPiv+jp-1) = jTemp
  end if
  ip = iWork(ipiPiv+i-1)
  jp = iWork(ipjPiv+i-1)
  diag = A(ip,jp)
  !Buf(i) = diag
  Work(ipBuf+i-1) = diag
  det = det*diag
  do k=i+1,n
    kp = iWork(ipiPiv+k-1)
    c = A(kp,jp)/diag
    A(kp,jp) = c
    do l=i+1,n
      lp = iWork(ipjPiv+l-1)
      A(kp,lp) = A(kp,lp)-c*A(ip,lp)
    end do
  end do
end do

! First resubstitution step.
do j=1,m
  do i=2,n
    ip = iWork(ipiPiv+i-1)
    sum = B(ip,j)
    do k=1,i-1
      sum = sum-A(ip,iWork(ipjPiv+k-1))*B(iWork(ipiPiv+k-1),j)
    end do
    B(ip,j) = sum
  end do
end do
!!
!!---- Second resubstitution step.
do j=1,m
  do i=n,1,-1
    ip = iWork(ipiPiv+i-1)
    sum = B(ip,j)
    do k=i+1,n
      sum = sum-A(ip,iWork(ipjPiv+k-1))*B(iWork(ipiPiv+k-1),j)
    end do
    B(ip,j) = sum/Work(ipBuf+i-1)
  end do
end do

! Reorganization part.
do j=1,m
  do i=1,n
    Work(ipBuf+i-1) = B(iWork(ipiPiv+i-1),j)
  end do
  do i=1,n
    B(iWork(ipjPiv+i-1),j) = Work(ipBuf+i-1)
  end do
end do

call GetMem('Buf','Free','Real',ipBuf,n)
call GetMem('iPiv','Free','Real',ipiPiv,n)
call GetMem('jPiv','Free','Real',ipjPiv,n)

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

implicit real*8(a-h,o-z)
integer n
real*8 A(n,n), S(n,n), C(n,n), D(n)
#include "WrkSpc.fh"

! Initialize.
nSqr = n**2
nSqrTri = n*(n+1)/2
!D write(6,*) 'SolveSecEq test prints.'
!D write(6,*) 'Matrix A:'
!D do i=1,n
!D   write(6,'(1x,5F16.8)') (A(i,j),j=1,n)
!D end do

! get memory for temporary matrices.
call GetMem('Scr','Allo','Real',ipScr,nSqrTri)
call GetMem('T','Allo','Real',ipT,nSqr)
call GetMem('Temp','Allo','Real',ipTemp,nSqr)
call GetMem('Asymm','Allo','Real',ipAsymm,nSqr)

! Transform S to lower packed storage in Scratch.
k = 1
do i=1,n
  do j=1,i
    Work(ipScr+k-1) = S(i,j)
    k = k+1
  end do
end do

! Turn T into a unit matrix.
!vv T = 0.0d0
call dcopy_(nSqr,[0.0d0],0,Work(ipT),1)
do i=1,n
  Work(ipT+i+n*(i-1)-1) = 1.0d0
end do

! Diagonalize Scratch and scale each column of T with the square
! root of the corresponding eigenvalue.
call Jacob(Work(ipScr),Work(ipT),n,n)
do j=1,n
  jj = j*(j+1)/2
  Scale = sqrt(Work(ipScr+jj-1))
  do i=1,n
    Work(ipT+i+n*(j-1)-1) = Work(ipT+i+n*(j-1)-1)*Scale
  end do
end do

! Make A symmetric and transform it to lower packed storage in
! Scratch.
call DGEMM_('N','N',n,n,n,1.0d0,A,n,Work(ipT),n,0.0d0,Work(ipTemp),n)
call DGEMM_('T','N',n,n,n,1.0d0,Work(ipT),n,Work(ipTemp),n,0.0d0,Work(ipAsymm),n)
k = 1
do i=1,n
  do j=1,i
    Work(ipScr+k-1) = Work(ipAsymm+i+n*(j-1)-1)
    k = k+1
  end do
end do

! Diagonalize Scratch.
call Jacob(Work(ipScr),Work(ipT),n,n)
call JacOrd(Work(ipScr),Work(ipT),n,n)

! Store the eigenvalues in array D.
do i=1,n
  ii = i*(i+1)/2
  D(i) = Work(ipScr+ii-1)
end do
!vv C = T
call dcopy_(nSqr,Work(ipT),1,C,1)

! Free memory of temporary matrices.
call GetMem('Scr','Free','Real',ipScr,nSqrTri)
call GetMem('T','Free','Real',ipT,nSqr)
call GetMem('Temp','Free','Real',ipTemp,nSqr)
call GetMem('Asymm','Free','Real',ipAsymm,nSqr)

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
#include "WrkSpc.fh"

! Initialize.
NrOfVar = nvar
nPolyTerm = nterm
myCoef2 = 1
!vv rhs = 0.0d0
call dcopy_(nterm,[0.0d0],0,rhs,1)
!vv equmat = 0.0d0
call dcopy_(nterm*nterm,[0.0d0],0,equmat,1)

! Set up weight vector.
call GetMem('weight','Allo','Real',ipweight,ndata)
if (use_weight) then
  e_min = yin(1)
  e_max = yin(1)
  do idata=2,ndata
    if (yin(idata) < e_min) e_min = yin(idata)
    if (yin(idata) > e_max) e_max = yin(idata)
  end do
  e_range = e_max-e_min
  do idata=1,ndata
    Work(ipweight+idata-1) = 1.0d0/(1.0d0+1000.0d0*((yin(idata)-e_min)/e_range))
  end do
else
  !vv weight = 1.0d0
  call dcopy_(ndata,[0.1d0],0,Work(ipweight),1)
end if

! Accumulate equation matrix and right-hand-side.
do idata=1,ndata
  ! Calculate powers of individual variable values.
  do ivar=1,NrOfVar
    pow = 1.0d0
    vpow(0,ivar) = 1.0d0
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
    rhs(iterm) = rhs(iterm)+yin(idata)*term(iterm)*Work(ipweight+idata-1)
    do jterm=1,iterm
      equmat(iterm,jterm) = equmat(iterm,jterm)+term(iterm)*term(jterm)*Work(ipweight+idata-1)
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
if (abs(det) == 0.0) write(6,*) 'WARNING!! Determinant=0 in PolFit'

! Calculate fitted result.
do idata=1,ndata
  ! Calculate powers of individual variable values.
  do ivar=1,NrOfVar
    pow = 1.0d0
    vpow(0,ivar) = 1.0d0
    do i=1,mxdeg
      pow = pow*var(idata,ivar)
      vpow(i,ivar) = pow
    end do
  end do
  ! Calculate value of each polynomial term. Add to yfit.
  pol = 0.0d0
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
sum = 0.0d0
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

call GetMem('weight','Free','Real',ipweight,ndata)

end subroutine PolFit
!####
subroutine factor(exponent,nder,rfactor)
!  Purpose:
!    Calculate coefficient for n'th derivative.

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
rfactor = dble(isum)

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

!implicit none
!VV: all calls use nxn
real*8 A(nd,nd), L(nd,nd), dd
integer n, j
integer iRow, jRow
#include "WrkSpc.fh"

! Initialize.
n = nd
n2 = n*n
call GetMem('D','Allo','Real',ipD,n)
call dcopy_(n2,A,1,L,1)
!vv L = A

! Take care of the n-1 last rows.
if (n > 1) then
  do iRow=n,2,-1
    work(ipd+iRow-1) = L(iRow,iRow)
    do j=1,n
      L(iRow,j) = L(iRow,j)/work(ipd+iRow-1)
    end do
    do jRow=iRow-1,1,-1
      do j=1,n
        L(jRow,j) = L(jRow,j)-L(jRow,iRow)*L(iRow,j)
      end do
    end do
  end do
end if

! Take care of the first row.
work(ipd) = L(1,1)
L(1,1) = 1.0d0

! Multiply Llow with the square root of d.
do iRow=1,n
  if (work(ipd+iRow-1) < 0.0) then
    write(6,*) 'Error in Cholesky!!! Matrix not positive definite.'
    call Abend()
  end if
  dd = sqrt(work(ipd+iRow-1))
  do j=1,n
    L(iRow,j) = L(iRow,j)*dd
  end do
end do
call GetMem('D','Free','Real',ipD,n)

!call dcopy_(n2,L,1,Llow,1)
!Llow = L

end subroutine Cholesky
