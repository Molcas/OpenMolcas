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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine get_pivot_idx_w(Diag,Wg,n,m,lu_A0,lu_A,iD_A,Scr,lScr,Thr)
!***********************************************************************
!                                                                      *
!     Author:  F. Aquilante                                            *
!                                                                      *
!       Note:  this routine differs from Get_pivot_idx because here    *
!              the pivoting/convergence is decided based on weighted   *
!              diagonals                                               *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(in) :: Wg(*), Thr
integer(kind=iwp), intent(in) :: n, lu_A0, lu_A, lScr
integer(kind=iwp), intent(out) :: m, iD_A(n)
real(kind=wp), intent(out) :: Scr(lScr)
#include "warnings.h"
integer(kind=iwp) :: i, iAddr, iD_Col, ij, is, istart, js, k, kAddr, kCol, ks, kScr, lindep, lmax, nMem_Col
real(kind=wp) :: Acc, XMax
integer(kind=iwp), allocatable :: List(:)

Acc = min(1.0e-12_wp,thr*1.0e-2_wp)
call mma_allocate(List,n,Label='List')
do i=1,n
  List(i) = i
end do

lmax = lScr-2*n
if (lmax < n) then
  call WarningMessage(2,'Error in Get_Pivot_idx_w')
  write(u6,*) ' Get_Pivot_idx_w: too little scratch space!! '
  call Quit(_RC_CHO_LOG_)
end if

nMem_Col = min(lmax/n,n)

kAddr = 0
is = 1+n
ij = n*nMem_Col
ks = is+ij
kScr = lScr-n-ij

m = 0
do kCol=1,n

  iD_Col = 0
  XMax = Zero
  do i=1,n
    if (abs(Diag(i)*Wg(i)) > xMax+Acc) then
      iD_Col = i
      xMax = abs(Diag(i))
    end if
  end do
  if ((iD_Col < 0) .or. (iD_Col > n)) then
    write(u6,*) 'Get_Pivot_idx_w: Index of MaxDiag out of bounds!'
    write(u6,*) 'iD_Col = ',iD_Col
    call Abend()
  else if (iD_Col == 0) then
    exit
  end if
  iD_A(kCol) = iD_Col ! set the mapping

  js = n*(kCol-1)+is ! overlay A and Z
  if (kCol > nMem_Col) js = 1

  kAddr = n*(iD_Col-1)
  call dDaFile(lu_A0,2,Scr(js),n,kAddr)

  call CHO_FACTOR(Diag,Scr(js),iD_A,kCol,n,Scr(is),nMem_Col,lu_A,Scr(ks),kScr,thr,lindep)

  if (lindep /= 0) exit

  list(iD_Col) = 0
  m = m+1

  iAddr = n*(kCol-1)
  if (kCol > nMem_Col) call dDaFile(lu_A,1,Scr(1),n,iAddr)

end do

iAddr = 0
call dDaFile(lu_A,1,Scr(is),ij,iAddr)

if (m < n) then
  istart = 1
  do k=m+1,n
    do i=istart,n
      if (list(i) /= 0) then
        iD_A(k) = i
        istart = i+1
        exit
      end if
    end do
  end do
else if (m > n) then
  write(u6,*) 'Get_Pivot_idx_w: m > n is not possible!'
  call Abend()
end if
call mma_deallocate(List)

return

end subroutine get_pivot_idx_w
