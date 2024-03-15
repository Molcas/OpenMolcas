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

subroutine get_pivot_idx(Diag,n,m,lu_A0,lu_A,iD_A,Scr,lScr,Thr)
!***********************************************************************
!                                                                      *
!     Author:  F. Aquilante                                            *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: Diag(*)
integer(kind=iwp), intent(in) :: n, lu_A0, lu_A, lScr
integer(kind=iwp), intent(out) :: m, iD_A(n)
real(kind=wp), intent(out) :: Scr(lScr)
real(kind=wp), intent(in) :: Thr
#include "warnings.h"
integer(kind=iwp) :: i, iAddr, iD_Col, ij, is, istart, js, k, kAddr, kCol, ks, kScr, lindep, lmax, nMem_Col
real(kind=wp) :: Acc, XMax
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: n_NegInpDiag
real(kind=wp) :: d_NegInpDiag
#endif
integer(kind=iwp), allocatable :: list(:)

#ifdef _DEBUGPRINT_
!-tbp: check diagonal for negative entries
n_NegInpDiag = 0
d_NegInpDiag = Zero
do i=1,n
  if (Diag(i) < Zero) then
    n_NegInpDiag = n_NegInpDiag+1
    if (Diag(i) < d_NegInpDiag) then
      d_NegInpDiag = Diag(i)
    end if
  end if
end do
write(u6,'(A,I10,A,I10)') 'GET_PIVOT_IDX: number of negative input diagonals:',n_NegInpDiag,' out of ',n
if (n_NegInpDiag > 0) then
  write(u6,'(A,ES12.4)') 'GET_PIVOT_IDX: most negative diagonal:          ',d_NegInpDiag
end if
#endif

Acc = min(1.0e-12_wp,thr*1.0e-2_wp)
call mma_Allocate(List,n,Label='List')
do i=1,n
  list(i) = i
end do

lmax = lScr-2*n
if (lmax < n) then
  call WarningMessage(2,'Error in Get_Pivot_idx')
  write(u6,*) ' Get_Pivot_idx: too little scratch space!! '
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
    if (abs(Diag(i)) > xMax+Acc) then
      iD_Col = i
      xMax = abs(Diag(i))
    end if
  end do
  if ((iD_Col < 0) .or. (iD_Col > n)) then
    write(u6,*) 'Get_Pivot_id: Index of Max Diag out of bounds!'
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
  if (kCol > nMem_Col) call dDaFile(lu_A,1,Scr,n,iAddr)

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
  write(u6,*) 'Get_Pivot_id: m > n is not possible!'
  call Abend()
end if
call mma_deallocate(List)

return

end subroutine get_pivot_idx
