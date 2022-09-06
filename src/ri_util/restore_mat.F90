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

subroutine Restore_mat(n,m,lu_A0,lu_A,iD_A,Scr,lScr,Add0s)
!***********************************************************************
!
!     Author:  F. Aquilante
!
!***********************************************************************

implicit real*8(a-h,o-z)
integer n, m, lu_A0, lu_A, iD_A(n), lScr
real*8 Scr(lScr)
logical Add0s
#include "warnings.h"

lmax = lScr-n
if (lmax < n) then
  call WarningMessage(2,'Error in Restore_mat')
  write(6,*) ' Restore_mat: too little scratch space!! '
  call Quit(_RC_CHO_LOG_)
end if

nMem_Col = m
mNeed = nMem_Col*(nMem_Col+1)/2
do while (mNeed > lmax)
  mNeed = mNeed-nMem_Col
  nMem_Col = nMem_Col-1
end do

kAddr = 0
ij = nMem_Col*(nMem_Col+1)/2
call dDaFile(lu_A0,2,Scr(1),ij,kAddr)

iOff = 0
do kCol=1,nMem_Col
  do i=1,kCol
    iCol = iD_A(i)
    iScr = ij+iCol
    jCol = iOff+i
    Scr(iScr) = Scr(jCol)
  end do
  do i=kCol+1,n
    iCol = iD_A(i)
    iScr = ij+iCol
    Scr(iScr) = 0.0d0
  end do
  iAddr = n*(kCol-1)
  call dDaFile(lu_A,1,Scr(ij+1),n,iAddr)
  !call RecPrt('QVec',' ',Scr(ij+1),1,n)
  iOff = iOff+kCol
end do

do kCol=nMem_Col+1,m
  call dDaFile(lu_A0,2,Scr(1),kCol,kAddr)
  do i=1,kCol
    iCol = iD_A(i)
    iScr = n+iCol
    Scr(iScr) = Scr(i)
  end do
  do i=kCol+1,n
    iCol = iD_A(i)
    iScr = n+iCol
    Scr(iScr) = 0.0d0
  end do
  iAddr = n*(kCol-1)
  call dDaFile(lu_A,1,Scr(n+1),n,iAddr)
  !call RecPrt('QVec',' ',Scr(n+1),1,n)
end do

if (Add0s) then
  do kCol=m+1,n   ! linearly dependent cols
    iAddr = n*(kCol-1)
    call FZero(Scr,n)
    call dDaFile(lu_A,1,Scr(1),n,iAddr)
  end do
end if

return

end subroutine Restore_mat
