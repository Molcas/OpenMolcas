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

subroutine Pivot_mat(n,m,lu_A0,lu_A,iD_A,Scr,lScr)
!***********************************************************************
!                                                                      *
!     Author:  F. Aquilante                                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, m, lu_A0, lu_A, iD_A(n), lScr
real(kind=wp), intent(out) :: Scr(lScr)
#include "warnings.h"
integer(kind=iwp) :: i, iAddr, iCol, ij, iScr, jCol, kAddr, kCol, lmax, mNeed, nMem_Col

lmax = lScr-n
if (lmax < n) then
  call WarningMessage(2,'Error in Pivot_mat')
  write(u6,*) ' Pivot_mat: too little scratch space !!'
  call Quit(_RC_CHO_LOG_)
end if

nMem_Col = m
mNeed = nTri_Elem(nMem_Col)
do while (mNeed > lmax)
  mNeed = mNeed-nMem_Col
  nMem_Col = nMem_Col-1
end do

iScr = n
do kCol=1,nMem_Col
  jCol = iD_A(kCol)
  iAddr = n*(jCol-1)
  call dDaFile(lu_A0,2,Scr,n,iAddr)
  do i=1,kCol
    iCol = iD_A(i)
    iScr = iScr+1
    Scr(iScr) = Scr(iCol)
  end do
end do
kAddr = 0
ij = nTri_Elem(nMem_Col)
call dDaFile(lu_A,1,Scr(n+1),ij,kAddr)

do kCol=nMem_Col+1,m
  jCol = iD_A(kCol)
  iAddr = n*(jCol-1)
  call dDaFile(lu_A0,2,Scr,n,iAddr)
  do i=1,kCol
    iCol = iD_A(i)
    iScr = n+i
    Scr(iScr) = Scr(iCol)
  end do
  call dDaFile(lu_A,1,Scr(n+1),kCol,kAddr)
end do

return

end subroutine Pivot_mat
