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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

function MemSO1(lOper,iCmp,jCmp,iShell,jShell,iAO,jAO)
!***********************************************************************
!  Object: to compile the number of SO block which will be generated   *
!          by the current shell doublet.                               *
!          "lOper" is the irreducible representation of which the      *
!          operator of the one-electron matrix element belongs.        *
! Called from: OneEl                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             February '90                                             *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden.                                         *
!             Modified to general non-symmetric operators January '91  *
!***********************************************************************

use SOAO_Info, only: iAOtSO
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: MemSO1
integer(kind=iwp), intent(in) :: lOper, iCmp, jCmp, iShell, jShell, iAO, jAO
integer(kind=iwp) :: j1, i1, j2, j12, jCmpMx, i2

MemSO1 = 0
do j1=0,nIrrep-1
  do i1=1,iCmp
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    do j2=0,nIrrep-1
      j12 = ieor(j1,j2)
      if (.not. btest(lOper,j12)) cycle
      jCmpMx = jCmp
      if ((iShell == jShell) .and. (j1 == j2)) jCmpMx = i1
      do i2=1,jCmpMx
        if (iAOtSO(jAO+i2,j2) < 0) cycle
        MemSO1 = MemSO1+1
      end do
    end do
  end do
end do

return

end function MemSO1
