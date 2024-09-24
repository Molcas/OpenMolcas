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
! Copyright (C) 1992, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine vOO2OV(v1,nOO,v2,mOV,nD,kOV)

use Constants, only: Zero

implicit none
integer nOO, mOV, nD, iSt, iEnd, iD
integer kOV(nD)
real*8 v1(nOO,nD), v2(mOV)
interface
  subroutine vOO2OV_inner(v1,n1,v2,n2,iD)
    integer n1, n2, iD
    real*8, target :: v1(n1), v2(n2)
  end subroutine vOO2OV_inner
end interface

iEnd = 0
v2(:) = Zero
do iD=1,nD
  iSt = iEnd+1
  iEnd = iEnd+kOV(iD)
  call vOO2OV_inner(v1(:,iD),nOO,v2(iSt:iEnd),kOV(iD),iD)
end do

end subroutine vOO2OV
