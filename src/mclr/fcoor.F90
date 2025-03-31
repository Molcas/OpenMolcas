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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine FCOOR(LUT,COOR)
!*******************************************************************
!                                                                  *
!      Transforms a symmetry adapted gradient to unsymmetric  form *
!                                                                  *
!       Written by Anders Bernhardsson                             *
!       960427                                                     *
!                                                                  *
!*******************************************************************

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep

implicit none
#include "Molcas.fh"
real*8 A(3), COOR(3,*), B(3)
integer :: LUT
integer :: iCnt, iCnttp, iCo, ii, kop, mdc
character*(LenIn) Lab

mdc = 0

write(LUT,'(A)') '*BEGIN COORDINATES'
write(LUT,'(A)') '*LABEL COORDINATES CHARGE'
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    call dcopy_(3,Coor(1,mdc),1,A,1)
    do iCo=0,nIrrep/dc(mdc)%nStab-1
      kop = dc(mdc)%iCoSet(iCo,0)
      call OA(kOp,A,B)
      ii = nint(dbsc(icnttp)%Charge)
      Lab = dc(mdc)%LblCnt(1:LenIn)
      call setLab(Lab,ico)
      write(LUT,'(1X,A,1X,3F20.10,1X,I3)') Lab,B(1:3),ii
    end do
  end do
end do
write(LUT,'(A)') '*END COORDINATES'

return

end subroutine FCOOR
