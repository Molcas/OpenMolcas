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
! Copyright (C) 2009, Roland Lindh                                     *
!               2010, Mickael G. Delcey                                *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

real*8 function dmwdot(nAt,mAt,A,B)

use Basis_Info

implicit real*8(a-h,o-z)
integer nAt, mAt
real*8 A(3,nAt), B(3,nAt)
logical Found
#include "real.fh"
#include "stdalloc.fh"
real*8 TMass
real*8, dimension(:), allocatable :: W

!***********************************************************************
!                                                                      *
!      Object : compute the weighted dot product of two vectors        *
!                                                                      *
!***********************************************************************
tmp = Zero
TMass = Zero
call Qpg_dArray('Weights',Found,nData)
if (Found .and. (nData >= mAt)) then
  call mma_allocate(W,nData,label='W')
  call Get_dArray('Weights',W,nData)
else
  call SysAbendMsg('dmwdot','No or wrong weights were found in the RUNFILE.','')
end if
iAt = 0
do iCnttp=1,nCnttp
  if (.not. (dbsc(iCnttp)%pChrg .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      iAt = iAt+1
      Fact = dble(iDeg(A(1,iAt)))
      xMass = Fact*W(iAt)
      TMass = TMass+xMass
      do i=1,3
        tmp = tmp+xMass*A(i,iAt)*B(i,iAt)
      end do
    end do
  end if
end do
call mma_deallocate(W)
dmwdot = tmp/TMass

return

end function dmwdot
