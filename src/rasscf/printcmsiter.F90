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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************

subroutine PrintCMSIter(iStep,Qnew,Qold,RMat,lRoots)

use CMS, only: iCMSOpt, NPosHess, LargestQaaGrad, NCMSScale

implicit none
integer iStep, lRoots
real*8 Qnew, Qold, Diff
real*8 RMat(lRoots**2)

!write(6,*) 'iteration information'
Diff = Qnew-Qold
if (iCMSOpt == 2) then

  if (lRoots == 2) then
    write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)') iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
  else
    write(6,'(6X,I4,2X,F14.8,2X,ES14.4E3)') iStep,Qnew,Diff
  end if

else

  !if (lRoots == 2) then
  !  write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)') iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
  !else
  if (NCMSScale > 0) then
    write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3,I1)') iStep,Qnew,Diff,nPosHess,LargestQaaGrad,'1E-',NCMSScale
  else
    write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3)') iStep,Qnew,Diff,nPosHess,LargestQaaGrad,'1.0'
  end if
  !end if

end if

return

end subroutine PrintCMSIter
