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

! This file contains simple codes called in CMSNewton. Complicated ones
! are written in files with the name as the subroutine name.

      Subroutine PrintCMSIter(iStep,Qnew,Qold,RMat,lRoots)
      use CMS, only: iCMSOpt,NPosHess,LargestQaaGrad,NCMSScale
      Implicit None
      INTEGER iStep,lRoots
      Real*8 Qnew,Qold,Diff
      Real*8 RMat(lRoots**2)

!      write(6,*) 'iteration information'
      Diff=Qnew-Qold
      IF(iCMSOpt.eq.2) THEN


       If(lRoots.eq.2) Then
        write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')                 &
     &  iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
       Else
         write(6,'(6X,I4,2X,F14.8,2X,ES14.4E3)')                        &
     &   iStep, Qnew,Diff
       End If


      ELSE


!       If(lRoots.eq.2) Then
!        write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
!     &  iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
!       Else
        if (NCMSScale.gt.0) then
      write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3,'//  &
     &        'I1)')                                                    &
     &   iStep, Qnew,Diff,nPosHess,LargestQaaGrad,'1E-',NCMSScale
        else
       write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3)')  &
     &   iStep, Qnew,Diff,nPosHess,LargestQaaGrad,'1.0'
        end if
!       End If


      END IF
      RETURN
      End Subroutine
