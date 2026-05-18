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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      Subroutine GradPrep(nState,UEFF,VECROT)

      use caspt2_global, only: iRoot1, iRoot2, jStLag
      use caspt2_module, only: IFMSCOUP, JSTATE
      use Constants, only: One, Half
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: nState
      real(kind=wp), intent(in) :: UEFF(nState,nState)
      real(kind=wp), intent(inout) :: VECROT(nState)

      real(kind=wp) :: TMP
      integer(kind=iwp) :: iState

!     If i  = j, UIi*dHij/dx*UJj
!     If i \= j, (UIi*UJj+UJi*UIj)*dHij/dx*0.5

      !! Construct the rotation vector
      If (IFMSCOUP) Then
        Do iState = 1, nState
          TMP = UEFF(iState,iRoot1)*UEFF(jState,iRoot2)                 &
     &        + UEFF(iState,iRoot2)*UEFF(jState,iRoot1)
          VECROT(iState) = TMP*Half
        End Do
      Else
!       write(u6,*) 'jState in gradprep: ',jstate
        VECROT(jState) = One
      End If
      jStLag = jState

      End Subroutine GradPrep
