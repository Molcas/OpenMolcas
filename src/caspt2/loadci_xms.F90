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

      Subroutine LoadCI_XMS(Bas,Mode,nConf,nState,CI,iState,U0)

      use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: IFXMS, IFRMS
      use Constants, only: Zero

      implicit none

      character(len=1), intent(in) :: Bas
      integer(kind=iwp), intent(in) :: Mode, nConf, nState, iState
      real(kind=wp), intent(inout) :: CI(Nconf)
      real(kind=wp), intent(in) :: U0(nState,nState)

!
!     MODE=0 is equivalent to LoadCI (XMS basis)
!     MODE=1 constructs the CI vector in CASSCF basis (back-transformed)
!     CSF in natural (Bas=N) or quasi-canonical (Bas=C) orbital basis
!
      If (Bas == 'N' .or. Bas == 'n') Then
        Call READ_CI(iState,IDCIEX,SIZE(IDCIEX))
!       ID = IDCIEX(1) !! natural
      Else If (Bas == 'C' .or. Bas == 'c') Then
        Call READ_CI(iState,IDTCEX,SIZE(IDTCEX))
!       ID = IDTCEX(1) !! quasi-canonical
      ELse
        write (u6,*)'the first argument in LoadCI_XMS should be either',&
     &              'N (natural) or C (quasi-canonical)'
        call abend()
      End If

      Contains

      Subroutine READ_CI(iState,IDEX,nIDEX)
      implicit none
      integer(kind=iwp), intent(in) :: iState, nIDEX, IDEX(nIDEX)

      integer(kind=iwp) :: ID, I
      real(kind=wp),allocatable :: WRK(:)

      If (Mode == 0 .or. (.not.IFXMS .and. .not.IFRMS)) Then
        ID = IDEX(iState)
        call ddafile(LUCIEX,2,CI,Nconf,ID)
      Else If (Mode == 1) Then
        CI(1:nconf) = Zero
        call mma_allocate(WRK,nconf,Label='WRK')
        do I = 1, nState
          ID = IDEX(I)
          call ddafile(LUCIEX,2,WRK,Nconf,ID)
          CI(1:nconf) = CI(1:nconf) + U0(iState,I)*WRK(1:nconf)
        end do
        call mma_deallocate(WRK)
      End If
      End Subroutine READ_CI

      End Subroutine LoadCI_XMS
