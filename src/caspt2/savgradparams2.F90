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
! Copyright (C) 2023, Yoshio Nishimoto                                 *
!***********************************************************************
!
! Save and restore many quantities that are used in CASPT2 gradient
! calculations using disk at present (hopefully)
! Save with Mode = 1, and restore with Mode = 2
! state-independent quantities
!
Subroutine SavGradParams2(Mode,UEFF,U0,H0,nState)
!
! It seems that values that are unchanged during the gradient loop
! have to be separately saved and restored
! If this subroutine is updated, the shift at the beginning of
! the SavGradParams subroutine should also be updated
!
  use caspt2_global, only: LUGRAD
  use definitions, only: iwp,wp
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: Energy, ERFSelf, nBTri, RFPert

  Implicit None

  integer(kind=iwp), intent(in) :: Mode, nState
  real(kind=wp)    , intent(inout) :: UEFF(nState,nState),U0(nState,nstate),H0(nSTate,nState)

  integer(kind=iwp) :: IORW,ID
  real(kind=wp), allocatable :: lTemp(:)

  !! Decide what to do
  If (Mode == 1) Then
    IORW = 1 !! Write
  Else If (Mode == 2) Then
    IORW = 2 !! Read
  End If

  ID = 0
  CALL DDAFILE(LUGRAD,IORW,ENERGY,NSTATE   ,ID)
  CALL DDAFILE(LUGRAD,IORW,UEFF  ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,U0    ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,H0    ,NSTATE**2,ID)

  if (RFpert) then
    Call mma_allocate(lTemp,NBTRI+1,Label='lTemp')
    if (Mode == 1) then
      Call Get_dScalar('RF Self Energy',lTemp(1+NBTRI))
      Call Get_dArray('Reaction field',lTemp,NBTRI)
      CALL DDAFILE(LUGRAD,IORW,lTemp,NBTRI+1,ID)
    else if (Mode == 2) then
      CALL DDAFILE(LUGRAD,IORW,lTemp,NBTRI+1,ID)
      Call Put_dScalar('RF Self Energy',lTemp(1+NBTRI))
      ERFSelf = lTemp(1+NBTRI)
      Call Put_dArray('Reaction field',lTemp,NBTRI)
    end if
    Call mma_deallocate(lTemp)
  end if

End Subroutine SavGradParams2
