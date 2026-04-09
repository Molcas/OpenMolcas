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
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This program computes start orbitals for use in SCF/RASSCF etc.      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Implementation notes:                                                *
!                                                                      *
!***********************************************************************

subroutine guessorb(iReturn,StandAlone)

use GuessOrb_Global, only: nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iReturn
logical(kind=iwp), intent(in) :: StandAlone
integer(kind=iwp) :: iRC, iUHF

!----------------------------------------------------------------------*
! Prologue                                                             *
!----------------------------------------------------------------------*
iReturn = 0
call InitGO()
if (StandAlone) call InpCtl_GuessOrb()
!----------------------------------------------------------------------*
! Select method to be used.                                            *
!----------------------------------------------------------------------*
call cre_gsswfn()
call FckByInt(iRC,StandAlone)
!if (iRC /= 0) then
if (.false.) then
  if (nSym == 1) then
    call Fmod1n()
  else
    call Fmod1s()
  end if
end if
call cls_gsswfn()
!----------------------------------------------------------------------*
! Produce MOLDEN input                                                 *
!----------------------------------------------------------------------*
iUHF = 0
if (iRC == 0) call Molden_Interface(iUHF,'GSSORB','MD_GSS')
!if (iRC == 0) call grid_driver(-1,'SEWARD','GSSORB',iRc)
!----------------------------------------------------------------------*
! Epilogue                                                             *
!----------------------------------------------------------------------*
if (StandAlone) call FastIO('STATUS')
iReturn = 0
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

end subroutine guessorb
