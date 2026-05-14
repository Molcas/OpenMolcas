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

subroutine check_caspt2(mode)
! Check the roots to be considered in CASPT2 gradient
!
! With mode = 0, this subroutine should first try to decide the root
! for which CASPT2 density and MCLR are performed. If we do not have
! CASPT2 density, i.e. iGo /= 3, call CASPT2 using the root
! specified by ALASKA (or 'Relax CASSCF root'). Othewise (iGo = 3),
! we are going to perform MCLR next. If ALASKA has not specified
! roots, perform MCLR for the roots specified by CASPT2 ('Relax
! original root' for gradients, or that and 'Relax CASSCF root' for
! NAC). If ALASKA has specified, the roots are obtained from 'MCLR
! Root', and CASPT2 is then called to compute the density etc.
!
! With mode = 1, this subroutine obtains the character in 'MCLR Root'
! and determine the roots for NAC calculation.

use MCLR_Data, only: IRLXROOT, ISNAC, NACSTATES, OVERRIDE
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Mode
#include "warnings.h"
integer(kind=iwp) :: iGo, iRlxRootPT2, iRoot1Com, iRoot1Req, iRoot2Com, iRoot2Req, iStatus, LuINPUT, LuSpool2
character(len=128) :: FileName
character(len=72) :: Line
character(len=16) :: mstate1, StdIn
logical(kind=iwp) :: Exists, NeedGrdt
integer(kind=iwp), external :: isFreeUnit, isStructure

iRlxRoot = 0
iRlxRootPT2 = 0
call Get_iScalar('SA ready',iGo)
! Requested root for gradient
call Get_iScalar('Relax CASSCF root',iRlxRoot)
! iGo=3 means that CASPT2 density has been computed
! Check the root of the density
if (iGo == 3) call Get_iScalar('Relax original root',iRlxRootPT2)

! Check this is NAC or not
isNAC = .false.
call Get_cArray('MCLR Root',mstate1,16)
if (index(mstate1,'@') /= 0) then
  ! Requested root for NAC by ALASKA
  read(mstate1,'(1X,I7,1X,I7)') NACStates(1),NACStates(2)
  if (NACStates(1) /= 0) isNAC = .true.
  if (NACStates(1) == 0) iRlxRoot = NACStates(2)
else if ((iGo == 3) .and. (iRlxRoot /= iRlxRootPT2)) then
  ! This means CASPT2 density has been computed for the states
  ! specified by the NAC option in &CASPT2 (either specified by
  ! the original input or the call below).
  ! In this case, perform MCLR anyway(?)
  ! The states can be different from those ALASKA requests.
  ! If different, ALASKA will call MCLR then CASPT2 again with
  ! the correct states.
  NACStates(1) = iRlxRoot
  NACStates(2) = iRlxRootPT2
  isNAC = .true.
  override = .true.
end if

! With mode = 1, just set NACStates
if (mode == 1) return

!write(u6,*) 'isnac = ',isnac
!if (isnac) then
!  write(u6,*) 'requested NAC:',nacstates(1),nacstates(2)
!  write(u6,*) 'computed  NAC:',irlxroot,irlxrootpt2
!else
!  write(u6,*) 'requested GRD:',irlxroot
!  write(u6,*) 'computed  GRD:',irlxrootpt2
!endif

! If CASPT2 density has been computed, and the root of
! the density is the desired one in ALASKA, go for MCLR
if (isNAC) then
  iRoot1req = max(NACStates(1),NACStates(2))
  iRoot2req = min(NACStates(1),NACStates(2))
  iRoot1com = max(iRlxRoot,iRlxRootPT2)
  iRoot2com = min(iRlxRoot,iRlxRootPT2)
  if ((iGo /= 0) .and. (iRoot1req == iRoot1com) .and. (iRoot2req == iRoot2com)) return
else
  if ((iGo /= 0) .and. (iRlxRoot == iRlxRootPT2)) return
end if

! Otherwise, compute CASPT2 density for the state specified by ALASKA (iRlxRoot)

LuInput = 11
LuInput = IsFreeUnit(LuInput)
call StdIn_Name(StdIn)
call Molcas_open(LuInput,StdIn)

write(LuInput,'(A)') '>ECHO OFF'
write(LuInput,'(A)') '>export MCLR_OLD_TRAP=$MOLCAS_TRAP'
write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

FileName = 'CASPTINP'
call f_inquire(Filename,Exists)

if (Exists) then
  LuSpool2 = 77
  LuSpool2 = IsFreeUnit(LuSpool2)
  call Molcas_Open(LuSpool2,Filename)

  NeedGrdt = (isStructure() /= 1)
  do
    read(LuSpool2,'(A)',iostat=istatus) Line
    !write(u6,'(a)') line
    if (istatus > 0) call Abend()
    if (istatus < 0) exit
    write(LuInput,'(A)') Line
    call UpCase(Line)
    if (Line(1:4) == 'GRDT') NeedGrdt = .false.
  end do

  if (NeedGrdt) then
    backspace LuInput
    write(LuInput,'(A)') 'GRDT'
  end if

  close(LuSpool2)
else
  write(u6,'(A)') 'CASPT2 gradient without &CASPT2?'
  write(u6,'(A)') 'this cannot happen, ig'
  call abend()
end if

FileName = 'MCLRINP'
call f_inquire(Filename,Exists)

! NAC states are obtained from "MCLR Roots"
if (Exists) then
  LuSpool2 = 77
  LuSpool2 = IsFreeUnit(LuSpool2)
  call Molcas_Open(LuSpool2,Filename)
  do
    read(LuSpool2,'(A)',iostat=istatus) Line
    !write(u6,'(a)') line
    if (istatus > 0) call Abend()
    if (istatus < 0) exit
    write(LuInput,'(A)') Line
  end do
  close(LuSpool2)
else
  write(LuInput,'(A)') ' &Mclr &End'
  write(LuInput,'(A)') 'End of Input'
end if

write(LuInput,'(A)') '>export MOLCAS_TRAP=$MCLR_OLD_TRAP'
write(LuInput,'(A)') '>ECHO ON'

!rewind(luinput)
!write(u6,*)
!write(u6,*) 'show luinput'
!write(u6,*)
!do
!  read(LuInput,'(A)',iostat=istatus) Line
!  write(u6,'(a)') line
!  if (istatus > 0) call Abend()
!  if (istatus < 0) exit
!end do

close(LuInput)

call Finish(_RC_INVOKED_OTHER_MODULE_)

return

end subroutine check_caspt2
