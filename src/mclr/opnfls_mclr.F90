!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine OpnFls_MCLR(iPL)
!***********************************************************************
!                                                                      *
!     Open files.                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MCLR_Data, only: FnMck, FnPT2, FnTemp, FnTwo, LuMck, LuTEMP, LuTwo
use input_mclr, only: ChIrr, McKinley, PT2
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPL
integer(kind=iwp) :: iRC, iOpt, iDum
logical(kind=iwp) :: DoCholesky, DoDirect, FoundTwoEls
character(len=8) :: Label

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!---  open the JOBIPH file --------------------------------------------*
!call DaName(LuJob,FnJob)
call DaName(LuTemp,FnTemp)
!---  open the ORDINT file --------------------------------------------*
call f_Inquire(FnTwo,FoundTwoEls)
call DecideOnDirect(.true.,FoundTwoEls,DoDirect,DoCholesky)
if (DoDirect) then
  write(u6,*) 'OpnFls: No direct option in MCLR'
  call Abend()
else
  if (.not. DoCholesky) then
    if (iPL >= 2) write(u6,*) 'Ordinary integral handling'
    iRc = -1
    iOpt = 0
    call OpnOrd(iRc,iOpt,FnTwo,LuTwo)
    if (iRc /= 0) then
      write(u6,*) 'OpnFls: Error opening ORDINT'
      call Abend()
    end if
  end if
end if
call f_Inquire(FnMCK,McKinley)
call f_Inquire(FnPT2,PT2)
if (McKinley) then
  !write(u6,*) 'Calculating response on perturbation from mckinley'
  iRc = -1
  iOpt = 0
  call OpnMck(iRc,iOpt,FnMck,LuMck)
  if (iRc /= 0) then
    write(u6,*) 'OpnFls: Error opening MCKINT'
    call Abend()
  end if
  iRc = -1
  idum = 0
  iOpt = 0
  Label = 'SYMOP'
  call cRdMck(irc,iopt,Label,idum,chirr(1),idum)
  if (iRc /= 0) then
    write(u6,*) 'OpnFls: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if

else if (PT2) then
  !if (iPL >= 2) write(u6,*) 'Calculating lagrange multipliers for CASPT2'
  !call DaName(LuPT2,FnPT2)
else
  if (iPL >= 2) then
    write(u6,*) 'No ',FnPT2,' or ',FNMCK,', I hope that is OK'
    write(u6,*) 'Seward mode is assumed, reading perturbation from ONEINT'
  end if
end if
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine OpnFls_MCLR
