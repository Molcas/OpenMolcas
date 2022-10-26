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

subroutine Cho_Alaska_RdInp(LuSpool)
!***********************************************************************
!
! Purpose: Read and process input for Cholesky section in alaska.
!
!***********************************************************************

use RI_glob, only: dmpK, nScreen
use Constants, only: One, Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
#include "chotime.fh"
integer(kind=iwp) :: istatus
real(kind=wp) :: dmpK_default
character(len=180) :: KWord, Key
character(len=*), parameter :: SECNAM = 'CHO_ALASKA_INPUT'
character(len=180), external :: Get_Ln

! Set defaults

dmpK = One
dmpK_default = dmpK
nScreen = 10

! Process the input
do
  Key = Get_Ln(LuSpool)
  Kword = Key
  call UpCase(Kword)
  if (KWord(1:1) == '*') cycle
  if (KWord(1:4) == '') cycle
  select case (KWord(1:4))
    case ('DMPK')
      !                                                                *
      !*** DMPK ********************************************************
      !                                                                *
      read(LuSpool,*,iostat=istatus) dmpK
      call Error()
      if (dmpK < Zero) then
        write(u6,*) 'OBS! Specified DMPK value is negative.'
        write(u6,*) 'Restoring Default!'
        dmpK = dmpK_default
      end if
    case ('SCRN')
      !                                                                *
      !*** SCRN ********************************************************
      !                                                                *
      read(LuSpool,*,iostat=istatus) nScreen
      call Error()
    case ('TIMI')
      !                                                                *
      !*** TIMI ********************************************************
      !                                                                *
      Timings = .true.
    case ('ENDC')
      !                                                                *
      !** ENDChoinput **************************************************
      !                                                                *
      exit
    case ('END ')
      exit
    case ('ENDO')
      exit
    case default
  end select
end do

return

contains

subroutine Error()
  if (istatus > 0) then
    write(u6,*) SECNAM,'Premature end of input file.'
    call Quit_onUserError()
  else if (istatus > 0) then
    write(u6,*) SECNAM,'Error while reading input file.'
    call Quit_onUserError()
  end if
end subroutine Error

end subroutine Cho_Alaska_RdInp
