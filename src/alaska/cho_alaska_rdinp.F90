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
!****************************************************************
!
! Purpose: Read and process input for Cholesky section in alaska.
!
!****************************************************************

implicit real*8(A-H,O-Z)
#include "exterm.fh"
character*180 KWord, Key, Get_Ln
external Get_Ln
character*16 SECNAM
parameter(SECNAM='CHO_ALASKA_INPUT')
real*8 dmpK
integer nScreen
#include "chotime.fh"

! Set defaults

dmpK = 1.0d0
dmpK_default = dmpK
nScreen = 10

! Process the input
1000 continue
Key = Get_Ln(LuSpool)
Kword = Key
call UpCase(Kword)
if (KWord(1:1) == '*') Go To 1000
if (KWord(1:4) == '') Go To 1000
if (KWord(1:4) == 'DMPK') Go To 100
if (KWord(1:4) == 'SCRN') Go To 110
if (KWord(1:4) == 'TIMI') Go To 120
if (KWord(1:4) == 'ENDC') Go To 998
if (KWord(1:4) == 'END ') Go To 998
if (KWord(1:4) == 'ENDO') Go To 998
!                                                             *
!*** DMPK *****************************************************
!                                                             *
100 continue
read(LuSpool,*,err=210,end=200) dmpK
if (dmpK < 0.0d0) then
  write(6,*) 'OBS! Specified DMPK value is negative.'
  write(6,*) 'Restoring Default!'
  dmpK = dmpK_default
end if
Go To 1000
!                                                             *
!*** SCRN *****************************************************
!                                                             *
110 continue
read(LuSpool,*,err=210,end=200) nScreen
Go To 1000
!                                                             *
!*** TIMI *****************************************************
!                                                             *
120 continue
Timings = .true.
Go To 1000
!                                                             *
!** ENDChoinput ***********************************************
!                                                             *
998 continue

return
!                                                             *
!**************************************************************
!                                                             *
call ErrTra
200 write(6,*) SECNAM,'Premature end of input file.'
call Quit_onUserError()
call ErrTra
210 write(6,*) SECNAM,'Error while reading input file.'
call Quit_onUserError()

end subroutine Cho_Alaska_RdInp
