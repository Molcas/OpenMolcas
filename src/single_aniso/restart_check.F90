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

subroutine restart_check(Ifrestart,input_to_read,input_file_name,nT,nH,nTempMagn,nDir,nDirZee,nMult,GRAD)
! this routine looks into the file "single_aniso.input" for the "RESTart" keyword

use Constants, only: Zero
use Definitions, only: wp, u5, u6

implicit none
integer :: linenr, input_to_read, nT, nH, nTempMagn
integer :: nDir, nDirZee, nMult, i
logical :: Ifrestart
logical :: GRAD
real(kind=wp) :: rdummy
character(len=280) :: line, tmp
character(len=180) :: input_file_name
integer :: ncut, nk, mg
real(kind=wp) :: encut_rate
logical :: KeyHEXP, KeyHINT, KeyTMAG, KeyMVEC, KeyZEEM, KeyNCUT, KeyENCU, KeyERAT
!logical :: KeyREST, KeyTEXP, KeyTINT, KeyMLTP, KeyGRAD, KeyDATA
logical :: DBG

DBG = .false.

nH = 0
nT = 0
nMult = 0
nDirZee = 0
nDir = 0
nk = 0
mg = 0
ncut = 0
encut_rate = Zero
nTempMagn = 0
input_file_name = 'aniso.input'
!origin_of_data_file = 'xxxxxxxx'

!KeyREST = .false.
!KeyTEXP = .false.
KeyHEXP = .false.
KeyHINT = .false.
!KeyTINT = .false.
KeyTMAG = .false.
KeyMVEC = .false.
KeyZEEM = .false.
!KeyMLTP = .false.
KeyNCUT = .false.
KeyENCU = .false.
KeyERAT = .false.
!KeyGRAD = .false.
!KeyDATA = .false.

!=========== End of default settings====================================
rewind(u5)
50 read(u5,'(A180)',end=998) LINE
call NORMAL(LINE)
if (LINE(1:7) /= '&SINGLE') go to 50
LINENR = 0
100 read(u5,'(A280)',end=998) line
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') go to 100
if (LINE == ' ') go to 100

if ((LINE(1:4) /= 'REST') .and. (LINE(1:4) /= 'TEXP') .and. (LINE(1:4) /= 'HEXP') .and. (LINE(1:4) /= 'END ') .and. &
    (LINE(1:4) /= '    ') .and. (LINE(1:4) /= 'HINT') .and. (LINE(1:4) /= 'TINT') .and. (LINE(1:4) /= 'TMAG') .and. &
    (LINE(1:4) /= 'MVEC') .and. (LINE(1:4) /= 'ZEEM') .and. (LINE(1:4) /= 'MLTP') .and. (LINE(1:4) /= 'NCUT') .and. &
    (LINE(1:4) /= 'ENCU') .and. (LINE(1:4) /= 'ERAT') .and. (LINE(1:4) /= 'GRAD') .and. (LINE(1:4) /= 'DATA')) go to 100
if ((LINE(1:4) == 'END ') .or. (LINE(1:4) == '    ')) go to 200

if (line(1:4) == 'REST') then
  Ifrestart = .true.
  !KeyREST = .true.
  read(u5,*) input_to_read
  input_file_name = 'aniso.input'
  if (DBG) write(u6,*) input_to_read
  if ((input_to_read == 2) .or. (input_to_read == 3) .or. (input_to_read == 4)) then
    backspace(u5)
    read(u5,*) input_to_read,tmp
    if (DBG) write(u6,*) tmp
    input_file_name = trim(tmp)
    if (DBG) write(u6,*) 'restart_check: REST, input_file_name='
    if (DBG) write(u6,*) input_file_name
  end if
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'DATA') then
  Ifrestart = .true.
  !KeyDATA = .true.
  read(u5,*) tmp
  input_file_name = trim(tmp)
  input_to_read = 6
  if (DBG) write(u6,*) 'restart_check: DATA, input_file_name='
  if (DBG) write(u6,*) input_file_name
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'TEXP') then
  read(u5,*) nT
  if (DBG) write(u6,*) 'restart_check: TEXP, nT=',nT
  !KeyTEXP = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'GRAD') then
  !KeyGRAD = .true.
  GRAD = .true.
  if (DBG) write(u6,*) 'restart_check:  GRAD = ',GRAD
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'HEXP') then
  read(u5,*) nTempMagn,(rdummy,i=1,nTempMagn)
  read(u5,*) nH
  if (DBG) write(u6,*) 'restart_check: HEXP, nH=',nH
  if (DBG) write(u6,*) 'restart_check: HEXP, nTempMagn=',nTempMagn
  KeyHEXP = .true.
  LINENR = LINENR+2
  go to 100
end if

if (line(1:4) == 'HINT') then
  read(u5,*) rdummy,rdummy,nH
  if (DBG) write(u6,*) 'restart_check: HINT, nH=',nH
  KeyHINT = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'TINT') then
  read(u5,*) rdummy,rdummy,nT
  if (DBG) write(u6,*) 'restart_check: HINT, nT=',nT
  !KeyTINT = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'TMAG') then
  read(u5,*) nTempMagn
  if (DBG) write(u6,*) 'restart_check: TMAG, nTempMagn=',nTempMagn
  KeyTMAG = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'MVEC') then
  read(u5,*) nDir
  if (DBG) write(u6,*) 'restart_check: MVEC, nDir=',nDir
  KeyMVEC = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'ZEEM') then
  read(u5,*) nDirZee
  if (DBG) write(u6,*) 'restart_check: ZEEM, nDirZee=',nDirZee
  KeyZEEM = .true.
  LINENR = LINENR+1
  go to 100
end if

if (line(1:4) == 'MLTP') then
  read(u5,*) nMult
  if (DBG) write(u6,*) 'restart_check: MLTP, nMult=',nMult
  !KeyMLTP = .true.
  LINENR = LINENR+1
  go to 100
end if

if (LINE(1:4) == 'NCUT') then
  read(u5,*) NCUT
  if (DBG) write(u6,*) 'restart_check: NCUT, NCUT=',NCUT
  KeyNCUT = .true.
  LINENR = LINENR+1
  go to 100
end if

if (LINE(1:4) == 'ENCU') then
  read(u5,*) NK,MG
  if (DBG) write(u6,*) 'restart_check: ENCU, NK, MG=',NK,MG
  LINENR = LINENR+1
  KeyENCU = .true.
  go to 100
end if

if (LINE(1:4) == 'ERAT') then
  read(u5,*) encut_rate
  if (DBG) write(u6,*) 'restart_check: ERAT, encut_rate=',encut_rate
  KeyERAT = .true.
  LINENR = LINENR+1
  go to 100
end if

200 continue
write(u6,'(5X,A)') 'restart_check: NO ERROR WAS LOCATED WHILE READING INPUT'

!write(u6,*) 'KeyREST=',KeyREST
!write(u6,*) 'KeyTEXP=',KeyTEXP
!write(u6,*) 'KeyHEXP=',KeyHEXP
!write(u6,*) 'KeyHINT=',KeyHINT
!write(u6,*) 'KeyTINT=',KeyTINT
!write(u6,*) 'KeyTMAG=',KeyTMAG
!write(u6,*) 'KeyMVEC=',KeyMVEC
!write(u6,*) 'KeyZEEM=',KeyZEEM
!write(u6,*) 'KeyMLTP=',KeyMLTP
!write(u6,*) 'KeyNCUT=',KeyNCUT
!write(u6,*) 'KeyENCU=',KeyENCU
!write(u6,*) 'KeyERAT=',KeyERAT
!write(u6,*) 'KeyGRAD=',KeyGRAD

!write(u6,*) 'LOGLINE=',KeyTMAG .or. KeyZEEM .or. KeyMVEC .or. KeyHINT .or. KeyHEXP .or. KeyNCUT .or. KeyENCU .or. KeyERAT

if (KeyTMAG .or. KeyZEEM .or. KeyMVEC .or. KeyHINT .or. KeyHEXP .or. KeyNCUT .or. KeyENCU .or. KeyERAT) then
  if (nTempMagn == 0) nTempMagn = 1
  if (nH == 0) nH = 21
end if
if (nT == 0) nT = 301

!write(u6,*) 'nTempMagn=',nTempMagn
!write(u6,*) 'nH       =',nH
!write(u6,*) 'nT       =',nT

go to 190
!------ errors ------------------------------
998 continue
write(u6,*) ' -- READIN: Unexpected End of input file.'

190 continue
return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_real(rdummy)
#endif

end subroutine restart_check
