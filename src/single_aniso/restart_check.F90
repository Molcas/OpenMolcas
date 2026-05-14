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
use Definitions, only: wp, iwp, u5, u6

implicit none
logical(kind=iwp), intent(out) :: Ifrestart, GRAD
integer(kind=iwp), intent(out) :: input_to_read, nT, nH, nTempMagn, nDir, nDirZee, nMult
character(len=180), intent(out) :: input_file_name
integer(kind=iwp) :: i, istatus, linenr, mg, ncut, nk
real(kind=wp) :: encut_rate, rdummy
logical(kind=iwp) :: KeyENCU, KeyERAT, KeyHEXP, KeyHINT, KeyMVEC, KeyNCUT, KeyTMAG, KeyZEEM
character(len=280) :: line, tmp

Ifrestart = .false.
input_to_read = 0
GRAD = .false.
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
do
  read(u5,'(A180)',iostat=istatus) LINE
  if (istatus < 0) then
    write(u6,*) ' -- READIN: Unexpected End of input file.'
    return
  end if
  call NORMAL(LINE)
  if (LINE(1:7) == '&SINGLE') exit
end do
LINENR = 0
do
  read(u5,'(A280)',iostat=istatus) line
  if (istatus < 0) then
    write(u6,*) ' -- READIN: Unexpected End of input file.'
    return
  end if
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle

  select case (LINE(1:4))
    case ('END ','    ')
      exit

    case ('REST')
      Ifrestart = .true.
      !KeyREST = .true.
      read(u5,*) input_to_read
      input_file_name = 'aniso.input'
#     ifdef _DEBUGPRINT_
      write(u6,*) input_to_read
#     endif
      if ((input_to_read == 2) .or. (input_to_read == 3) .or. (input_to_read == 4)) then
        backspace(u5)
        read(u5,*) input_to_read,tmp
#       ifdef _DEBUGPRINT_
        write(u6,*) tmp
#       endif
        input_file_name = trim(tmp)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'restart_check: REST, input_file_name='
        write(u6,*) input_file_name
#       endif
      end if
      LINENR = LINENR+1

    case ('DATA')
      Ifrestart = .true.
      !KeyDATA = .true.
      read(u5,*) tmp
      input_file_name = trim(tmp)
      input_to_read = 6
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: DATA, input_file_name='
      write(u6,*) input_file_name
#     endif
      LINENR = LINENR+1

    case ('TEXP')
      read(u5,*) nT
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: TEXP, nT=',nT
#     endif
      !KeyTEXP = .true.
      LINENR = LINENR+1

    case ('GRAD')
      !KeyGRAD = .true.
      GRAD = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check:  GRAD = ',GRAD
#     endif
      LINENR = LINENR+1

    case ('HEXP')
      read(u5,*) nTempMagn,(rdummy,i=1,nTempMagn)
      read(u5,*) nH
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: HEXP, nH=',nH
      write(u6,*) 'restart_check: HEXP, nTempMagn=',nTempMagn
#     endif
      KeyHEXP = .true.
      LINENR = LINENR+2

    case ('HINT')
      read(u5,*) rdummy,rdummy,nH
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: HINT, nH=',nH
#     endif
      KeyHINT = .true.
      LINENR = LINENR+1

    case ('TINT')
      read(u5,*) rdummy,rdummy,nT
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: HINT, nT=',nT
#     endif
      !KeyTINT = .true.
      LINENR = LINENR+1

    case ('TMAG')
      read(u5,*) nTempMagn
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: TMAG, nTempMagn=',nTempMagn
#     endif
      KeyTMAG = .true.
      LINENR = LINENR+1

    case ('MVEC')
      read(u5,*) nDir
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: MVEC, nDir=',nDir
#     endif
      KeyMVEC = .true.
      LINENR = LINENR+1

    case ('ZEEM')
      read(u5,*) nDirZee
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: ZEEM, nDirZee=',nDirZee
#     endif
      KeyZEEM = .true.
      LINENR = LINENR+1

    case ('MLTP')
      read(u5,*) nMult
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: MLTP, nMult=',nMult
#     endif
      !KeyMLTP = .true.
      LINENR = LINENR+1

    case ('NCUT')
      read(u5,*) NCUT
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: NCUT, NCUT=',NCUT
#     endif
      KeyNCUT = .true.
      LINENR = LINENR+1

    case ('ENCU')
      read(u5,*) NK,MG
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: ENCU, NK, MG=',NK,MG
#     endif
      LINENR = LINENR+1
      KeyENCU = .true.

    case ('ERAT')
      read(u5,*) encut_rate
#     ifdef _DEBUGPRINT_
      write(u6,*) 'restart_check: ERAT, encut_rate=',encut_rate
#     endif
      KeyERAT = .true.
      LINENR = LINENR+1

  end select
end do

#include "macros.fh"
unused_var(rdummy)

! These variables are read, but never used for anything
unused_var(nk)
unused_var(mg)
unused_var(ncut)
unused_var(encut_rate)

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

return

end subroutine restart_check
