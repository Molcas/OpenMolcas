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

subroutine CNF2TXT(IFORM,NORB,NCLS,NOPN,ICONF,LENGTH,TEXT)

use definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IFORM, NORB, NCLS, NOPN
integer(kind=iwp), intent(in) :: ICONF(*)
integer(kind=iwp), intent(out) :: LENGTH
character(len=*), intent(out) :: TEXT
character(len=1) SEP
integer(kind=iwp) MXWR, NWR, IWORD, NOCC, IORB, IOCC, IREST, IW

MXWR = len(TEXT)
NWR = 1
TEXT(1:1) = '('
IWORD = 0 ! dummy initialize
if ((IFORM == 1) .or. (IFORM == 3)) then
  NOCC = NCLS+NOPN
  if (NCLS == 0) then
    NWR = 2
    TEXT(2:2) = ';'
  end if
  if (IFORM == 1) then
    do IOCC=1,NOCC
      IORB = ICONF(IOCC)
      SEP = ','
      if (IOCC == NCLS) SEP = ';'
      if (IORB < 10) then
        NWR = min(MXWR,NWR+2)
        write(TEXT(NWR-1:NWR),'(I1,A1)') IORB,SEP
      else if (IORB < 100) then
        NWR = min(MXWR,NWR+3)
        write(TEXT(NWR-2:NWR),'(I2,A1)') IORB,SEP
      else
        NWR = min(MXWR,NWR+4)
        write(TEXT(NWR-3:NWR),'(I3,A1)') IORB,SEP
      end if
    end do
  else
    do IOCC=1,NOCC
      IW = (3+IOCC)/4
      IREST = (3+IOCC)-4*IW
      if (IREST == 0) IWORD = ICONF(IW)
      IORB = mod(IWORD,256)
      IWORD = IWORD/256
      SEP = ','
      if (IOCC == NCLS) SEP = ';'
      if (IORB < 10) then
        NWR = min(MXWR,NWR+2)
        write(TEXT(NWR-1:NWR),'(I1,A1)') IORB,SEP
      else if (IORB < 100) then
        NWR = min(MXWR,NWR+3)
        write(TEXT(NWR-2:NWR),'(I2,A1)') IORB,SEP
      else
        NWR = min(MXWR,NWR+4)
        write(TEXT(NWR-3:NWR),'(I3,A1)') IORB,SEP
      end if
    end do
  end if
else if ((IFORM == 2) .or. (IFORM == 4)) then
  if (IFORM == 2) then
    do IORB=1,NORB
      IOCC = ICONF(IORB)
      NWR = min(MXWR,NWR+1)
      write(TEXT(NWR:NWR),'(I1)') IOCC
    end do
  else
    do IORB=1,NORB
      IW = (IORB+14)/15
      IREST = IORB+14-15*IW
      if (IREST == 0) IWORD = ICONF(IW)
      IOCC = mod(IWORD,4)
      IWORD = IWORD/4
      NWR = min(MXWR,NWR+1)
      write(TEXT(NWR:NWR),'(I1)') IOCC
    end do
  end if
end if
TEXT(NWR:NWR) = ')'
LENGTH = NWR

end subroutine CNF2TXT
