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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine ALSO(ISTOP,LPERMA,IRC1,ISMAX)

use cpf_global, only: IFIRST, IPASS, IROW, JBUF, KBUF, KBUFF1, LBUF, LIC, LN, LW, MADR, NORBT, NTIBUF, NVIRT
use Definitions, only: iwp, u6, RtoI

implicit none
integer(kind=iwp) :: ISTOP, LPERMA, IRC1, ISMAX
integer(kind=iwp) :: INOV, INVT5, JBUF1, KBUF1, LBUF1, LICX, LICXX, LIM, LIMT, LPERMX, LSTO3, LSTO4, MAX11, NOB2, NOT2, NOV, NOVT, &
                     NVT, NVT5

! DYNAMICAL ALLOCATION FOR SORTING
! ADDRESSES LW(11)-LW(25)
NVT = IROW(NVIRT+1)
! BUFFER FOR INTEGRALS AT LPERMA = LW(7) ON FPS
LPERMX = LPERMA+NTIBUF
LICX = LIC-LPERMX
! DYNAMICAL ALLOCATION FOR SORTING AIBJ
NOB2 = IROW(NORBT+1)
NOT2 = IROW(LN+1)
NOV = 3*NOT2
INOV = NOV
LSTO3 = LICX-2*INOV-3*NOB2
! LBUF1: Nr of available reals per bin
LBUF1 = LSTO3/NOV
!PAM96 LBUF = (LBUF/2)*2
!PAM96 ! *** FPS ***
!PAM96 !LBUF = (LBUF1-2)/2
! LBUF: Nr of items per bin
LBUF = (RTOI*LBUF1-2)/(RTOI+1)
if (LBUF > 998) LBUF = 998
LBUF = ((LBUF+2)/RTOI)*RTOI-2
LBUF1 = (LBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!PAM96 write(u6,150) NOV,MADR,LBUF
if (LBUF < 20) then
  ISTOP = 3
  write(u6,*) 'ALSO: Impossibly small buffers, too many bins,'
  write(u6,*) 'for sorting AIBJ. Program will have to stop.'
end if
! SORTING AREA , BUFOUT AND INDOUT
! ALSO HDIAG IN DIAG
LW(11) = LPERMX
! CORE ADDRESSES , ICAD
!RL LW(12) = LW(11)+NOV*(2*LBUF+2)
!PAM96 LW(12) = LW(11)+NOV*(LBUF+LBUF/2+2)
LW(12) = LW(11)+NOV*((LBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
! BUFFER-COUNTER , IBUFL
LW(13) = LW(12)+INOV
! FOCK-MATRIX
LW(14) = LW(13)+INOV
! MAXIMUM LENGTH OF HDIAG
MAX11 = NVT
if (MAX11 < IRC1) MAX11 = IRC1
if (LW(14) < LW(11)+MAX11) LW(14) = LW(11)+MAX11
! IIJJ-INTEGRALS
LW(15) = LW(14)+NOB2
! IJIJ-INTEGRALS
LW(16) = LW(15)+NOB2
LIM = LW(16)+NOB2-1
!PAM96 write(u6,404) LIM
if (LIM > LIC) then
  ISTOP = 1
  write(u6,*) 'ALSO: Too much storage needed for AIBJ.'
  write(u6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR SORTING ABCD
JBUF = 1
NOVT = 0
if (IFIRST == 0) then
  IPASS = 1
  do
    NVT5 = (NVT-1)/IPASS+1
    INVT5 = NVT5
    LICXX = LICX-KBUFF1
    LSTO4 = LICXX-2*INVT5
    !RL JBUF1 = LSTO4/NVT5
    ! JBUF1: Nr of available reals per bin
    JBUF1 = LSTO4/NVT5-1
    !PAM96 JBUF = 2*(JBUF1-1)/3
    !PAM96 JBUF = (JBUF/2)*2
    !PAM96 ! *** FPS ***
    !PAM96 !RL JBUF = (JBUF1-2)/2
    ! JBUF: Nr of items per bin
    JBUF = (RTOI*JBUF1-2)/(RTOI+1)
    if (JBUF > 800) exit
    IPASS = IPASS+1
    if (IPASS > 5) exit
  end do
  if (JBUF > 998) JBUF = 998
  NOVT = NOV+NVT
  JBUF = ((JBUF+2)/RTOI)*RTOI-2
  JBUF1 = (JBUF*(RTOI+1)+2+(RTOI-1))/RTOI
  !PAM96 write(u6,150) NOVT,MADR,JBUF
  !PAM96 write(u6,160) IPASS
  if (JBUF < 20) then
    ISTOP = 3
    write(u6,*) 'ALSO: Impossibly small buffers, too many bins,'
    write(u6,*) 'for sorting ABCD. Program will have to stop.'
  end if
  ! BUFACBD
  LW(96) = LPERMX
  ! SORTING AREA , BUFOUT AND INDOUT
  LW(17) = LW(96)+KBUFF1
  ! CORE ADDRESSES , ICAD
  !RL LW(18) = LW(17)+NVT5*(2*JBUF+2)
  !PAM96 LW(18) = LW(17)+NVT5*(JBUF+JBUF/2+2)
  LW(18) = LW(17)+NVT5*((JBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
  ! BUFFER-COUNTER , IBUFL
  LW(19) = LW(18)+INVT5
  LIM = LW(19)+INVT5-1
  ! ACBDS
  ! (NOTE: FIRST PART OF BUFOUT (=2*JBUF+2) USED DURING CONSTRUCTION O
  ! ACBDS AND ACBDT VECTORS)
  !RL LW(94) = LW(17)+2*JBUF+2
  !PAM96 LW(94) = LW(17)+JBUF+JBUF/2+2
  LW(94) = LW(17)+(JBUF*(RTOI+1)+2+(RTOI-1))/RTOI
  ! ACBDT
  LW(95) = LW(94)+ISMAX
  LIMT = LW(95)+ISMAX
  if (LIMT > LW(18)) then
    ISTOP = 1
    write(u6,*) 'ALSO: Too much storage needed for ABCD.'
    write(u6,'(1X,A,2I10)') 'LIMT,LW(18):',LIMT,LW(18)
  end if
  !PAM96 write(u6,402)
  !PAM96 write(u6,405) LIM
  if (LIM > LIC) then
    ISTOP = 1
    write(u6,*) 'ALSO: Too much storage needed for ABCD.'
    write(u6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
  end if
  ! DYNAMICAL ALLOCATION FOR SORTING ABCI
  NOV = LN*NVIRT+1
else
  NOV = 1
end if
INOV = NOV
LSTO4 = LICX-2*INOV
! KBUF1: Nr of available reals per bin
KBUF1 = LSTO4/NOV
!PAM96 KBUF = 2*(KBUF1-1)/3
!PAM96 KBUF = (KBUF/2)*2
!PAM96 ! *** FPS ***
!PAM96 !RL KBUF = (KBUF1-2)/2
! KBUF: Nr of items per bin
KBUF = (RTOI*KBUF1-2)/(RTOI+1)
if (KBUF > 998) KBUF = 998
NOVT = NOVT+NOV
KBUF = ((KBUF+2)/RTOI)*RTOI-2
KBUF1 = (KBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!PAM96 write(u6,150) NOVT,MADR,KBUF
if (KBUF < 20) then
  ISTOP = 3
  write(u6,*) 'ALSO: Impossibly small buffers, too many bins,'
  write(u6,*) 'for sorting ABCI. Program will have to stop.'
end if
! SORTING AREA , BUFOUT AND INDOUT
LW(20) = LPERMX
! CORE ADDRESSES , ICAD
!RL LW(21) = LW(20)+NOV*(2*KBUF+2)
!PAM96 LW(21) = LW(20)+NOV*(KBUF+KBUF/2+2)
LW(21) = LW(20)+NOV*((KBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
! BUFFER-COUNTER , IBUFL
LW(22) = LW(21)+INOV
LIM = LW(22)+max(INOV-1,25000)
! BIAC
! (NOTE: FIRST PART OF BUFOUT(=2*KBUF+2) USED WHEN READING THE SORTE
!RL LW(23) = LW(20)+2*KBUF+2
!PAM96 LW(23) = LW(20)+KBUF+KBUF/2+2
LW(23) = LW(20)+(KBUF*(RTOI+1)+2+(RTOI-1))/RTOI
! BICA
LW(24) = LW(23)+ISMAX
! BUFBI,INDBI
LW(25) = LW(24)+ISMAX
LIMT = LW(25)+KBUFF1+2
if (LIMT >= LIM) then
  ISTOP = 1
  write(u6,*) 'ALSO: Too much storage needed.'
  write(u6,'(1X,A,2I10)') 'LIMT,LIM:',LIMT,LIM
end if
if (LIM > LIC) then
  ISTOP = 1
  write(u6,*) 'ALSO: Too much storage needed for ABCI.'
  write(u6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
!PAM96 411 write(u6,410) LIM
if (NOVT >= MADR) then
  ISTOP = 2
  write(u6,*) 'ALSO: Too much storage needed.'
  write(u6,'(1X,A,2I10)') 'NOVT,MADR:',NOVT,MADR
end if

return

!PAM96 150 format(6X,'NUMBER OF CHAINS ON DRUM',I7,/,6X,'PRESENT LIMIT',I18,/,6X,'BUFFERT FOR SORTING',I13,/,6X,'PRESENT LIMIT', &
!PAM96            16X,'20')
!PAM96 160 format(6X,'NUMBER OF PASSES',I15)
!PAM96 402 format(6X,'NOT ENOUGH STORAGE IN SORTB')
!PAM96 404 format(6X,'STORAGE FOR SORTING AIBJ',I7)
!PAM96 405 format(6X,'STORAGE FOR SORTING ABCD',I7)
!PAM96 410 format(6X,'STORAGE FOR SORTING ABCI',I7)

end subroutine ALSO
