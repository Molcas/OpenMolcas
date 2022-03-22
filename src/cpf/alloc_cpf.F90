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

subroutine ALLOC_CPF()

use cpf_global, only: IFIRST, ILIM, IPASS, IRC, IROW, JBUF, JMAX, KBUF, KBUFF1, LBUF, LIC, LN, MADR, MAX11, MX1, MX2, NNS, NORBT, &
                      NOV, NOV1, NSYM, NTIBUF, NTMAX, NVIR, NVIRT, NVMAX, NVT5
use guga_util_global, only: IAD10
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6, RtoI

implicit none
integer(kind=iwp) :: I, IPF, IPOF(9), ISTOP, J, JBUF1, KBUF1, LBUF1, LICX, LICXX, LPERMX, LSTO3, LSTO4, NOB2, NOT2, NOVT, NVT

ISTOP = 0
MX1 = 0
MX2 = 0
NVMAX = 0
do I=1,NSYM
  call IPO_CPF(IPOF,NVIR,MUL,NSYM,I,-1)
  if (IPOF(NSYM+1) > MX1) MX1 = IPOF(NSYM+1)
  if (NVIR(I) > NVMAX) NVMAX = NVIR(I)
  do J=1,NSYM
    IPF = IPOF(J+1)-IPOF(J)
    if (IPF > MX2) MX2 = IPF
  end do
end do

! DYNAMICAL ALLOCATION FOR SORTING
NVT = IROW(NVIRT+1)
! BUFFER FOR INTEGRALS
LPERMX = NTIBUF
LICX = LIC-LPERMX
! DYNAMICAL ALLOCATION FOR SORTING AIBJ
NOB2 = IROW(NORBT+1)
NOT2 = IROW(LN+1)
NOV = 3*NOT2
LSTO3 = LICX-2*NOV-3*NOB2
! LBUF1: Nr of available reals per bin
LBUF1 = LSTO3/NOV
!PAM96 LBUF = (LBUF/2)*2
!PAM96 ! *** FPS ***
!PAM96 !LBUF = (LBUF1-2)/2
! LBUF: Nr of items per bin
LBUF = (RTOI*LBUF1-2)/(RTOI+1)
if (LBUF > 998) LBUF = 998
LBUF = ((LBUF+2)/RTOI)*RTOI-2
!PAM96 write(u6,150) NOV,MADR,LBUF
if (LBUF < 20) then
  ISTOP = 3
  write(u6,*) 'ALLOC_CPF: Impossibly small buffers, too many bins,'
  write(u6,*) 'for sorting AIBJ. Program will have to stop.'
end if
! MAXIMUM LENGTH OF HDIAG
MAX11 = max(NVT,IRC(1))
if (MAX11 < IRC(1)) MAX11 = IRC(1)
NOV1 = NOV
! DYNAMICAL ALLOCATION FOR SORTING ABCD
JBUF = 1
NOVT = 0
if (IFIRST == 0) then
  IPASS = 1
  do
    NVT5 = (NVT-1)/IPASS+1
    LICXX = LICX-KBUFF1
    LSTO4 = LICXX-2*NVT5
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
  !PAM96 write(u6,150) NOVT,MADR,JBUF
  !PAM96 write(u6,160) IPASS
  if (JBUF < 20) then
    ISTOP = 3
    write(u6,*) 'ALLOC_CPF: Impossibly small buffers, too many bins,'
    write(u6,*) 'for sorting ABCD. Program will have to stop.'
  end if
  ! DYNAMICAL ALLOCATION FOR SORTING ABCI
  NOV = LN*NVIRT+1
else
  NVT5 = 0
  NOV = 1
end if
LSTO4 = LICX-2*NOV
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
!PAM96 write(u6,150) NOVT,MADR,KBUF
if (KBUF < 20) then
  ISTOP = 3
  write(u6,*) 'ALLOC_CPF: Impossibly small buffers, too many bins,'
  write(u6,*) 'for sorting ABCI. Program will have to stop.'
end if
if (NOVT >= MADR) then
  ISTOP = 2
  write(u6,*) 'ALLOC_CPF: Too much storage needed.'
  write(u6,'(1X,A,2I10)') 'NOVT,MADR:',NOVT,MADR
end if

JMAX = IAD10(1)
if (IFIRST /= 0) JMAX = 0
NTMAX = 0
do I=1,NSYM
  if (NVIR(I) > NTMAX) NTMAX = NVIR(I)
  if (NNS(I) > NTMAX) NTMAX = NNS(I)
end do
if (IRC(ILIM) > NTMAX) NTMAX = IRC(ILIM)
if (ISTOP /= 0) then
  write(u6,*) 'ALLOC: Too little memory available.'
  write(u6,*) 'Program stops here.'

  call Abend()
end if

return

!PAM96 150 format(6X,'NUMBER OF CHAINS ON DRUM',I7,/,6X,'PRESENT LIMIT',I18,/,6X,'BUFFERT FOR SORTING',I13,/,6X,'PRESENT LIMIT', &
!PAM96            16X,'20')
!PAM96 160 format(6X,'NUMBER OF PASSES',I15)

end subroutine ALLOC_CPF
