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

subroutine ALLOC_CPF(ISMAX,LPERMA)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension IPOF(9)

ILIM = 4
if (IFIRST /= 0) ILIM = 2
!PAM96 write(6,50)
!PAM96 50 format(/,6X,'DYNAMICAL CORE STORAGE INFORMATION',/)
ISTOP = 0
MAX1 = 0
MAX2 = 0
NVMAX = 0
do I=1,NSYM
  call IPO_CPF(IPOF,NVIR,MUL,NSYM,I,-1)
  if (IPOF(NSYM+1) > MAX1) MAX1 = IPOF(NSYM+1)
  if (NVIR(I) > NVMAX) NVMAX = NVIR(I)
  do J=1,NSYM
    IPF = IPOF(J+1)-IPOF(J)
    if (IPF > MAX2) MAX2 = IPF
  end do
end do
!PAM97 call ALSO(LW,IROW,LIC,NORBT,NVIRT,LN,ISTOP,IFIRST,MADR,LPERMA,LBUF,KBUF,JBUF,IPASS,IRC(1),ISMAX,KBUFF1)
call ALSO(ISTOP,LPERMA,IRC(1),ISMAX)
! VECTORS PERMANENTLY IN CORE DURING CI ITERATIONS
! CI-VECTOR
! C
LW(26) = LPERMA
! S (SIGMA)
LW(27) = LW(26)+JSC(ILIM)
! W
LW(28) = LW(27)+JSC(ILIM)
! ITHE (TETA)
LW(29) = LW(28)+JSC(ILIM)
if ((ICPF == 1) .or. (ISDCI == 1) .or. (INCPF == 1)) LW(29) = LW(28)
! TPQ (ONE ROW)
LW(30) = LW(29)+IRC(ILIM)*IRC(ILIM)
if ((ICPF == 1) .or. (ISDCI == 1) .or. (INCPF == 1)) LW(30) = LW(29)
! ENP (NP)
LW(31) = LW(30)+IRC(ILIM)
! EPP (EPRIME)
LW(32) = LW(31)+IRC(ILIM)
! BST (STORED BIJ MATRIX)
LW(33) = LW(32)+IRC(ILIM)
! ADDRESSES NOT USED
LW(34) = LW(33)+(MAXIT+1)**2
LW(35) = LW(34)
LPERMB = LW(35)
! DYNAMICAL ALLOCATION FOR AIBJ
! MATRIX ABIJ
LW(36) = LPERMB
! MATRIX AIBJ
LW(37) = LW(36)+MAX1
! MATRIX AJBI
LW(38) = LW(37)+MAX1
! BUFIN , IBUFIN
LW(39) = LW(38)+MAX1
! A
!RL LW(40) = LW(39)+2*LBUF+2
!PAM96 LW(40) = LW(39)+LBUF+LBUF/2+2
LW(40) = LW(39)+((RTOI+1)*LBUF+2+(RTOI-1))/RTOI
! B
LW(41) = LW(40)+MAX2
! F
LW(42) = LW(41)+MAX2
! FSEC
LW(43) = LW(42)+MAX1
! ADDRESSES NOT USED
LW(44) = LW(43)+MAX1
LW(45) = LW(44)
LIM = LW(45)
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for AIBJ.'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR IJKL
! FIJKL
LW(46) = LPERMB
NIJ = IROW(LN+1)
IJKL = (NIJ*(NIJ+1))/2
! BUFIN , IBUFIN
LW(47) = LW(46)+IJKL
! ADDRESSES NOT USED
LW(48) = LW(47)+KBUFF1+2
LW(49) = LW(48)
LIM = LW(49)
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for IJKL.'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR ABCI
! BMN
LW(50) = LPERMB
! IBMN
JMAX = IAD10(1)
if (IFIRST /= 0) JMAX = 0
LW(51) = LW(50)+JMAX
! BIAC
LW(52) = LW(51)+JMAX
! BICA
LW(53) = LW(52)+ISMAX
! BUFIN
LW(54) = LW(53)+ISMAX
! ADDRESSES NOT USED
LW(55) = LW(54)+KBUFF1
LW(56) = LW(55)
LIM = LW(56)
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for ABCI.'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR ABCD
! ACBDS
LW(57) = LPERMB
! ACBDT
LW(58) = LW(57)+ISMAX
! BUFIN
LW(59) = LW(58)+ISMAX
! ADDRESSES NOT USED
LW(60) = LW(59)+KBUFF1
LW(61) = LW(60)
LIM = LW(61)
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for ABCD.'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR FIJ, AI AND AB
! FC
LW(62) = LPERMB
NOB2 = IROW(NORBT+1)
! BUFIN , IBUFIN
LW(63) = LW(62)+NOB2
! A
!RL LW(64) = LW(63)+2*LBUF+2
!PAM96 LW(64) = LW(63)+LBUF+LBUF/2+2
LW(64) = LW(63)+((RTOI+1)*LBUF+2+(RTOI-1))/RTOI
! B
LW(65) = LW(64)+MAX2
! FK IN AI AND F IN AB
LW(66) = LW(65)+MAX2
! DBK
LW(67) = LW(66)+NVMAX
LIM = LW(67)+NVMAX
LIM1 = LW(66)+MAX1
if (LIM < LIM1) LIM = LIM1
! ADDRESSES NOT USED
LW(68) = LIM
LW(69) = LW(68)
LW(70) = LW(69)
LW(71) = LW(70)
LIM = LW(71)
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for FIJ,AI,AB'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! DYNAMICAL ALLOCATION FOR NPSET
! TEMP
LW(72) = LPERMB
! ADDRESSES NOT USED
LW(73) = LW(72)+IRC(ILIM)
LW(74) = LW(73)
LIM = LW(74)
!PAM96 write(6,358) LIM
!PAM96 358 format(6X,'STORAGE FOR NPSET',I14)
! DYNAMICAL ALLOCATION FOR CPFCTL
! EPB (EPBIS)
LW(75) = LPERMB
! AP (APPRIME)
LW(76) = LW(75)+IRC(ILIM)
! BIJ
LW(77) = LW(76)+IRC(ILIM)
! CN
LW(78) = LW(77)+(MAXIT+1)**2
! TEMP1
LW(79) = LW(78)+MAXIT+1
NTMAX = 0
do I=1,NSYM
  if (NVIR(I) > NTMAX) NTMAX = NVIR(I)
  if (NNS(I) > NTMAX) NTMAX = NNS(I)
end do
if (IRC(ILIM) > NTMAX) NTMAX = IRC(ILIM)
! TEMP2
LW(80) = LW(79)+NTMAX
LIM = LW(80)+NTMAX
if (LIM > LIC) then
  ISTOP = 1
  write(6,*) 'ALLOC: Too much storage needed for CPFCTL'
  write(6,'(1X,A,2I10)') 'LIM,LIC:',LIM,LIC
end if
! NATURAL ORBITALS
! DENSITY MATRIX AT LPERMB , D
! MOLECULAR ORBITALS ALL SYMMETRIES , CM
LW(87) = LPERMB+NOB2
LCIN = 0
NBMAX = 0
do I=1,NSYM
  NBMAX = max(NBMAX,NBAS(I))
  LCIN = LCIN+NBAS(I)*NBAS(I)
end do
! NATURAL ORBITALS ONE SYMMETRY , CMO
LW(88) = LW(87)+LCIN
! DENSITY MATRIX FOR ONE SYMMETRY , DSYM
! MOLECULAR ORBITAL IN AO-BASIS for one symmetry, CAO
LW(89) = LW(88)+NBMAX*NBMAX
! OCCUPATION NUMBERS FOR ALL POSSIBLE ORBITALS
NTOT = 0
do I=1,NSYM
  NTOT = NTOT+NBAS(I)
end do
LW(90) = LW(89)+NBMAX*NBMAX
! OVERLAP MATRIX IN CHARGE
LW(91) = LW(90)+NTOT
! ADDRESSES NOT USED
LW(92) = LW(91)
LW(93) = LW(92)
LIM = LW(93)
!PAM96 write(6,349) LIM
!PAM96 349 format(6X,'STORAGE FOR DENS',I15)
if (IPRINT >= 2) then
  ! LW(94), LW(95) AND LW(96) USED IN SORTING ABCD
  write(6,450)
  write(6,451) (LW(I),I=1,96)
450 format(//,6X,'DYNAMICAL STORAGE ADDRESSES LW:',/)
451 format(6X,5I10)
end if
if (ISTOP == 0) return
write(6,*) 'ALLOC: Too little memory available.'
write(6,*) 'Program stops here.'

call Abend()

end subroutine ALLOC_CPF
