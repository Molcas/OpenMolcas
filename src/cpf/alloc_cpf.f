************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE ALLOC_CPF(ISMAX,LPERMA)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
      DIMENSION IPOF(9)
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
CPAM96      WRITE(6,50)
CPAM9650    FORMAT(/,6X,'DYNAMICAL CORE STORAGE INFORMATION',/)
      ISTOP=0
      MAX1=0
      MAX2=0
      NVMAX=0
      DO 10 I=1,NSYM
      CALL IPO(IPOF,NVIR,MUL,NSYM,I,-1)
      IF(IPOF(NSYM+1).GT.MAX1)MAX1=IPOF(NSYM+1)
      IF(NVIR(I).GT.NVMAX)NVMAX=NVIR(I)
      DO 20 J=1,NSYM
      IPF=IPOF(J+1)-IPOF(J)
      IF(IPF.GT.MAX2)MAX2=IPF
20    CONTINUE
10    CONTINUE
CPAM97      CALL ALSO(LW,IROW,LIC,NORBT,NVIRT,LN,ISTOP,IFIRST,
CPAM97     *MADR,LPERMA,LBUF,KBUF,JBUF,IPASS,IRC(1),ISMAX,KBUFF1)
      CALL ALSO(ISTOP,LPERMA,IRC(1),ISMAX)
C     VECTORS PERMANENTLY IN CORE DURING CI ITERATIONS
C     CI-VECTOR
C     C
      LW(26)=LPERMA
C     S (SIGMA)
      LW(27)=LW(26)+JSC(ILIM)
C     W
      LW(28)=LW(27)+JSC(ILIM)
C     ITHE (TETA)
      LW(29)=LW(28)+JSC(ILIM)
      IF(ICPF.EQ.1.OR.ISDCI.EQ.1.OR.INCPF.EQ.1)LW(29)=LW(28)
C     TPQ (ONE ROW)
      LW(30)=LW(29)+IRC(ILIM)*IRC(ILIM)
      IF(ICPF.EQ.1.OR.ISDCI.EQ.1.OR.INCPF.EQ.1)LW(30)=LW(29)
C     ENP (NP)
      LW(31)=LW(30)+IRC(ILIM)
C     EPP (EPRIME)
      LW(32)=LW(31)+IRC(ILIM)
C     BST (STORED BIJ MATRIX)
      LW(33)=LW(32)+IRC(ILIM)
C     ADDRESSES NOT USED
      LW(34)=LW(33)+(MAXIT+1)**2
      LW(35)=LW(34)
      LPERMB=LW(35)
C     DYNAMICAL ALLOCATION FOR AIBJ
C     MATRIX ABIJ
      LW(36)=LPERMB
C     MATRIX AIBJ
      LW(37)=LW(36)+MAX1
C     MATRIX AJBI
      LW(38)=LW(37)+MAX1
C     BUFIN , IBUFIN
      LW(39)=LW(38)+MAX1
C     A
CRL   LW(40)=LW(39)+2*LBUF+2
CPAM96      LW(40)=LW(39)+LBUF+LBUF/2+2
      LW(40)=LW(39)+((RTOI+1)*LBUF+2+(RTOI-1))/RTOI
C     B
      LW(41)=LW(40)+MAX2
C     F
      LW(42)=LW(41)+MAX2
C     FSEC
      LW(43)=LW(42)+MAX1
C     ADDRESSES NOT USED
      LW(44)=LW(43)+MAX1
      LW(45)=LW(44)
      LIM=LW(45)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for AIBJ.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     DYNAMICAL ALLOCATION FOR IJKL
C     FIJKL
      LW(46)=LPERMB
      NIJ=IROW(LN+1)
      IJKL=(NIJ*(NIJ+1))/2
C     BUFIN , IBUFIN
      LW(47)=LW(46)+IJKL
C     ADDRESSES NOT USED
      LW(48)=LW(47)+KBUFF1+2
      LW(49)=LW(48)
      LIM=LW(49)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for IJKL.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     DYNAMICAL ALLOCATION FOR ABCI
C     BMN
      LW(50)=LPERMB
C     IBMN
      JMAX=IAD10(1)
      IF(IFIRST.NE.0)JMAX=0
      LW(51)=LW(50)+JMAX
C     BIAC
      LW(52)=LW(51)+JMAX
C     BICA
      LW(53)=LW(52)+ISMAX
C     BUFIN
      LW(54)=LW(53)+ISMAX
C     ADDRESSES NOT USED
      LW(55)=LW(54)+KBUFF1
      LW(56)=LW(55)
      LIM=LW(56)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for ABCI.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     DYNAMICAL ALLOCATION FOR ABCD
C     ACBDS
      LW(57)=LPERMB
C     ACBDT
      LW(58)=LW(57)+ISMAX
C     BUFIN
      LW(59)=LW(58)+ISMAX
C     ADDRESSES NOT USED
      LW(60)=LW(59)+KBUFF1
      LW(61)=LW(60)
      LIM=LW(61)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for ABCD.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     DYNAMICAL ALLOCATION FOR FIJ, AI AND AB
C     FC
      LW(62)=LPERMB
      NOB2=IROW(NORBT+1)
C     BUFIN , IBUFIN
      LW(63)=LW(62)+NOB2
C     A
CRL   LW(64)=LW(63)+2*LBUF+2
CPAM96      LW(64)=LW(63)+LBUF+LBUF/2+2
      LW(64)=LW(63)+((RTOI+1)*LBUF+2+(RTOI-1))/RTOI
C     B
      LW(65)=LW(64)+MAX2
C     FK IN AI AND F IN AB
      LW(66)=LW(65)+MAX2
C     DBK
      LW(67)=LW(66)+NVMAX
      LIM=LW(67)+NVMAX
      LIM1=LW(66)+MAX1
      IF(LIM.LT.LIM1)LIM=LIM1
C     ADDRESSES NOT USED
      LW(68)=LIM
      LW(69)=LW(68)
      LW(70)=LW(69)
      LW(71)=LW(70)
      LIM=LW(71)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for FIJ,AI,AB'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     DYNAMICAL ALLOCATION FOR NPSET
C     TEMP
      LW(72)=LPERMB
C     ADDRESSES NOT USED
      LW(73)=LW(72)+IRC(ILIM)
      LW(74)=LW(73)
      LIM=LW(74)
CPAM96      WRITE(6,358)LIM
CPAM96358   FORMAT(6X,'STORAGE FOR NPSET',I14)
C     DYNAMICAL ALLOCATION FOR CPFCTL
C     EPB (EPBIS)
      LW(75)=LPERMB
C     AP (APPRIME)
      LW(76)=LW(75)+IRC(ILIM)
C     BIJ
      LW(77)=LW(76)+IRC(ILIM)
C     CN
      LW(78)=LW(77)+(MAXIT+1)**2
C     TEMP1
      LW(79)=LW(78)+MAXIT+1
      NTMAX=0
      DO 357 I=1,NSYM
      IF(NVIR(I).GT.NTMAX)NTMAX=NVIR(I)
      IF(NNS(I).GT.NTMAX)NTMAX=NNS(I)
357   CONTINUE
      IF(IRC(ILIM).GT.NTMAX)NTMAX=IRC(ILIM)
C     TEMP2
      LW(80)=LW(79)+NTMAX
      LIM=LW(80)+NTMAX
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALLOC: Too much storage needed for CPFCTL'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
C     NATURAL ORBITALS
C     DENSITY MATRIX AT LPERMB , D
C     MOLECULAR ORBITALS ALL SYMMETRIES , CM
      LW(87)=LPERMB+NOB2
      LCIN=0
      NBMAX=0
      DO 350 I=1,NSYM
         NBMAX = MAX(NBMAX,NBAS(I))
         LCIN=LCIN+NBAS(I)*NBAS(I)
350   CONTINUE
C     NATURAL ORBITALS ONE SYMMETRY , CMO
      LW(88)=LW(87)+LCIN
C     DENSITY MATRIX FOR ONE SYMMETRY , DSYM
C     MOLECULAR ORBITAL IN AO-BASIS for one symmetry, CAO
      LW(89)=LW(88)+NBMAX*NBMAX
C     OCCUPATION NUMBERS FOR ALL POSSIBLE ORBITALS
      NTOT=0
      DO 333 I=1,NSYM
      NTOT=NTOT+NBAS(I)
333   CONTINUE
      LW(90)=LW(89)+NBMAX*NBMAX
C     OVERLAP MATRIX IN CHARGE
      LW(91)=LW(90)+NTOT
*     Scratch for property calculation
      LW(92)=LW(91)+LCIN
C     ADDRESSES NOT USED
      LW(93)=LW(92)
      LIM=LW(93)
CPAM96      WRITE(6,349)LIM
CPAM96349   FORMAT(6X,'STORAGE FOR DENS',I15)
      IF(IPRINT.GE.2) THEN
C LW(94), LW(95) AND LW(96) USED IN SORTING ABCD
        WRITE(6,450)
        WRITE(6,451)(LW(I),I=1,96)
450     FORMAT(//,6X,'DYNAMICAL STORAGE ADDRESSES LW:',/)
451     FORMAT(6X,5I10)
      END IF
      IF(ISTOP.EQ.0)RETURN
      WRITE(6,*)'ALLOC: Too little memory available.'
      WRITE(6,*)'Program stops here.'
      CALL Abend
      END
