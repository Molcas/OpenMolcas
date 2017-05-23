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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE ProtoSD(NPELA,NPELB,NPSDSZ,IPSDMS)
      DIMENSION IPSDMS(NPELA+NPELB,NPSDSZ)
      DIMENSION ITMP(50)
      INTEGER ASPIN, BSPIN
      PARAMETER (ASPIN=1,BSPIN=0)
C Given NPELA and NPELB, returns a table of all possible
C Slater determinants with NPELA alpha electrons and NPELB
C beta electrons in NPELA+NPELB orbitals.
C The ordering of the resulting SD is consistent with the
C index function   Index=1+Sum(j) noverm(k-1,j), where the sum
C is over only the alpha spins, enumerated j=1..NPELA,
C while k is the orbital with the j-th alpha electron.

      IF(NPELA.LT.0) GOTO 999
      IF(NPELB.LT.0) GOTO 999
      NPORB=NPELA+NPELB
      IF(NPORB.EQ.0) RETURN
      DO K=1,NPELA
        ITMP(K)=K
        IPSDMS(K,1)=ASPIN
      END DO
      IF(NPELA.EQ.NPORB) RETURN
      DO K=NPELA+1,NPORB
        IPSDMS(K,1)=BSPIN
      END DO
      IF(NPELA.EQ.0) RETURN

      NDET=NOVERM(NPORB,NPELA)
      IF(NDET.GT.NPSDSZ) GOTO 998
CPAM      write(*,*)'PROTOSD. Nr of orbitals:',NPORB
CPAM      write(*,*)'         Nr of determin:',NDET
CPAM      write(*,*)'         Need array siz:',NPORB*NDET

      ITMP(NPELA+1)=NPORB+1
      NDET=1

  10  CONTINUE
      K=0
  20  CONTINUE
      K=K+1
      IF(K.GT.NPELA) GOTO 30
      IF(ITMP(K+1).EQ.1+ITMP(K)) GOTO 20
      ITMP(K)=ITMP(K)+1
      DO L=1,K-1
       ITMP(L)=L
      END DO
      NDET=NDET+1
      IF(NDET.GT.NPSDSZ) GOTO 997
CPAM      CALL GETMEM('E0','Chec','Dummy',LDUM,NDUM)
      DO J=1,NPORB
        IPSDMS(J,NDET)=BSPIN
      END DO
CPAM      CALL GETMEM('E1','Chec','Dummy',LDUM,NDUM)
      DO L=1,NPELA
        IPSDMS(ITMP(L),NDET)=ASPIN
      END DO
CPAM      CALL GETMEM('E2','Chec','Dummy',LDUM,NDUM)
      GOTO 10

  30  CONTINUE
      RETURN
 997  CONTINUE
      WRITE(6,*)' Serious error in PROTOSD. '//
     &            'Too many SD''s are produced.'
      CALL ABEND()
 998  CONTINUE
      WRITE(6,*)' Too small space allocated in PROTOSD. Input:'
      WRITE(6,'(1x,a,3i6)')' NPELA,NPELB,NPSDSZ:',NPELA,NPELB,NPSDSZ
      WRITE(6,'(1x,a,i12)')' Required NPSDSZ is',NDET
      CALL ABEND()
 999  CONTINUE
      WRITE(6,*)' Invalid input to ProtoSD.'
      WRITE(6,*)'  NPELA,NPELB:',NPELA,NPELB
      CALL ABEND()
      END
