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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
      SUBROUTINE ProtoSD(NPELA,NPELB,NPSDSZ,IPSDMS)
      DIMENSION IPSDMS(NPELA+NPELB,NPSDSZ)
      DIMENSION ITMP(50)
      INTEGER ASPIN, BSPIN
      PARAMETER (ASPIN=1,BSPIN=0)
! Given NPELA and NPELB, returns a table of all possible
! Slater determinants with NPELA alpha electrons and NPELB
! beta electrons in NPELA+NPELB orbitals.
! The ordering of the resulting SD is consistent with the
! index function   Index=1+Sum(j) noverm(k-1,j), where the sum
! is over only the alpha spins, enumerated j=1..NPELA,
! while k is the orbital with the j-th alpha electron.

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
!PAM      write(*,*)'PROTOSD. Nr of orbitals:',NPORB
!PAM      write(*,*)'         Nr of determin:',NDET
!PAM      write(*,*)'         Need array siz:',NPORB*NDET

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
      DO J=1,NPORB
        IPSDMS(J,NDET)=BSPIN
      END DO
      DO L=1,NPELA
        IPSDMS(ITMP(L),NDET)=ASPIN
      END DO
      GOTO 10

  30  CONTINUE
      RETURN
 997  CONTINUE
      WRITE(6,*)' Serious error in PROTOSD. '//                         &
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
