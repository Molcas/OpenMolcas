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
      SUBROUTINE PROTOCSF(NPEL,MLTPL,NPCSFSZ,IPCSFCP)
      DIMENSION IPCSFCP(NPEL,NPCSFSZ)
      INTEGER UPCPL,DWNCPL
      PARAMETER (UPCPL=1,DWNCPL=0)
C Return a table with all possible CSF's with NPEL electrons
C coupled to give a total spin multiplicity MLTPL.
C The order of the resulting CSF's is consistent with the
C index function,
C  Index=1+Sum(j) NGENE(j-1,2*S_j+2)
C where the sum is over only the up-coupled orbitals j,
C S_j is the accumulated spin, summed over orbitals <= j,
C and NGENE(N,2*S+1) is in general the number of genealogical
C couplings of N electrons to obtain spin S.

      IF(NPEL.EQ.0) RETURN
      ISP2=MLTPL-1
      IF(ISP2.LT.0) RETURN
      IF(ISP2.GT.NPEL) RETURN
      NPELU=(NPEL+ISP2)/2
      NPELD=(NPEL-ISP2)/2
      IF(NPELU.LT.NPELD) RETURN
      IF(NPELU+NPELD.NE.NPEL) RETURN
      DO N=1,NPELU
        IPCSFCP(N,1)=UPCPL
      END DO
      IF(NPELU.EQ.NPEL) RETURN
      DO N=NPELU+1,NPEL
        IPCSFCP(N,1)=DWNCPL
      END DO

      NPCSF=NGENE(NPEL,MLTPL)
      IF(NPCSF.GT.NPCSFSZ) GOTO 998

      NPCSF=1
      IF(NPEL.LE.2) RETURN
  10  CONTINUE
      N=0
      NU=0
  20  CONTINUE
      N=N+1
      IF(N.GT.NPEL) GOTO 30
      IF(IPCSFCP(N,NPCSF).EQ.UPCPL) NU=NU+1
      ND=N-NU
      IF(IPCSFCP(N,NPCSF).EQ.UPCPL .OR. NU.EQ.ND) GOTO 20
      DO L=1,NU-1
       IPCSFCP(L,NPCSF+1)=UPCPL
      END DO
      DO L=NU,N-1
       IPCSFCP(L,NPCSF+1)=DWNCPL
      END DO
      IPCSFCP(N,NPCSF+1)=UPCPL
      DO L=N+1,NPEL
       IPCSFCP(L,NPCSF+1)=IPCSFCP(L,NPCSF)
      END DO
      NPCSF=NPCSF+1
      GOTO 10

  30  CONTINUE
      RETURN
 998  CONTINUE
      WRITE(6,*)' Too small space allocated in PROTOCSF. Input:'
      WRITE(6,'(1x,a,3i6)')' NPEL,MLTPL,NPCSFSZ:',NPEL,MLTPL,NPCSFSZ
      WRITE(6,'(1x,a,i12)')' Required NPCSFSZ is',NPCSF
      CALL ABEND()
      END
