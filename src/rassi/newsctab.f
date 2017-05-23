************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      INTEGER FUNCTION NEWSCTAB(MINOP,MAXOP,MLTPL,MS2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ASPIN,BSPIN
C     INTEGER ID,ICP
      PARAMETER(ASPIN=1, BSPIN=0)
      INTEGER UPCPL, DWNCPL
      PARAMETER (UPCPL=1,DWNCPL=0)
      PARAMETER (NULLPTR=-1)
#include "WrkSpc.fh"

      NEWSCTAB=0 ! dummy initialize
      IF(MLTPL-MS2.LT.1) GOTO 900
      IF(MLTPL+MS2.LT.1) GOTO 900
C Run through construction loop twice. First get size of table
C and total nr of transformation coefficients:
      NSPCPL=0
      NSPDET=0
      NTRANS=0
      IBLK=0
      DO IOPEN=MINOP,MAXOP
        IBLK=IBLK+1
        NCP=NGENE(IOPEN,MLTPL)
        IF(NCP.EQ.0) GOTO 100
         NA=(IOPEN+MS2)/2
         ND=NOVERM(IOPEN,NA)
         NSPCPL=NSPCPL+IOPEN*NCP
         NSPDET=NSPDET+IOPEN*ND
         NTRANS=NTRANS+ND*NCP
 100    CONTINUE
      END DO
      NBLK=IBLK
      NTAB=8+6*NBLK+NSPCPL+NSPDET
C The spincoupling table will be NTAB integers long, and there
C will be NTRANS spin-coupling coefficients (real*8).
C The table consists of 8 header words, then an array (6,NBLK)
C with pointers and sizes to the spin coupling, spin determinant
C and spin-coupling coefficient arrays, and finally the
C spin coupling and spin determinant arrays themselves.
C The transformation coefficients are real*8 data and stored
C in a separate array.
      CALL GETMEM('SpnCplTb','Allo','Inte',LTAB,NTAB)
      CALL GETMEM('SpnCplCf','Allo','Real',LTRANS,NTRANS)
      KSPCPL=9+6*NBLK
      KSPDET=KSPCPL+NSPCPL
C Table size
      IWORK(LTAB+0)=NTAB
C Table type identifier
      IWORK(LTAB+1)=47
C Spin multiplicity
      IWORK(LTAB+2)=MLTPL
C Spin projection
      IWORK(LTAB+3)=MS2
C Min and max nr of open shells
      IWORK(LTAB+4)=MINOP
      IWORK(LTAB+5)=MAXOP
C Associated workspace array for Re*8 data (transf matrices)
      IWORK(LTAB+6)=LTRANS
      IWORK(LTAB+7)=NTRANS
C Individual information for each separate nr of open shells:
      NTAB=6
      NTRANS=0
      IBLK=0
      DO IOPEN=MINOP,MAXOP
        IBLK=IBLK+1
        NCP=NGENE(IOPEN,MLTPL)
        NA=(IOPEN+MS2)/2
        NB= IOPEN-NA
CTEST      write(*,*)' In the loop. IOPEN=',IOPEN
CTEST      write(*,*)'              MLTPL=',MLTPL
CTEST      write(*,*)'                NCP=',NCP
CTEST      write(*,*)'                MS2=',MS2
CTEST      write(*,*)'                 NA=',NA
CTEST      write(*,*)'                 NB=',NB
        IF(NCP.GT.0) THEN
CPAM        CALL GETMEM('A','Chec','Dummy',LDUM,NDUM)
CTEST        write(*,*)' Call NOVERM for ND.'
          ND=NOVERM(IOPEN,NA)
CTEST        write(*,*)' Back from NOVERM. ND=',ND
CPAM        CALL GETMEM('A','Chec','Dummy',LDUM,NDUM)
          IWORK(LTAB+ 8+(IBLK-1)*6)=IOPEN
          IWORK(LTAB+ 9+(IBLK-1)*6)=NCP
          IWORK(LTAB+10+(IBLK-1)*6)=ND
C Compute spin couplings:
CPAM        CALL GETMEM('B','Chec','Dummy',LDUM,NDUM)
CTEST        write(*,*)' Call PROTOCSF. KSPCPL=',KSPCPL
          CALL PROTOCSF(IOPEN,MLTPL,NCP,IWORK(LTAB-1+KSPCPL))
CTEST        write(*,*)' Back from PROTOCSF.'
CPAM        CALL GETMEM('B','Chec','Dummy',LDUM,NDUM)
          IWORK(LTAB+11+(IBLK-1)*6)=KSPCPL
C Compute spin determinants:
CPAM        CALL GETMEM('C','Chec','Dummy',LDUM,NDUM)
CTEST        write(*,*)' Call PROTOSD.'
CTEST        write(*,*)' Assumed available start:',KSPDET
CTEST        write(*,*)'                     end:',
CTEST     &                                     KSPDET+IOPEN*ND-1
          CALL PROTOSD(NA,NB,ND,IWORK(LTAB-1+KSPDET))
CTEST        write(*,*)' Back from PROTOSD.'
CPAM        CALL GETMEM('C','Chec','Dummy',LDUM,NDUM)
          IWORK(LTAB+12+(IBLK-1)*6)=KSPDET
C Compute spin coupling coefficients:
CPAM        CALL GETMEM('D','Chec','Dummy',LDUM,NDUM)
CTEST        write(*,*)' Call PROTOT.'
          CALL PROTOT(IOPEN,ND,IWORK(LTAB-1+KSPDET),NCP,
     &                      IWORK(LTAB-1+KSPCPL),WORK(LTRANS))
CTEST        write(*,*)' Back from PROTOT.'
CTEST        write(*,*)' Spin transformation matrix at LTRANS=',LTRANS
CTEST        do id=1,nd
CTEST          write(*,'(1x,5f16.8)')
CTEST     &                (work(ltrans-1+id+nd*(icp-1)),icp=1,ncp)
CTEST        end do
CPAM        CALL GETMEM('D','Chec','Dummy',LDUM,NDUM)
          IWORK(LTAB+13+(IBLK-1)*6)=LTRANS
          KSPCPL=KSPCPL+IOPEN*NCP
          KSPDET=KSPDET+IOPEN*ND
          LTRANS=LTRANS+ND*NCP
        ELSE
CTEST        write(*,*)' ELSE clause'
          IWORK(LTAB+ 8+(IBLK-1)*6)=IOPEN
          IWORK(LTAB+ 9+(IBLK-1)*6)=0
          IWORK(LTAB+10+(IBLK-1)*6)=0
          IWORK(LTAB+11+(IBLK-1)*6)=NULLPTR
          IWORK(LTAB+12+(IBLK-1)*6)=NULLPTR
          IWORK(LTAB+13+(IBLK-1)*6)=NULLPTR
        END IF
CPAM        write(*,*)' Checking memory now...'
CPAM        CALL GETMEM('LoopEnd','Chec','Dummy',LDUM,NDUM)
      END DO
      NEWSCTAB=LTAB
      RETURN
 900  CONTINUE
      WRITE(6,*)'NewSCTab: Contradictory values of MLTPL vs. MS2.'
      WRITE(6,*)'The function was invoked with the following arguments:'
      WRITE(6,'(1X,A,I9)')' MINOP:',MINOP
      WRITE(6,'(1X,A,I9)')' MAXOP:',MAXOP
      WRITE(6,'(1X,A,I9)')' MLTPL:',MLTPL
      WRITE(6,'(1X,A,I9)')' MS2  :',MS2
      CALL ABEND()
      END
