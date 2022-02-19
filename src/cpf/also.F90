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
      SUBROUTINE ALSO(ISTOP,LPERMA,IRC1,ISMAX)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "cpfmcpf.fh"
!     DYNAMICAL ALLOCATION FOR SORTING
!     ADDRESSES LW(11)-LW(25)
      NVT=IROW(NVIRT+1)
!     BUFFER FOR INTEGRALS AT LPERMA = LW(7) ON FPS
      LPERMX=LPERMA+NTIBUF
      LICX=LIC-LPERMX
!     DYNAMICAL ALLOCATION FOR SORTING AIBJ
      NOB2=IROW(NORBT+1)
      NOT2=IROW(LN+1)
      NOV=3*NOT2
      INOV=NOV
      LSTO3=LICX-2*INOV-3*NOB2
! LBUF1: Nr of available reals per bin
      LBUF1=LSTO3/NOV
!PAM96      LBUF=(LBUF/2)*2
!PAM96C     *** FPS ***
!PAM96C     LBUF=(LBUF1-2)/2
! LBUF: Nr of items per bin
      LBUF=(RTOI*LBUF1-2)/(RTOI+1)
      IF(LBUF.GT.998)LBUF=998
      LBUF=((LBUF+2)/RTOI)*RTOI-2
      LBUF1=(LBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!PAM96      WRITE(6,150)NOV,MADR,LBUF
!PAM96150   FORMAT(6X,'NUMBER OF CHAINS ON DRUM',I7,
!PAM96     */,6X,'PRESENT LIMIT',I18,
!PAM96     */,6X,'BUFFERT FOR SORTING',I13,
!PAM96     */,6X,'PRESENT LIMIT',16X,'20')
      IF(LBUF.LT.20) THEN
        ISTOP=3
        WRITE(6,*)'ALSO: Impossibly small buffers, too many bins,'
        WRITE(6,*)'for sorting AIBJ. Program will have to stop.'
      END IF
!     SORTING AREA , BUFOUT AND INDOUT
!     ALSO HDIAG IN DIAG
      LW(11)=LPERMX
!     CORE ADDRESSES , ICAD
!RL   LW(12)=LW(11)+NOV*(2*LBUF+2)
!PAM96      LW(12)=LW(11)+NOV*(LBUF+LBUF/2+2)
      LW(12)=LW(11)+NOV*((LBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
!     BUFFER-COUNTER , IBUFL
      LW(13)=LW(12)+INOV
!     FOCK-MATRIX
      LW(14)=LW(13)+INOV
!     MAXIMUM LENGTH OF HDIAG
      MAX11=NVT
      IF(MAX11.LT.IRC1)MAX11=IRC1
      IF(LW(14).LT.LW(11)+MAX11)LW(14)=LW(11)+MAX11
!     IIJJ-INTEGRALS
      LW(15)=LW(14)+NOB2
!     IJIJ-INTEGRALS
      LW(16)=LW(15)+NOB2
      LIM=LW(16)+NOB2-1
!PAM96      WRITE(6,404)LIM
!PAM96404   FORMAT(6X,'STORAGE FOR SORTING AIBJ',I7)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALSO: Too much storage needed for AIBJ.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
!     DYNAMICAL ALLOCATION FOR SORTING ABCD
      JBUF=1
      NOVT=0
      IF(IFIRST.NE.0)GO TO 35
      IPASS=1
110   NVT5=(NVT-1)/IPASS+1
      INVT5=NVT5
      LICXX=LICX-KBUFF1
      LSTO4=LICXX-2*INVT5
!RL   JBUF1=LSTO4/NVT5
! JBUF1: Nr of available reals per bin
      JBUF1=LSTO4/NVT5-1
!PAM96      JBUF=2*(JBUF1-1)/3
!PAM96      JBUF=(JBUF/2)*2
!PAM96C     *** FPS ***
!PAM96CRL   JBUF=(JBUF1-2)/2
! JBUF: Nr of items per bin
      JBUF=(RTOI*JBUF1-2)/(RTOI+1)
      IF(JBUF.GT.800)GO TO 120
      IPASS=IPASS+1
      IF(IPASS.GT.5)GO TO 120
      GO TO 110
120   IF(JBUF.GT.998)JBUF=998
      NOVT=NOV+NVT
      JBUF=((JBUF+2)/RTOI)*RTOI-2
      JBUF1=(JBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!PAM96      WRITE(6,150)NOVT,MADR,JBUF
!PAM96      WRITE(6,160)IPASS
!160   FORMAT(6X,'NUMBER OF PASSES',I15)
      IF(JBUF.LT.20) THEN
        ISTOP=3
        WRITE(6,*)'ALSO: Impossibly small buffers, too many bins,'
        WRITE(6,*)'for sorting ABCD. Program will have to stop.'
      END IF
!     BUFACBD
      LW(96)=LPERMX
!     SORTING AREA , BUFOUT AND INDOUT
      LW(17)=LW(96)+KBUFF1
!     CORE ADDRESSES , ICAD
!RL   LW(18)=LW(17)+NVT5*(2*JBUF+2)
!PAM96      LW(18)=LW(17)+NVT5*(JBUF+JBUF/2+2)
      LW(18)=LW(17)+NVT5*((JBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
!     BUFFER-COUNTER , IBUFL
      LW(19)=LW(18)+INVT5
      LIM=LW(19)+INVT5-1
!     ACBDS
!     (NOTE: FIRST PART OF BUFOUT (=2*JBUF+2) USED DURING CONSTRUCTION O
!     ACBDS AND ACBDT VECTORS)
!RL   LW(94)=LW(17)+2*JBUF+2
!PAM96      LW(94)=LW(17)+JBUF+JBUF/2+2
      LW(94)=LW(17)+(JBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!     ACBDT
      LW(95)=LW(94)+ISMAX
      LIMT=LW(95)+ISMAX
      IF(LIMT.GT.LW(18)) THEN
        ISTOP=1
        WRITE(6,*)'ALSO: Too much storage needed for ABCD.'
        WRITE(6,'(1X,A,2I10)')'LIMT,LW(18):',LIMT,LW(18)
      END IF
!PAM96      WRITE(6,402)
!PAM96402   FORMAT(6X,'NOT ENOUGH STORAGE IN SORTB')
!PAM96401   WRITE(6,405)LIM
!PAM96405   FORMAT(6X,'STORAGE FOR SORTING ABCD',I7)
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALSO: Too much storage needed for ABCD.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
!     DYNAMICAL ALLOCATION FOR SORTING ABCI
      NOV=LN*NVIRT+1
35    IF(IFIRST.NE.0)NOV=1
      INOV=NOV
      LSTO4=LICX-2*INOV
! KBUF1: Nr of available reals per bin
      KBUF1=LSTO4/NOV
!PAM96      KBUF=2*(KBUF1-1)/3
!PAM96      KBUF=(KBUF/2)*2
!PAM96C     *** FPS ***
!PAM96CRL   KBUF=(KBUF1-2)/2
! KBUF: Nr of items per bin
      KBUF=(RTOI*KBUF1-2)/(RTOI+1)
      IF(KBUF.GT.998)KBUF=998
      NOVT=NOVT+NOV
      KBUF=((KBUF+2)/RTOI)*RTOI-2
      KBUF1=(KBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!PAM96      WRITE(6,150)NOVT,MADR,KBUF
      IF(KBUF.LT.20) THEN
        ISTOP=3
        WRITE(6,*)'ALSO: Impossibly small buffers, too many bins,'
        WRITE(6,*)'for sorting ABCI. Program will have to stop.'
      END IF
!     SORTING AREA , BUFOUT AND INDOUT
      LW(20)=LPERMX
!     CORE ADDRESSES , ICAD
!RL   LW(21)=LW(20)+NOV*(2*KBUF+2)
!PAM96      LW(21)=LW(20)+NOV*(KBUF+KBUF/2+2)
      LW(21)=LW(20)+NOV*((KBUF*(RTOI+1)+2+(RTOI-1))/RTOI)
!     BUFFER-COUNTER , IBUFL
      LW(22)=LW(21)+INOV
      LIM=LW(22)+MAX(INOV-1,25000)
!     BIAC
!     (NOTE: FIRST PART OF BUFOUT(=2*KBUF+2) USED WHEN READING THE SORTE
!RL   LW(23)=LW(20)+2*KBUF+2
!PAM96      LW(23)=LW(20)+KBUF+KBUF/2+2
      LW(23)=LW(20)+(KBUF*(RTOI+1)+2+(RTOI-1))/RTOI
!     BICA
      LW(24)=LW(23)+ISMAX
!     BUFBI,INDBI
      LW(25)=LW(24)+ISMAX
      LIMT=LW(25)+KBUFF1+2
      IF(LIMT.GE.LIM) then
        ISTOP=1
        WRITE(6,*)'ALSO: Too much storage needed.'
        WRITE(6,'(1X,A,2I10)')'LIMT,LIM:',LIMT,LIM
      END IF
      IF(LIM.GT.LIC) THEN
        ISTOP=1
        WRITE(6,*)'ALSO: Too much storage needed for ABCI.'
        WRITE(6,'(1X,A,2I10)')'LIM,LIC:',LIM,LIC
      END IF
!PAM96411   WRITE(6,410)LIM
!PAM96410   FORMAT(6X,'STORAGE FOR SORTING ABCI',I7)
      IF(NOVT.GE.MADR) THEN
        ISTOP=2
        WRITE(6,*)'ALSO: Too much storage needed.'
        WRITE(6,'(1X,A,2I10)')'NOVT,MADR:',NOVT,MADR
      END IF
      RETURN
      END
