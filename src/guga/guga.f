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
************************************************************************
************************************************************************
*                                                                      *
* PER SIEGBAHN                                                         *
* DEPARTMENT OF THEORETICAL PHYSICS                                    *
* UNIVERSITY OF STOCKHOLM                                              *
* SWEDEN                                                               *
*                                                                      *
* UPDATED FOR MOLCAS-4 BY P-A MALMQVIST & NIGEL MORIARTY 1996          *
************************************************************************
      SUBROUTINE GUGA(IRETURN)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
      COMMON/D/JNDX(500 000)
*      DIMENSION JSYM(30000),SO(1 000 000),JSY(3000)
      DIMENSION JSYM(30000),JSY(3000)
      DIMENSION L0(4*MXVERT),L1(4*MXVERT),L2(4*MXVERT),L3(4*MXVERT)
*
*     Prologue
*
* Allocate workspace through GETMEM:
* PAM Aug -06: Get rid of fixed upper limit of workspace
* PAM      NCOR=1 000 000
* PAM      Call GETMEM('SOArr','Allo','Real',LSOArr,NCOR)
* Replace by: Find max possible allocatable
      Call GETMEM('SOArr','Max','Inte',LDummy,NCOR)
* Grab almost all of it, but leave a little to be safe:
      NCOR=NCOR-100000
      Call GETMEM('SOArr','Allo','Inte',LSOArr,NCOR)
*     Call SetQue('Trace=On')
      Call qEnter('GUGA')
      CALL SETTIM
      CALL XUFLOW
      CALL ERRSET(208,256,-1,1,1,208)
      CALL ERRSET(151,256,-1,1,1,151)
      CALL JTIME(IST)
*
*     Print program header
*
*
*     Initialize files
*
      Lu_11=11
      CALL DANAME_wa(Lu_11,'TEMP01')
      Lu_10=10
      CALL DANAME(Lu_10,'CIGUGA')
      DO 5 I=1,9
         IAD10(I)=0
5     CONTINUE
      IADD10=0
      CALL iDAFILE(Lu_10,1,IAD10,9,IADD10)
*
*     Read input
*
      IO=5
      IW=6
      ISPA=NCOR
      LIX=500 000
      NBUF=600
CPAM96: Use variable MCOP, size of buffers:
      MCOP=NBUF*RTOI+NBUF+1
      D0=0.0d0
      D1=1.0d0
      D2=2.0d0
      CALL INPUT_GUGA(iWork(LSOArr),JSYM,JSY,L0,L1,L2,L3,ISPAC)
*
*     Main body
*
C     SORT ALLOCATION , ISPAC WORDS TO BE SORTED IN NCOR CORE SPACE
C     TWO BUFFERS OF LENGTH NBINS NEEDED
C     NCOR IS IN UNITS OF FLOATING-POINT WORDS, e.g. REAL*8
C     NBINS=ISPAC/(NCORX-2*NBINS)+1
      NCORX=NCOR-(MCOP+1)
      A=D2
      B=-NCORX-2
      C=ISPAC+NCORX
      NBINS=INT((-B-SQRT(B*B-D2*D2*A*C))/(D2*A))
C     NUMBER OF WORDS IN EACH BIN
      NTPB=(ISPAC-1)/NBINS+1
C     SPACE IN CORE FOR EACH BIN
      KB=RTOI*(NCORX-2*NBINS)/NBINS
      KBUF=(KB-1)/(RTOI+1)
      KBUF=(KBUF/2)*2
      IF(KBUF.GT.600)KBUF=600
      KBUF2=KBUF*RTOI+KBUF+2
      IF(IPRINT.GE.2) WRITE(IW,10)KBUF,NBINS,NTPB,NCOR,ISPAC
10    FORMAT(/6X,'SORTING INFORMATION',
     */6X,'KBUF=',I7,/6X,'NBINS=',I6,/6X,'NTPB=',I7,
     */6X,'NCOR=',I7,/6X,'ISPAC=',I6)
C     STORAGE FOR NBINS BINS EACH OF SIZE KBUF2 IN AIAI
      LSTO=NBINS*KBUF2
C     ALSO SPACE FOR NTPB WORDS IN EMPTY
      IF(NTPB.GT.LSTO)LSTO=NTPB
      LW1=LSTO+1
      LW2=LW1+NBINS
      LIM=LW2+NBINS
      If (LIM.GT.NCOR) Then
         Write (6,*) 'Guga: LIM.GT.NCOR, position 1'
         Write (6,*) 'LIM,NCOR=',LIM,NCOR
         Call QTrace
         Call Abend
      End If
      LIM=LW2+KBUF2
      If (LIM.GT.NCOR) Then
         Write (6,*) 'Guga: LIM.GT.NCOR, position 2'
         Write (6,*) 'LIM,NCOR=',LIM,NCOR
         Call QTrace
         Call Abend
      End If
      CALL ICOPY(LW1,0,0,iWork(LSOArr),1)
      CALL CI_SELECT(iWork(LSOArr),iWork(LSOArr-1+LW1),
     &               iWork(LSOArr-1+LW2),
     &               L0,L1,L2,L3,KBUF,NTPB,NBINS)
      IADD10=0
      CALL iDAFILE(Lu_10,1,IAD10,9,IADD10)
*
*     Epilogue, end
      Call GETMEM('SOArr','Free','Inte',LSOArr,NCOR)
*
*                                                                      *
************************************************************************
*                                                                      *
*     Close open dafiles
*
      Call DaClos(Lu_10)
      Call DaClos(Lu_11)

*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('GUGA')
      ireturn=0
      Return
      End
