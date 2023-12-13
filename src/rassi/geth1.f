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
* Copyright (C) 1989, Per Ake Malmqvist                                *
************************************************************************
*****************************************************************
*  PROGRAM RASSI        PER-AAKE MALMQVIST
*  SUBROUTINE GETH1     IBM-3090 RELEASE 89 01 30
*  READ THE ONE-ELECTRON HAMILTONIAN MATRIX ELEMENTS AND RETURN
*  IT AS HONEAO IN SYMMETRY-BLOCKED SQUARED FORMAT.
*  Also reads and adds reaction field contribution.
*  Also ERFNUC, reaction field contribution to nuclear repulsion.
*****************************************************************
      SUBROUTINE GETH1_RASSI(HONEAO)
      use OneDat, only: sNoNuc, sNoOri
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION HONEAO(NBSQ)
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "WrkSpc.fh"
      Character*8 OneLbl
      Logical Found
*
      CALL GETMEM('H1    ','ALLO','REAL',LH1,NBTRI)
      iRc=-1
      iOpt=ibset(ibset(0,sNoOri),sNoNuc)
      iCmp=1
      iSyLab=1
      OneLbl='OneHam  '
      Call RdOne(iRc,iOpt,OneLbl,iCmp,Work(LH1),iSyLab)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  GETH1 ***'
        WRITE(6,*)'   BARE NUCLEI HAMILTONIAN IS NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
      ERFNuc=0.0D0
      If ( RFpert ) then
         Call f_Inquire('RUNOLD',Found)
         If (Found) Call NameRun('RUNOLD')
         Call GetMem('RFFLD','Allo','Real',ipTmp,nBtri)
         Call Get_dScalar('RF Self Energy',ERFNuc)
         Call Get_dArray('Reaction field',Work(ipTmp),nBtri)
         If (Found) Call NameRun('#Pop')
         Call Daxpy_(nBtri,1.0D0,Work(ipTmp),1,WORK(LH1),1)
         Call GetMem('RFFLD','Free','Real',ipTmp,nBtri)
      End If
      IBUF=LH1
      ISTQ=0
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        if (NB == 0) cycle
        DO IP=1,NB
          DO IQ=1,IP
            IPQ=NB*(IP-1)+IQ+ISTQ
            IQP=NB*(IQ-1)+IP+ISTQ
            HONEAO(IPQ)=WORK(IBUF)
            HONEAO(IQP)=WORK(IBUF)
            IBUF=IBUF+1
          END DO
        END DO
       ISTQ=ISTQ+NB**2
      END DO
      CALL GETMEM('      ','FREE','REAL',LH1,NBTRI)
      RETURN
      END
