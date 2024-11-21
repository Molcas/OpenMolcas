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
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: ERFNuc, RFPert
      use Symmetry_Info, only: nSym=>nIrrep
      IMPLICIT None
#include "rassi.fh"
      Real*8 HONEAO(NBSQ)
      Character(LEN=8) OneLbl
      Logical Found
      Real*8, Allocatable:: H1(:), Tmp(:)
      Integer IRC, IOPT, ICMP, iSyLab, iBuf, ISTQ, ISYM, NB, IP, IQ,
     &        IPQ, IQP
*
      CALL mma_allocate(H1,NBTRI,Label='H1')
      iRc=-1
      iOpt=ibset(ibset(0,sNoOri),sNoNuc)
      iCmp=1
      iSyLab=1
      OneLbl='OneHam  '
      Call RdOne(iRc,iOpt,OneLbl,iCmp,H1,iSyLab)
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
         Call mma_allocate(Tmp,nBtri,Label='Tmp')
         Call Get_dScalar('RF Self Energy',ERFNuc)
         Call Get_dArray('Reaction field',Tmp,nBtri)
         If (Found) Call NameRun('#Pop')
         Call Daxpy_(nBtri,1.0D0,Tmp,1,H1,1)
         Call mma_deallocate(Tmp)
      End If
      IBUF=1
      ISTQ=0
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        if (NB == 0) cycle
        DO IP=1,NB
          DO IQ=1,IP
            IPQ=NB*(IP-1)+IQ+ISTQ
            IQP=NB*(IQ-1)+IP+ISTQ
            HONEAO(IPQ)=H1(IBUF)
            HONEAO(IQP)=H1(IBUF)
            IBUF=IBUF+1
          END DO
        END DO
       ISTQ=ISTQ+NB**2
      END DO
      Call mma_deallocate(H1)

      END SUBROUTINE GETH1_RASSI
