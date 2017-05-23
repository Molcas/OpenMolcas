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
      SUBROUTINE GETINC_ABS_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IKSM,JLSM,INTLST,IJKLOF,NSMOB,
     &                  ICTL)
*
* Obtain integrals
* ICOUL = 0 :      XINT(IK,JL) = (IJ!KL) for IXCHNG = 0
*                              = (IJ!KL)-(IL!KJ) for IXCHNG = 1
* ICOUL = 1 :      XINT(IJ,KL) = (IJ!KL)
*
* Version for integrals stored in INTLST
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "detdim.fh"
#include "orbinp_mclr.fh"
      Real * 8 Intlst(*)
      Dimension IJKLof(NsmOB,NsmOb,NsmOB)
      DIMENSION XINT(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      iOrb=NTSOB(ITP,ISM)
      jOrb=NTSOB(JTP,JSM)
      kOrb=NTSOB(KTP,KSM)
      lOrb=NTSOB(LTP,LSM)
      iOff=IBTSOB(ITP,ISM)
      jOff=IBTSOB(JTP,JSM)
      kOff=IBTSOB(KTP,KSM)
      lOff=IBTSOB(LTP,LSM)
      ntash=nacob
*
*     Collect Coulomb terms
*
      IF(ICTL.EQ.1) THEN
      iint=1
      Do lBas=lOff,lOff+lOrb-1
        jMin=jOff
        If ( JLSM.ne.0 ) jMin=lBas
        Do jBas=jMin,jOff+jOrb-1
          Do kBas=kOff,kOff+kOrb-1
            iMin=iOff
            If(IKSM.ne.0) iMin = kBas
            Do iBas=iMin,iOff+iOrb-1
              ij = jbas+ntash*(ibas-1)
              kl = lbas+ntash*(kbas-1)
              ijkl=itri(ij,kl)
              Xint(iInt) = Intlst(ijkl)
              iInt=iInt+1
            End Do
          End Do
        End Do
      End Do
      Else If (ICTL.eq.4) Then
*
*     Collect Coulomb-Exchange terms
*
        iint=1
        Do lBas=lOff,lOff+lOrb-1
          jMin=jOff
          If ( JLSM.ne.0 ) jMin=lBas
          Do jBas=jMin,jOff+jOrb-1
            Do kBas=kOff,kOff+kOrb-1
              iMin=iOff
              If(IKSM.ne.0) iMin = kBas
              Do iBas=iMin,iOff+iOrb-1
                 il = ibas+ntash*(lbas-1)
                 jk = kbas+ntash*(jbas-1)
                 iljk=itri(il,jk)
                 ij = ibas+ntash*(jbas-1)
                 kl = kbas+ntash*(lbas-1)
                 ijkl=itri(ij,kl)
                 XInt(iInt)=Intlst(ijkl)-Intlst(iljk)
                 iInt=iInt+1
              End Do
            End Do
          End Do
        End Do
      Else
       Call Abend()
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(IJKLOF)
      End
