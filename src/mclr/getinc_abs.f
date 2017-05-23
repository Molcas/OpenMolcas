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
      SUBROUTINE GETINC_ABS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
     &                  IXCHNG,IKSM,JLSM,INTLST,IJKLOF,NSMOB,
     &                  ICOUL )
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
      iOrb=NTSOB(ITP,ISM)
      jOrb=NTSOB(JTP,JSM)
      kOrb=NTSOB(KTP,KSM)
      lOrb=NTSOB(LTP,LSM)
      iOff=IBTSOB(ITP,ISM)
      jOff=IBTSOB(JTP,JSM)
      kOff=IBTSOB(KTP,KSM)
      lOff=IBTSOB(LTP,LSM)
*
*     Collect Coulomb terms
*
      IF(ICOUL.EQ.0) THEN
      iint=1
      Do lBas=lOff,lOff+lOrb-1
        Do jBas=jOff,jOff+jOrb-1
          Do kBas=kOff,kOff+kOrb-1
            Do iBas=iOff,iOff+iOrb-1
              jINT = (lBas-1)*nACOB**3
     &              +(kBas-1)*nACOB**2
     &              +(jBas-1)*nACOB
     &              + iBas
              Xint(iInt) = Intlst(jint)
              iInt=iInt+1
            End Do
          End Do
        End Do
      End Do
*
*     Collect Exchange terms
*
      If ( IXCHNG.ne.0 ) Then
        iint=1
        Do lBas=lOff,lOff+lOrb-1
          Do jBas=joff,jOff+jOrb-1
            Do kBas=kOff,kOff+kOrb-1
              Do iBas=ioff,iOff+iOrb-1
              jINT = (jBas-1)*nACOB**3
     &              +(kBas-1)*nACOB**2
     &              +(lBas-1)*nACOB
     &              + iBas
                XInt(iInt)=XInt(iInt)-Intlst(jint)
                iInt=iInt+1
              End Do
            End Do
          End Do
        End Do
      End If
      ELSE IF(ICOUL.NE.0) THEN
        iint=0
        Do lBas=lOff,lOff+lOrb-1
          Do kBas=kOff,kOff+kOrb-1
           Do jBas=joff,jOff+jOrb-1
            Do iBas=ioff,iOff+iOrb-1
              JINT = (LBAS-1)*nACOB**3
     &              + (KBAS-1)*nACOB**2
     &              + (JBAS-1)*nACOB
     &              +  IBAS
              Xint(iInt) = Intlst(jint)
              iInt=iint+1
            End Do
          End Do
        End Do
      End Do
      END IF
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(IKSM)
        Call Unused_integer(JLSM)
        Call Unused_integer(IJKLOF)
      End If
      End
