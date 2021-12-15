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
      SUBROUTINE SGMONE(ISGSTRUCT,ICISTRUCT,IXSTRUCT,
     &                  IP,IQ,CPQ,ISYCI,CI,SGM)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SGMONE')
      DIMENSION CI(*),SGM(*)
#include "Struct.fh"
      Dimension iSGStruct(nSGSize)
      Dimension iCIStruct(nCISize)
      Dimension iXStruct (nXSize)
#include "WrkSpc.fh"




      nSym  =iSGStruct(1)
      lISm  =iSGStruct(3)
      lNOCSF=iCIStruct(6)
      lIOCSF=iCIStruct(7)
      lNOW  =iCIStruct(3)
      lIOW  =iCIStruct(4)
      lNOCP =iXStruct(2)
      lIOCP =iXStruct(3)
      lICoup=iXStruct(5)
      lVTab =iXStruct(7)
      lMVL  =iXStruct(8)
      lMVR  =iXStruct(9)
      nMidV =iCIStruct(1)
      MxEO  =iXStruct(1)
      nVTab =iXStruct(6)
      nICoup=iXStruct(4)
      CALL SIGMA_1(ISGSTRUCT,ICISTRUCT,IXSTRUCT,
     &             NMIDV,MXEO,NVTAB,NICOUP,IWORK(LISM),
     &             IP,IQ,CPQ,ISYCI,CI,SGM,IWORK(LNOCSF),
     &             IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

      return
      end
