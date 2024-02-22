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
      SUBROUTINE SGMONE(SGS,ICISTRUCT,CIS,IXSTRUCT,
     &                  IP,IQ,CPQ,ISYCI,CI,SGM)
      use Struct, only: nCISize, nXSize, SGStruct, CIStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CI(*),SGM(*)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Dimension iCIStruct(nCISize)
      Dimension iXStruct (nXSize)
#include "WrkSpc.fh"

      lNOCSF=iCIStruct(6)
      lIOCSF=iCIStruct(7)
      lNOW  =iCIStruct(3)
      lIOW  =iCIStruct(4)

      lNOCSF=CIS%lNOCSF
      lIOCSF=CIS%lIOCSF
      lNOW  =CIS%lNOW
      lIOW  =CIS%lIOW

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
      CALL SIGMA_1(SGS,ICISTRUCT,CIS,IXSTRUCT,
     &             NMIDV,MXEO,NVTAB,NICOUP,SGS%ISM,
     &             IP,IQ,CPQ,ISYCI,CI,SGM,IWORK(LNOCSF),
     &             IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

      end SUBROUTINE SGMONE
