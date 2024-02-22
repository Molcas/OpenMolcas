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
      SUBROUTINE SGMONE(SGS,CIS,IXSTRUCT,
     &                  IP,IQ,CPQ,ISYCI,CI,SGM)
      use Struct, only: nXSize, SGStruct, CIStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CI(*),SGM(*)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Dimension iXStruct (nXSize)
#include "WrkSpc.fh"

      nMidV =CIS%nMidV
      lNOCSF=CIS%lNOCSF
      lIOCSF=CIS%lIOCSF

      MxEO  =iXStruct(1)
      lNOCP =iXStruct(2)
      lIOCP =iXStruct(3)
      nICoup=iXStruct(4)
      lICoup=iXStruct(5)
      nVTab =iXStruct(6)
      lVTab =iXStruct(7)
      lMVL  =iXStruct(8)
      lMVR  =iXStruct(9)
      CALL SIGMA_1(SGS,CIS,IXSTRUCT,
     &             NMIDV,MXEO,NVTAB,NICOUP,SGS%ISM,
     &             IP,IQ,CPQ,ISYCI,CI,SGM,IWORK(LNOCSF),
     &             IWORK(LIOCSF),CIS%NOW,CIS%IOW,
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

      end SUBROUTINE SGMONE
