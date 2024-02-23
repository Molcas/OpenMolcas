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
      SUBROUTINE SGMONE(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM)
      use Struct, only: SGStruct, CIStruct, EXStruct
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CI(*),SGM(*)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS
#include "WrkSpc.fh"

      nMidV =CIS%nMidV

      MxEO  =EXS%MxEO
      nICoup=EXS%nICoup
      nVTab =EXS%nVTab
      lVTab =EXS%lVTab
      CALL SIGMA_1(SGS,CIS,EXS,NMIDV,MXEO,NVTAB,NICOUP,SGS%ISM,
     &             IP,IQ,CPQ,ISYCI,CI,SGM,CIS%NOCSF,
     &             CIS%IOCSF,CIS%NOW,CIS%IOW,
     &             EXS%NOCP,EXS%IOCP,EXS%ICOUP,
     &             WORK(LVTAB),EXS%MVL,EXS%MVR)

      end SUBROUTINE SGMONE
