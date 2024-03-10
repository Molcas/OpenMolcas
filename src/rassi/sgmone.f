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
      use symmul_data, only: nSym
      IMPLICIT REAL*8 (A-H,O-Z)
      Real*8 CI(*),SGM(*)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Type (EXStruct) EXS
      Integer nMidV

      nMidV =CIS%nMidV

      CALL SIGMA1(SGS,CIS,EXS,
     &            IP,IQ,CPQ,ISYCI,CI,SGM,CIS%NOCSF,
     &            CIS%IOCSF,CIS%NOW,CIS%IOW,
     &            NMIDV,nSym)

      end SUBROUTINE SGMONE
