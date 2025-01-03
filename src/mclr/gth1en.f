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
      REAL*8 FUNCTION GTH1EN(IORB,ITP,ISM,JORB,JTP,JSM)
      use Arrays, only: KAIN1, pInt1
      use MCLR_Data, only: IBsO,IBTSOB,IREOTS
*
* One-electron integral for active
* orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
*
      IMPLICIT None
      Integer IORB,ITP,ISM,JORB,JTP,JSM
      Real*8, External:: GTH1ES_MCLR

      GTH1EN = GTH1ES_MCLR(IREOTS(1),pInt1,KAIN1,IBSO,
     &                     IBTSOB,IORB,ITP,ISM,JORB,JTP,JSM)

      END FUNCTION GTH1EN
