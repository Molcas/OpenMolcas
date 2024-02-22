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
      SUBROUTINE PRWF(SGS,CIS,ISYCI,CI,CITHR)
      use Struct, only: SGStruct, CIStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CI(*)
      Type (SGStruct) SGS
      Type (CIStruct) CIS
      Integer, Allocatable:: ICS(:)
#include "WrkSpc.fh"

      NLEV  =SGS%nLev

      NMIDV =CIS%nMidV
      LNOCSF=CIS%lNOCSF
      LIOCSF=CIS%lIOCSF

      CALL mma_allocate(ICS,NLEV,Label='ICS')
      CALL PRWF1(SGS,CIS,NLEV,NMIDV,
     &           SGS%ISM,ICS,IWORK(LNOCSF),
     &           IWORK(LIOCSF),CIS%NOW,CIS%IOW,
     &           ISYCI,CI,CITHR)
      CALL mma_deallocate(ICS)

      END SUBROUTINE PRWF
