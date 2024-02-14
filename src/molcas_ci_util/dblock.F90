!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE DBLOCK(D)
!
!     RASSCF program version IBM-3090: SX section
!
!     Purpose: To symmetry-block a full matrix d over the active orbital
!              the result is overlaid on the input matrix.
!
!     ********** IBM-3090 release 88 10 10 **********
!
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
      Real*8 D(*)

      Integer IA, NTU, ISYM, NA, NAT, ITU
!
      IA=NASH(1)
      NTU=ITRI(IA+1)
      DO ISYM=2,NSYM
       NA=NASH(ISYM)
       IF(NA.EQ.0) Cycle
       DO NAT=1,NA
        DO NAU=1,NAT
         NTU=NTU+1
         ITU=ITRI(NAT+IA)+NAU+IA
         D(NTU)=D(ITU)
        END DO
       END DO
       IA=IA+NASH(ISYM)
      END DO
!
      END SUBROUTINE DBLOCK
