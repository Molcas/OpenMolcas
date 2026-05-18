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
!
! WRAPPERS FOR PARALLEL S AND B MATRIX ROUTINES
!
      FUNCTION PSBMAT_FPRINT(lg_M,NM)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: GA_arrays
      use definitions, only: iwp, wp
      IMPLICIT NONE
      real(kind=wp) PSBMAT_FPRINT
      INTEGER(kind=iwp) lg_M, NM

      INTEGER(kind=iwp) nTri
      REAL(kind=wp), External::DNRM2_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        PSBMAT_FPRINT=SQRT(GA_DDOT(lg_M,lg_M))
      ELSE
#endif
        nTri=(NM*(NM+1))/2
        PSBMAT_FPRINT=DNRM2_(nTri,GA_Arrays(lg_M)%A(:),1)
#ifdef _MOLCAS_MPP_
      END IF
#endif
      END FUNCTION PSBMAT_FPRINT
