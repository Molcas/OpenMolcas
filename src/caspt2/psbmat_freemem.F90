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
      SUBROUTINE PSBMAT_FREEMEM(lg_M)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_ga, only: Deallocate_GA_Array
      use definitions, only: iwp
      IMPLICIT NONE
!SVC2010: destroy square global array S/B for symmetry iSYM
! with integer handle lg_M or if replicate or serial, free the
! tridiagonal local array at Work(lg_M)
      Integer(kind=iwp) lg_M

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      LOGICAL(kind=iwp) bStat
#endif


#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        bStat = GA_Destroy(lg_M)
      ELSE
#endif
        Call Deallocate_GA_Array(lg_M)
#ifdef _MOLCAS_MPP_
      END IF
#include "macros.fh"
      unused_var(bStat)
#endif
      END SUBROUTINE PSBMAT_FREEMEM
