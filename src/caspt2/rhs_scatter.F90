!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

      SUBROUTINE RHS_SCATTER (LDW,lg_W,Buff,idxW,nBuff)
!SVC: this routine scatters + adds values of a buffer array into the RHS
!     array at positions given by the buffer index array.
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
      use stdalloc, only: mma_allocate, mma_deallocate
#endif
      use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
      use constants, only: One
#endif
      use definitions, only: iwp, wp
      IMPLICIT None
      Integer(kind=iwp), intent(in):: LDW,lg_W,nBuff
      real(kind=wp), intent(in):: Buff(nBuff)
      Integer(kind=iwp), intent(in):: idxW(nBuff)

      Integer(kind=iwp) I
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
      Integer(kind=iwp), allocatable:: TMPW1(:), TMPW2(:)
#else
#include "macros.fh"
      unused_var(LDW)
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
!SVC: global array RHS matrix expects 2 index buffers
        CALL mma_allocate(TMPW1,nBuff,Label='TMPW1')
        CALL mma_allocate(TMPW2,nBuff,Label='TMPW2')
        DO I=1,nBuff
          TMPW2(I)=(idxW(I)-1)/LDW+1
          TMPW1(I)=idxW(I)-LDW*(TMPW2(I)-1)
        END DO
        CALL GA_Scatter_Acc (lg_W,Buff,TMPW1,TMPW2,nBuff,One)
        CALL mma_deallocate(TMPW1)
        CALL mma_deallocate(TMPW2)
      ELSE
#endif
        DO I=1,nBuff
          GA_Arrays(lg_W)%A(idxW(I)) =                                  &
     &      GA_Arrays(lg_W)%A(idxW(I)) + BUFF(I)
        END DO
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE RHS_SCATTER
