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

subroutine RHS_SCATTER(LDW,lg_W,Buff,idxW,nBuff)
!SVC: this routine scatters + adds values of a buffer array into the RHS
!     array at positions given by the buffer index array.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
#endif
use fake_GA, only: GA_Arrays
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LDW, lg_W, nBuff, idxW(nBuff)
real(kind=wp), intent(in) :: Buff(nBuff)
integer(kind=iwp) :: I
#ifdef _MOLCAS_MPP_
integer(kind=iwp), allocatable :: TMPW1(:), TMPW2(:)
#include "global.fh"
#include "mafdecls.fh"
#else
#include "macros.fh"
unused_var(LDW)
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  !SVC: global array RHS matrix expects 2 index buffers
  call mma_allocate(TMPW1,nBuff,Label='TMPW1')
  call mma_allocate(TMPW2,nBuff,Label='TMPW2')
  do I=1,nBuff
    TMPW2(I) = (idxW(I)-1)/LDW+1
    TMPW1(I) = idxW(I)-LDW*(TMPW2(I)-1)
  end do
  call GA_Scatter_Acc(lg_W,Buff,TMPW1,TMPW2,nBuff,One)
  call mma_deallocate(TMPW1)
  call mma_deallocate(TMPW2)
else
#endif
  do I=1,nBuff
    GA_Arrays(lg_W)%A(idxW(I)) = GA_Arrays(lg_W)%A(idxW(I))+BUFF(I)
  end do
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_SCATTER
