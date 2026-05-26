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

! WRAPPER FOR PARALLEL S AND B MATRIX ROUTINES
subroutine PSBMAT_GETMEM(cNAME,lg_M,nSize)
!SVC2010: create square global array S/B for symmetry iSYM
! with integer handle lg_M or if replicate or serial, create
! tridiagonal local array at Work(lg_M)

use Index_Functions, only: nTri_Elem
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_ga, only: Allocate_GA_Array, GA_arrays
use Constants, only: Zero
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: cNAME
integer(kind=iwp), intent(out) :: lg_M
integer(kind=iwp), intent(in) :: nSize

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_CREATE_STRIPED('H',nSize,nSize,cNAME,LG_M)
  call GA_ZERO(LG_M)
else
#endif
  lg_M = Allocate_GA_Array(nTri_Elem(nSize),cName)
  GA_Arrays(lg_M)%A(:) = Zero
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine PSBMAT_GETMEM
