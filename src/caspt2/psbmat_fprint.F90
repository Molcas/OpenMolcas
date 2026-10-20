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
function PSBMAT_FPRINT(lg_M,NM)

use Index_Functions, only: nTri_Elem
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: GA_DDot
#endif
use fake_ga, only: GA_arrays
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: PSBMAT_FPRINT
integer(kind=iwp), intent(in) :: lg_M, NM
real(kind=wp), external :: DNRM2_

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  PSBMAT_FPRINT = sqrt(GA_DDOT(lg_M,lg_M))
else
#endif
  PSBMAT_FPRINT = DNRM2_(nTri_Elem(NM),GA_Arrays(lg_M)%A,1)
#ifdef _MOLCAS_MPP_
end if
#endif

end function PSBMAT_FPRINT
