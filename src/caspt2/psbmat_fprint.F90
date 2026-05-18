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

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_ga, only: GA_arrays
use definitions, only: iwp, wp

implicit none
real(kind=wp) PSBMAT_FPRINT
integer(kind=iwp) lg_M, NM
integer(kind=iwp) nTri
real(kind=wp), external :: DNRM2_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  PSBMAT_FPRINT = sqrt(GA_DDOT(lg_M,lg_M))
else
#endif
  nTri = (NM*(NM+1))/2
  PSBMAT_FPRINT = DNRM2_(nTri,GA_Arrays(lg_M)%A(:),1)
#ifdef _MOLCAS_MPP_
end if
#endif

end function PSBMAT_FPRINT
