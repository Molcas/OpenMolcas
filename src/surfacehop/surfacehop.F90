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
! Copyright (C) 2015, Luis Manuel Frutos                               *
!               2015, Ignacio Fdez. Galvan                             *
!               2015, Alessio Valentini                                *
!***********************************************************************

subroutine surfacehop(rc)

use Tully_variables
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: rc
#include "warnings.fh"

integer(kind=iwp) :: NSTATE, LUIPH, IAD, ITOC15(15), NCI, IDISK, I
real(kind=wp), allocatable :: CIBigArray(:)

call initial_surfacehop()
call rdinp_surfacehop()

LUIPH=20
CALL DANAME(LUIPH,"JOBIPH")
IAD=0
call IDAFILE(LUIPH,2,ITOC15,15,IAD)
call getIphInfo(LUIPH,NCI,NSTATE,ITOC15)
call mma_allocate(CIBigArray,NCI*NSTATE)
CIBigArray(:)=0.0_wp

IDISK=ITOC15(4)
do I=1,NSTATE
  call DDAFILE(LUIPH,2,CIBigArray(NCI*(I-1)+1),NCI,IDISK)
end do
call DACLOS(LUIPH)

!call recprt('CI coefficients','',CIBigArray,NCI,NSTATE)

call tully(CIBigArray,NSTATE,NCI)

call mma_deallocate(CIBigArray)
rc=_RC_ALL_IS_WELL_

return

end subroutine surfacehop
