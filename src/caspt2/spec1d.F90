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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine SPEC1D(IFC,FACT,X,nX,Y,nY)
! If IFC=0, compute
! X(tt1,ia) <- X(tt1,ia)+FACT*Y(ia), else
! the conjugate expression (summing into Y, values from X).

use SUPERINDEX, only: KTU
use caspt2_module, only: nAshT, nASup, nISup
use constants, only: Zero
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: IFC, nX, nY
real(kind=wp), intent(in) :: FACT
real(kind=wp), intent(inout) :: X(nX), Y(nY)
integer(kind=iwp) NAS, NIS, ITQ, ITT

NIS = NISUP(1,5)
if (NIS == 0) return
NAS = NASUP(1,5)
if (IFC == 0) then
  do ITQ=1,NASHT
    ITT = KTU(ITQ,ITQ)
    call DAXPY_(NIS,FACT,Y,1,X(ITT),NAS)
  end do
else
  Y(1:NIS) = Zero
  do ITQ=1,NASHT
    ITT = KTU(ITQ,ITQ)
    call DAXPY_(NIS,FACT,X(ITT),NAS,Y,1)
  end do
end if

end subroutine SPEC1D
