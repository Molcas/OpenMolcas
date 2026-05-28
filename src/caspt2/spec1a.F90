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

subroutine SPEC1A(IFC,FACT,ISYM,X,nX,Y,nY)
! If IFC=0, compute
! X(tuu,i) <- X(tuu,i)+FACT*Y(t,i), else
! the conjugate expression (summing into Y, values from X).

use SUPERINDEX, only: KTUV
use caspt2_module, only: nAES, nAsh, nAshT, nIsh, nTUV, nTUVES
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IFC, ISYM, nX, nY
real(kind=wp), intent(in) :: FACT
real(kind=wp), intent(inout) :: X(nX), Y(nY)
integer(kind=iwp) :: IT, ITQ, ITUU, IUQ, NAS, NI, NT

NI = NISH(ISYM)
if (NI <= 0) return
NAS = NTUV(ISYM)
NT = NASH(ISYM)
do IT=1,NT
  ITQ = IT+NAES(ISYM)
  do IUQ=1,NASHT
    ITUU = KTUV(ITQ,IUQ,IUQ)-NTUVES(ISYM)
    if (IFC == 0) then
      call DAXPY_(NI,FACT,Y(IT),NT,X(ITUU),NAS)
    else
      call DAXPY_(NI,FACT,X(ITUU),NAS,Y(IT),NT)
    end if
  end do
end do

end subroutine SPEC1A
