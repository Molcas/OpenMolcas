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

function GETH1I_MCLR(I,J)
! Obtain one-electron integral H(IORB,JOB)

use MCLR_Lucia_Data, only: IBSO, IREOTS, ISMFTO, KAIN1, pInt1
use MCLR_Lucia_Data, only: nIsh => L_nIsh, nOrb => L_nOrb
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GETH1I_MCLR
integer(kind=iwp), intent(in) :: I, J
integer(kind=iwp) :: I1, IJ, ISM, J1, JSM

ISM = ISMFTO(I)
JSM = ISMFTO(J)
I1 = IREOTS(I)-IBSO(ISM)+nISH(ISM)+1
J1 = IREOTS(J)-IBSO(JSM)+nISH(JSM)+1

IJ = pInt1(iSM)-1+(J1-1)*NORB(ISM)+I1

GETH1I_MCLR = KAIN1(IJ)

end function GETH1I_MCLR
