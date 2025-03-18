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

real*8 function GETH1I_MCLR(I,J)
! Obtain one -electron integral H(IORB,JOB)

use Arrays, only: FIMO
use MCLR_Data, only: ipCM
use MCLR_Data, only: IBSO, IREOTS, ISMFTO
use input_mclr, only: nIsh, nOrb

implicit none
integer I, J
integer ISM, JSM, I1, J1, IJ

ISM = ISMFTO(I)
JSM = ISMFTO(J)
I1 = IREOTS(I)-IBSO(ISM)+nISH(ISM)+1
J1 = IREOTS(J)-IBSO(JSM)+nISH(JSM)+1

IJ = ipCM(iSM)-1+(J1-1)*NORB(ISM)+I1

GETH1I_MCLR = FIMO(IJ)

end function GETH1I_MCLR
