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

function GTH1ES_MCLR(IORB,ITP,ISM,JORB,JTP,JSM)
! This is a fantastic solution of the problem hehehehe!!
!  EAW
!
! one electron integral between orbitals (iorb,itp,ism,jorb,jsm,jtp)
!
! correct combination of row and column symmetry is assumed

use MCLR_Data, only: IBsO, IBTSOB, IREOTS, KAIN1, pInt1
use input_mclr, only: nIsh, nOrb
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GTH1ES_MCLR
integer(kind=iwp), intent(in) :: IORB, ITP, ISM, JORB, JTP, JSM
integer(kind=iwp) :: I1, I_ABS, IJ, IREO, J1, J_ABS, JREO

I_ABS = IORB+IBTSOB(ITP,ISM)-1
IREO = IREOTS(I_ABS)
J_ABS = JORB+IBTSOB(JTP,JSM)-1
JREO = IREOTS(J_ABS)
I1 = IREO-IBSO(ISM)+1+nISH(ISM)
J1 = JREO-IBSO(JSM)+1+nISH(JSM)
IJ = pInt1(ISM)-1+(J1-1)*NORB(ISM)+I1
GTH1ES_MCLR = KAIN1(IJ)

end function GTH1ES_MCLR
