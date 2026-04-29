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

subroutine sg2symg(CI,lCI,imode,pState_Sym)

use sguga, only: CIS, EXS, SGS
use input_mclr, only: iSpin, nActEl, nConf, nCSF, nElec3, nHole1, nRS1, nRS2, nRS3, nSym, State_Sym
use Definitions, only: iwp, wp

integer(kind=iwp), intent(in) :: lCI, imode, pState_Sym
real(kind=wp), intent(inout) :: CI(lCI)

! Transformation of CI vector to symmetric group from GUGA pepresentation, or the reverse

call GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,CI,imode,pState_Sym,State_Sym)
NCSF(1:nSym) = CIS%NCSF(1:nSym)
NCONF        = CIS%NCSF(State_Sym)
call SGUGA_Free(SGS,CIS,EXS)

end subroutine sg2symg
