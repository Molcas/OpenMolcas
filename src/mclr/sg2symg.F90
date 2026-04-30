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

use sguga, only: CIS, EXS, SGS, SG_Free
use Str_Info, only: CFTP, CNSM
use input_mclr, only: iSpin, nActEl, nConf, nCSF, nElec3, nHole1, nRS1, nRS2, nRS3, nSym, State_Sym
use Definitions, only: iwp, wp

integer(kind=iwp), intent(in) :: lCI, imode, pState_Sym
real(kind=wp), intent(inout) :: CI(lCI)
integer(kind=iwp) iss

! Transformation of CI vector to symmetric group from GUGA pepresentation, or the reverse

call SG_Init_MCLR(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,pState_Sym)

NCSF(1:nSym) = CIS%NCSF(1:nSym)
NCONF        = CIS%NCSF(pState_Sym)

iss = 1
if (pState_sym /= state_sym) iss = 2

#ifdef _DEBUGPRINT_
Block
use Definitions, only: u6
real(kind=wp), parameter :: PRWTHR = 0.05_wp
write(u6,101)
101 format(/,6X,100('-'),/,6X,29X,'Wave function printout: Split Graph format',/, &
           6X,8X,'in parenthesis: midvertex, upper-walk symmetry upper- and lower-walk serial numbers',/,6X,100('-'),/)
write(u6,102) PRWTHR
102 format(6X,'printout of CI-coefficients larger than',F6.2)
call SG_PrWF(SGS,CIS,pState_sym,PRWTHR,SGS%iSpin,CI,nConf,.false.,-99)
write(u6,103)
103 format(/,6X,100('-'),/)
End Block
#endif

call REORD(SGS,CIS,EXS,NCONF,iMode,CNSM(iss)%ICONF,CFTP,pState_Sym,CI)

#ifdef _DEBUGPRINT_
Block
real(kind=wp), parameter :: PRWTHR = 0.05_wp
call SG_PrWF(SGS,CIS,pState_sym,PRWTHR,SGS%iSpin,CI,nConf,.false.,-99)
ENd Block
#endif


call SG_Free(SGS,CIS,EXS)

end subroutine sg2symg
