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

use sguga, only: SG_Free
use Str_Info, only: CFTP, CNSM
use input_mclr, only: nConf, nCSF, nSym, State_Sym, CIS, EXS, SGS
use stdalloc, only: mma_allocate, mma_deallocate

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lCI, imode, pState_Sym
real(kind=wp), intent(inout) :: CI(lCI)

integer(kind=iwp) :: iss
real(kind=wp), allocatable :: CINEW(:)
#ifdef _DEBUGPRINT_
real(kind=wp), parameter :: PRWTHR = 0.05_wp
#endif

! Transformation of CI vector to symmetric group from GUGA pepresentation, or the reverse

call SG_Setup_MCLR(pState_Sym)

NCSF(1:nSym) = CIS%NCSF(1:nSym)
NCONF        = CIS%NCSF(pState_Sym)

iss = 1
if (pState_sym /= state_sym) iss = 2

#ifdef _DEBUGPRINT_
write(u6,101)
write(u6,102) PRWTHR
call SG_PrWF(SGS,CIS,pState_sym,PRWTHR,SGS%iSpin,CI,nConf,.false.,-99)
write(u6,103)
101 format(/,6X,100('-'),/,6X,29X,'Wave function printout: Split Graph format',/,6X,8X, &
           'in parenthesis: midvertex, upper-walk symmetry upper- and lower-walk serial numbers',/,6X,100('-'),/)
102 format(6X,'printout of CI-coefficients larger than',F6.2)
103 format(/,6X,100('-'),/)
#endif

Call mma_allocate(CINEW,nConf,Label='CINEW')
call REORD(SGS,EXS,iMode,CNSM(iss)%ICONF,CFTP,pState_Sym,nConf,CI,CINEW)
CI(1:nConf)=CINEW(1:nConf)
Call mma_deallocate(CINEW)

#ifdef _DEBUGPRINT_
call SG_PrWF(SGS,CIS,pState_sym,PRWTHR,SGS%iSpin,CI,nConf,.false.,-99)
#endif

call SG_Free(SGS,CIS,EXS)

end subroutine sg2symg
