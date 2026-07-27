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

<<<<<<< HEAD
use sguga, only: CIS, EXS, SGS, SG_Free
use Str_Info, only: CFTP, CNSM
use input_mclr, only: nConf, nCSF, nSym, State_Sym
use Definitions, only: iwp, wp

integer(kind=iwp), intent(in) :: lCI, imode, pState_Sym
real(kind=wp), intent(inout) :: CI(lCI)
integer(kind=iwp) iss
=======
use sguga, only: SG_Free
use Str_Info, only: CNSM, CFTP_MCLR=>CFTP
use lucia_data, only: CONF_OCC, CFTP
use input_mclr, only: nConf, nCSF, nSym, State_Sym
use sguga_states, only: CIS, EXS, SGS
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
integer(kind=iwp), Parameter:: istate=1
#ifdef _DEBUGPRINT_
real(kind=wp), parameter :: PRWTHR = 0.05_wp
#endif
>>>>>>> upstream-openmolcas/master

! Transformation of CI vector to symmetric group from GUGA pepresentation, or the reverse

call SG_Setup_MCLR(pState_Sym)

<<<<<<< HEAD
NCSF(1:nSym) = CIS%NCSF(1:nSym)
NCONF        = CIS%NCSF(pState_Sym)
=======
NCSF(1:nSym) = CIS(istate)%NCSF(1:nSym)
NCONF        = CIS(istate)%NCSF(pState_Sym)
>>>>>>> upstream-openmolcas/master

iss = 1
if (pState_sym /= state_sym) iss = 2

#ifdef _DEBUGPRINT_
<<<<<<< HEAD
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
=======
write(u6,101)
write(u6,102) PRWTHR
call SG_PrWF(SGS(istate),CIS(istate),pState_sym,PRWTHR,SGS(istate)%iSpin,CI,nConf,.false.,-99)
write(u6,103)
101 format(/,6X,100('-'),/,6X,29X,'Wave function printout: Split Graph format',/,6X,8X, &
           'in parenthesis: midvertex, upper-walk symmetry upper- and lower-walk serial numbers',/,6X,100('-'),/)
102 format(6X,'printout of CI-coefficients larger than',F6.2)
103 format(/,6X,100('-'),/)
#endif

Call mma_allocate(CINEW,nConf,Label='CINEW')
Call mma_allocate(Conf_Occ(pState_Sym)%A,SIZE(CNSM(iss)%ICONF),Label='CINEW')
Conf_Occ(pState_Sym)%A(:)=-CNSM(iss)%ICONF
Call mma_allocate(CFTP,SIZE(CFTP_MCLR),Label='CFTP')
CFTP(:)=CFTP_MCLR(:)

call SG_REORD(SGS(istate),EXS(istate),pState_Sym,iMode,nConf,CI,CINEW)
CI(1:nConf)=CINEW(1:nConf)

Call mma_deallocate(CFTP)
Call mma_deallocate(Conf_Occ(pState_Sym)%A)
Call mma_deallocate(CINEW)

#ifdef _DEBUGPRINT_
call SG_PrWF(SGS(istate),CIS(istate),pState_sym,PRWTHR,SGS(istate)%iSpin,CI,nConf,.false.,-99)
#endif

call SG_Free(SGS(istate),CIS(istate),EXS(istate))
>>>>>>> upstream-openmolcas/master

end subroutine sg2symg
