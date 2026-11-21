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

use CI_Info, only: CI_St
use CI_States, only: CI_ReOrd_Orig
use Str_Info, only: CFTP_MCLR => CFTP, CNSM
use lucia_data, only: CFTP, CONF_OCC, MINOP, NCNFTP, NCSFTP, NTYP
use MCLR_Lucia_Data, only: MINOP_MCLR => MINOP, NCNFTP_MCLR => NCNFTP, NCSFTP_MCLR => NCSFTP, NTYP_MCLR => NTYP
use Orbital_Spaces, only: nSym
use input_mclr, only: NCSF, NCONF
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use CI_States, only: CI_PrWF
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lCI, imode, pState_Sym
real(kind=wp), intent(inout) :: CI(lCI)
integer(kind=iwp) :: iss
real(kind=wp), allocatable :: CINEW(:)
integer(kind=iwp), parameter :: istate = 1
#ifdef _DEBUGPRINT_
real(kind=wp), parameter :: PRWTHR = 0.05_wp
#endif

! Transformation of CI vector to symmetric group from GUGA pepresentation, or the reverse

! (the CI state of the run, SG_Setup_MCLR in RdJobIph)

NCSF(1:nSym) = CI_St(istate)%nCSF(1:nSym)
NCONF = CI_St(istate)%nCSF(pState_Sym)

iss = 1
if (pState_sym /= CI_St(iState)%STSym) iss = 2

#ifdef _DEBUGPRINT_
write(u6,101)
write(u6,102) PRWTHR
call CI_PrWF(istate,pState_sym,PRWTHR,CI_St(istate)%iSpin,CI,nConf,.false.,-99)
write(u6,103)
101 format(/,6X,100('-'),/,6X,29X,'Wave function printout: Split Graph format',/,6X,8X, &
           'in parenthesis: midvertex, upper-walk symmetry upper- and lower-walk serial numbers',/,6X,100('-'),/)
102 format(6X,'printout of CI-coefficients larger than',F6.2)
103 format(/,6X,100('-'),/)
#endif

call mma_allocate(CINEW,nConf,Label='CINEW')
call mma_allocate(Conf_Occ(pState_Sym)%A,size(CNSM(iss)%ICONF),Label='CINEW')
Conf_Occ(pState_Sym)%A(:) = -CNSM(iss)%ICONF
call mma_allocate(CFTP,size(CFTP_MCLR),Label='CFTP')
CFTP(:) = CFTP_MCLR(:)
! MCLR's configuration types, for the reordering in the CI interface
! (sg_lucia_bridge reads Lucia's tables)
MINOP = MINOP_MCLR
NTYP = NTYP_MCLR
NCNFTP(:,:) = NCNFTP_MCLR(:,:)
NCSFTP(:) = NCSFTP_MCLR(:)

call CI_ReOrd_Orig(istate,pState_Sym,iMode,nConf,CI,CINEW)

CI(1:nConf) = CINEW(:)

call mma_deallocate(CFTP)
call mma_deallocate(Conf_Occ(pState_Sym)%A)
call mma_deallocate(CINEW)

#ifdef _DEBUGPRINT_
call CI_PrWF(iState,pState_sym,PRWTHR,CI_St(istate)%iSpin,CI,nConf,.false.,-99)
#endif

end subroutine sg2symg
