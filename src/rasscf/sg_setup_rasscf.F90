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

subroutine SG_Setup_RASSCF(SkipGUGA)

use general_data, only: iSpin, nActel, nConf, nSym, STSYM, NLEV, Level, &
                        nRas,nRasEl,nRsPrt, NSM
use sguga, only: CIS, SGS, SG_init
#ifdef _DMRG_
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in):: SkipGUGA
real(kind=wp) :: dum1, dum2, dum3, Eterna_1, Eterna_2
integer(kind=iwp), parameter :: istate=1

if (SkipGUGA) Return

! Construct the Guga tables

call Timing(Eterna_1,dum1,dum2,dum3)
call SG_Init(iState,nSym,nActEl,iSpin,                    &
             nRas,nRasEl,nRsPrt,                           &
             xLevel=Level,xL2Act=Level,xNLEV=NLEV,xNSM=NSM)

if (SGS(istate)%NVERT0 == 0) then
  CIS(istate)%NCSF(STSYM) = 0
else
  if (NActEl == 0) CIS(istate)%NCSF(STSYM) = 1
end if

call SETSXCI()
NCONF = CIS(istate)%NCSF(STSYM)

call Timing(Eterna_2,dum1,dum2,dum3)

end subroutine sg_setup_rasscf
