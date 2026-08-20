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

use Molcas, only: MxLev
use fciqmc, only: DoNECI
use fcidump, only: DumpOnly
use CC_CI_mod, only: Do_CC_CI
use rasscf_global, only: NSM
use general_data, only: iSpin, nActel, nConf, nElec3, nHole1, nRs1, nRs2, nRs3, nSym, STSYM, NLEV, Level, &
                        iDoGAS, NGAS, NGSSH
use sguga, only: CIS, SGS, SG_init
#ifdef _DMRG_
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use general_data, only: nRas,nRasEl,nRsPrt
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(inout):: SkipGUGA
integer(kind=iwp) :: IGAS, iq, ISYM, nRs1T, NSTA
real(kind=wp) :: dum1, dum2, dum3, Eterna_1, Eterna_2
integer(kind=iwp), parameter :: istate=1

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do

Level(1:MxLev)=[(iq,iq=1,MxLev)]

if (nHole1+nElec3 /= 0) then
   nRsPrt=3
   nRas(:,1)=nRs1(:)
   nRas(:,2)=nRs2(:)
   nRas(:,3)=nRs3(:)
   nRs1T = sum(nRs1(1:nSym))
   nRasEl(1)=2*nRs1T-nHole1
   nRasEl(2)=nActel-nElec3
   nRasEl(3)=nActel
else
   nRsPrt=1
   nRas(:,1)=nRs2(:)
   nRasEl(1)=nActel
end if

! Construct the Guga tables

if (.not. (DoNECI .or. Do_CC_CI .or. DumpOnly .or. SkipGUGA .or. iDoGAS)) then

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

end if

end subroutine sg_setup_rasscf
