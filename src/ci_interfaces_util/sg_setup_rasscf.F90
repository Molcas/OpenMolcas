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

subroutine SG_Setup_RASSCF()

use general_data, only: iDoGAS, iSpin, Level, nActel, NLEV, nRas, nRasEl, nRsPrt, NSM, nSym, STSYM
use rasscf_global, only: nConf
use CI_States, only: CI_Setup_SG
use CI_Info, only: CI_St, SG_Rel_State, Space_CI
use Molcas, only: MxSym
use stdalloc, only: mma_allocate, mma_deallocate
#ifdef _DMRG_
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: dum1, dum2, dum3, Eterna_1, Eterna_2
integer(kind=iwp), parameter :: istate = 1
integer(kind=iwp) :: k, nP
integer(kind=iwp), allocatable :: nRasG(:,:), MinR(:), MaxR(:)

! Construct the Guga tables

call Timing(Eterna_1,dum1,dum2,dum3)

if (iDoGAS) then
  ! GAS: the partitions and the cumulative lower and upper bounds of the
  ! common CI space description
  call mma_allocate(nRasG,MxSym,Space_CI%nPart,Label='nRasG')
  nRasG(:,:) = 0
  nRasG(1:nSym,:) = Space_CI%nOrb(1:nSym,1:Space_CI%nPart)
  call CI_Setup_SG(iState,nSym,nActEl,iSpin,nRasG,Space_CI%MinEl,Space_CI%nPart,Level,Level,NLEV,NSM,'FULL', &
                   nRasMx=Space_CI%MaxEl)
  ! the relaxed space: bounds widened by one electron at every boundary
  ! (a single excitation changes the cumulative counts by at most one)
  nP = Space_CI%nPart
  call mma_allocate(MinR,nP,Label='MinR')
  call mma_allocate(MaxR,nP,Label='MaxR')
  do k=1,nP
    if (k < nP) then
      MinR(k) = max(0,Space_CI%MinEl(k)-1)
      MaxR(k) = min(nActEl,Space_CI%MaxEl(k)+1)
    else
      MinR(k) = Space_CI%MinEl(k)
      MaxR(k) = Space_CI%MaxEl(k)
    end if
  end do
  call CI_Setup_SG(SG_Rel_State,nSym,nActEl,iSpin,nRasG,MinR,nP,Level,Level,NLEV,NSM,'FULL',nRasMx=MaxR)
  CI_St(iState)%iRel = SG_Rel_State
  call mma_deallocate(MaxR)
  call mma_deallocate(MinR)
  call mma_deallocate(nRasG)
else
  call CI_Setup_SG(iState,nSym,nActEl,iSpin,nRas,nRasEl,nRsPrt,Level,Level,NLEV,NSM,'FULL')
end if

call SETSXCI()

NCONF = CI_St(istate)%nCSF(STSYM)

call Timing(Eterna_2,dum1,dum2,dum3)

end subroutine sg_setup_rasscf
