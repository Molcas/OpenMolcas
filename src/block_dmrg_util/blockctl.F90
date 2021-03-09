!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2014, Naoki Nakatani                                   *
!***********************************************************************

subroutine BlockCtl(LW1,TUVX,IFINAL,IRST)
!***********************************************************************
!                                                                      *
!     Description                                                      *
!                                                                      *
!     This is an alternative of DavCtl for DMRG-CASSCF in DMRGCtl      *
!                                                                      *
!***********************************************************************
!                                                                      *
!     DMRG control section                                             *
!                                                                      *
!     calling arguments:                                               *
!     LW1     : array of real*8                                        *
!               Memory pointer to active Fock matrix                   *
!     TUVX    : array of real*8                                        *
!               two-electron integrals (tu!vx)                         *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!     IRST    : integer                                                *
!               restart flag of DMRG                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     N. Nakatani, Hokkaido University, Japan, 2014                    *
!                                                                      *
!***********************************************************************

use rasscf_data, only: Do3RDM, ENER, HFOCC, iOrbTyp, ITER, lRoots, mxSym, MxDMRG, NAC, ROTMAX, THRE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Five, Ten
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: LW1(*), TUVX(*)
integer(kind=iwp), intent(in) :: IFINAL, IRST
#include "general.fh"
integer(kind=iwp) :: iChMolpro(8), iOper, iOrb, iSigma, iSym, jOrb, JRST, lSymMolpro, nIrrep, NRDM_ORDER
real(kind=wp) :: ThDMRG, ThNoise
character(len=3) :: Label
integer(kind=iwp), allocatable :: OrbSym(:)

! Load symmetry info from RunFile
call Get_iScalar('NSYM',nIrrep)
call Get_iArray('Symmetry operations',iOper,nIrrep)
call Get_iScalar('Rotational Symmetry Number',iSigma)

! Get character table to convert MOLPRO symmetry format
call MOLPRO_ChTab(nSym,Label,iChMolpro)

! Convert orbital symmetry into MOLPRO format
call mma_allocate(OrbSym,NAC,label='OrbSym')
iOrb = 1
do iSym=1,nSym
  do jOrb=1,NASH(iSym)
    OrbSym(iOrb) = iChMolpro(iSym)
    iOrb = iOrb+1
  end do
end do
lSymMolpro = iChMolpro(stSym)

NRDM_ORDER = 2
if (NACTEL == 1) NRDM_ORDER = 1

! Default setting for the first iteraction
JRST = IRST
ThDMRG = THRE*Ten
ThNoise = 1.0e-4_wp

if (JRST /= 0) then
  if (IFINAL == 0) then
    if (abs(ROTMAX) >= 1.0e-1_wp) then
      ! No restart
      JRST = 0
      ThDMRG = THRE*Ten
      ThNoise = 1.0e-4_wp
    else if (abs(ROTMAX) >= 0.5e-1_wp) then
      ! Full-restart with noise
      JRST = 2
      ThDMRG = THRE*Five
      ThNoise = 1.0e-5_wp
    else if (abs(ROTMAX) >= 1.0e-2_wp) then
      ! Full-restart with smaller noise
      JRST = 2
      ThDMRG = THRE
      ThNoise = 1.0e-6_wp
    else
      ! Restart without noise
      JRST = 1
      ThDMRG = THRE
      ThNoise = Zero
    end if
  else if (IFINAL == 1) then
    ! Full-restart for CI-only
    JRST = 2
    ThDMRG = THRE
    ThNoise = 1.0e-6_wp
  else
    ! Full-restart for the final wavefunction
    if (iOrbTyp == 2) then
      ! with OutOrb = Canonical, this will be problematic for LMO-DMRG calc.
      JRST = 2
      ThDMRG = THRE
      ThNoise = 1.0e-6_wp
    else
      ! by default, just restarting
      JRST = 1
      ThDMRG = THRE
      ThNoise = Zero
    end if
  end if
end if

! In the following calls the extra array HFOCC has been passed
! from Molcas to Block to have a user customized reference det.
! This change must be done accordingly into the block code.
! In particular, in file molcas/block_calldmrg.C:
! line~23:  const char* hf_occ)
! line~24:  block_calldmrg(*Restart, *N_roots, *N_act, *N_elec, *M_s, Sym, *iSym, OrbSym, *E_core, h0, tuvx, *M_state, *N_pdm, *T_sweep, *T_noise, E_sweep, hf_occ);
! line~63:  const char* hf_occ)
! line~186: fcon << "hf_occ " << hf_occ << endl;
!
! In molcas/block_calldmrg.h:
! line~25:  const char* hf_occ);
! line~46:  const char* hf_occ);

! Compute DMRG
call block_calldmrg(JRST,lRoots,NAC,NACTEL,ISPIN-1,Label,lSymMolpro,OrbSym,Zero,LW1,TUVX,MxDMRG,NRDM_ORDER,ThDMRG, &
                    ThNoise,ENER(1,ITER),HFOCC,NRS2T)

if (IFINAL == 2 .and. Do3RDM .and. NACTEL > 2) then
  ! Compute 3RDM for DMRG-cu4-CASPT2
  call block_calldmrg(1,lRoots,NAC,NACTEL,ISPIN-1,Label,lSymMolpro,OrbSym,Zero,LW1,TUVX,MxDMRG,3,THRE,Zero,ENER(1,ITER), &
                      HFOCC,NRS2T)
end if

call mma_deallocate(OrbSym)

return

end subroutine BlockCtl
