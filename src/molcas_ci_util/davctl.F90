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
! Copyright (C) 1989, Bjorn O. Roos                                    *
!               1989, Per Ake Malmqvist                                *
!               1991, Jeppe Olsen                                      *
!               1991,1996, Markus P. Fuelscher                         *
!               2000, Thorstein Thorsteinsson                          *
!***********************************************************************

subroutine DavCtl(nFMO,FMO,nTUVX,TUVX,IFINAL)
!***********************************************************************
!                                                                      *
!     CI control section                                               *
!                                                                      *
!     calling arguments:                                               *
!     FMO     : active Fock matrix                                     *
!               array of real                                          *
!     TUVX    : array of real                                          *
!               two-electron integrals (tu|vx)                         *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     B.O. Roos and P.Aa. Malmqvist                                    *
!     University of Lund, Sweden, 1989                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     - updated to use determinant based CI-procedures                 *
!       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
!     - updated for MOLCAS version 3                                   *
!       J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
!     - updated for integral direct and reaction field calculations    *
!       M.P. Fuelscher, University of Lund, Sweden, 1996               *
!     - various modifications                                          *
!       T. Thorsteinsson, University of Lund, Sweden, 2000             *
!                                                                      *
!***********************************************************************

use rasscf_global, only: Conv, Emy, Ener, hRoots, IADR15, ICICH, iCIOnly, Iter, ITERCI, kTight, lRoots, MAXJT, n_Keep, NAC, &
                         ThFact, ThrEn
use rasscf_files, only: JOBIPH, LUDAVID
use rasscf_global, only: NCONF, NDET, NSEL
use CI_Lib_Info, only: CI_Lib_Space
use CI_Input, only: HEXS_DEXS_Inp
use CI_Space_Def, only: Space_Removals
use Rank0_Sync, only: Rank0_Active, Rank0_Bcast
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: nFMO, nTUVX
real(kind=wp), intent(in) :: FMO(nFMO), TUVX(nTUVX)
integer(kind=iwp), intent(in) :: IFINAL
integer(kind=iwp) :: iDisk, ItLimit, jRoot, m_Sel, mSel, nMaxSel
real(kind=wp) :: ESize, Threshold, ThrRule
integer(kind=iwp), allocatable :: iSel(:)
real(kind=wp), allocatable :: CI_conv(:,:,:), CIVEC(:), ExplE(:), ExplV(:,:), Vec(:)

!-----------------------------------------------------------------------
! INITIALIZE THE DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------

lRoots = lRoots+hroots
call Ini_David(lRoots,nConf,nDet,nSel,n_keep,nAc,LuDavid)

!-----------------------------------------------------------------------
! COMPUTE THE DIAGONAL ELEMENTS OF THE HAMILTONIAN
!-------------------------------------------------------------------

! CIVEC: TEMPORARY CI VECTOR IN CSF BASIS

call mma_allocate(CIVEC,NCONF,label='CIVEC')
if (NAC > 0) call CIDIA(NCONF,CIVEC,LUDAVID)

!-----------------------------------------------------------------------
! OBTAIN STARTING VECTORS
!-----------------------------------------------------------------------

mSel = nSel
m_Sel = mSel
if (NAC == 0) m_Sel = 0
call mma_allocate(iSel,m_Sel,label='iSel')
call mma_allocate(ExplE,m_Sel,label='ExplE')
call mma_allocate(ExplV,m_Sel,mSel,label='ExplV')
! with HEXS/DEXS the explicit subspace holds at most the kept CSFs
nMaxSel = nConf
if (Space_Removals(HEXS_DEXS_Inp)) nMaxSel = CI_Lib_Space(iState)%nCSF_Kept

call CStart(CIVEC,FMO,TUVX,iSel,ExplE,ExplV,nMaxSel,IFINAL)

! Under MPI all processes start from rank 0's explicit Hamiltonian, H
! diagonal and start vectors: their own diagonalizations need not agree
! (degenerate states), nor the CI vectors on their own JOBIPH, which
! were transformed with their own orbitals (see Rank0_Sync)
call Rank0_Bcast(nSel)
if (m_Sel > 0) then
  call Rank0_Bcast(iSel)
  call Rank0_Bcast(ExplE)
  call Rank0_Bcast(ExplV)
end if
if ((NAC > 0) .and. Rank0_Active()) then
  call mma_allocate(Vec,nConf,label='Vec')
  call Load_H_diag(nConf,Vec,LuDavid)
  call Rank0_Bcast(Vec)
  call Save_H_diag(nConf,Vec,LuDavid)
  ! (CStart leaves them in the tmp slots if nSel = nMaxSel)
  do jRoot=1,lRoots
    if (nSel == nMaxSel) then
      call Load_tmp_CI_vec(jRoot,nConf,Vec,LuDavid)
      call Rank0_Bcast(Vec)
      call Save_tmp_CI_vec(jRoot,nConf,Vec,LuDavid)
    else
      call Load_CI_vec(jRoot,nConf,Vec,LuDavid)
      call Rank0_Bcast(Vec)
      call Save_CI_vec(jRoot,nConf,Vec,LuDavid)
    end if
  end do
  call mma_deallocate(Vec)
end if

!-----------------------------------------------------------------------
! DIAGONALIZATION SECTION
!-----------------------------------------------------------------------

! PAM Jun 2006: Gradual lowering of threshold in first 3 iterations
!if (Iter > 1) Threshold = THFACT*abs(CONV(4,ITER-1))
! Energy threshold for CI convergence criterion:
if (Iter == 1) then
  Threshold = THREN
else if ((ITER > 1) .and. (ITER <= 3)) then
  ThrRule = THFACT*abs(CONV(4,ITER-1))
  Threshold = (real(4-ITER,kind=wp)*THREN+real(ITER,kind=wp)*ThrRule)*Quart
else
  Threshold = THFACT*abs(CONV(4,ITER-1))
end if
! End of new rule, PAM Jun 2006
Threshold = max(Threshold,1.0e-9_wp)
if (NAC == 0) then
  ESize = abs(EMY)
else
  ESize = abs(ExplE(1))
end if
Threshold = max(Threshold,ESize*1.0e-14_wp)

! CI_conv: CONVERGENCE PARAMETERS
call mma_allocate(CI_conv,2,lRoots,MAXJT,label='CI_conv')
ITERCI = 1
if (NAC == 0) then
  ENER(1,ITER) = EMY
else
  if ((nSel == nConf) .or. (Space_Removals(HEXS_DEXS_Inp) .and. (nSel == CI_Lib_Space(iState)%nCSF_Kept))) then
    do jRoot=1,lRoots-hRoots
      ENER(jRoot,ITER) = ExplE(jRoot)
    end do
  else
    ! PAM Jun 2006: Limit the number of CI iterations in the
    ! first few macroiterations:
    if (KTIGHT == 0) ItLimit = min(12*ITER,MAXJT)
    if (KTIGHT == 1) ItLimit = MAXJT
    ! PAM Oct 2006: Full precision if this is final CI.
    if ((ICIONLY == 1) .or. (IFINAL == 2)) then
      Threshold = max(1.0e-9_wp,ESize*1.0e-14_wp)
      ITLIMIT = MAXJT
    end if
    ! PAM Feb 2009: New code in david5.

    call David5(nDet,ItLimit,IterCI,CI_conv,Threshold,iSel,ExplE,ExplV,nFMO,FMO,nTUVX,TUVX)

    do jRoot=1,lRoots-hRoots
      ENER(jRoot,ITER) = CI_conv(1,jRoot,ITERCI)
    end do
  end if
end if
call mma_deallocate(CI_conv)
call mma_deallocate(iSel)
call mma_deallocate(ExplE)
call mma_deallocate(ExplV)
nSel = mSel
lRoots = lRoots-hroots

!-----------------------------------------------------------------------
! CLEANUP AFTER THE DAVIDSON DIAGONALIZATION
!-----------------------------------------------------------------------

iDisk = IADR15(4)
call Term_David(ICICH,ITERCI,lRoots,nConf,CIVEC,JOBIPH,LuDavid,iDisk)
call mma_deallocate(CIVEC)

end subroutine DavCtl
