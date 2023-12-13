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

subroutine DavCtl(LW1,TUVX,IFINAL)
!***********************************************************************
!                                                                      *
!     CI control section                                               *
!                                                                      *
!     calling arguments:                                               *
!     LW1     : active Fock matrix                                     *
!               array of real*8                                        *
!     TUVX    : array of real*8                                        *
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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Quart
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: LW1(*), TUVX(*)
integer(kind=iwp), intent(in) :: IFINAL
integer(kind=iwp) :: iDisk, ItLimit, jRoot, m_Sel, mSel, nMaxSel
real(kind=wp) :: ESize, Threshold, ThrRule
integer(kind=iwp), allocatable :: iSel(:)
real(kind=wp), allocatable :: CI_conv(:,:,:), CIVEC(:), ExplE(:), ExplV(:,:)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "lucia_ini.fh"

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
if (NAC > 0) call CIDIA_CI_UTIL(NCONF,STSYM,CIVEC,LUDAVID)

!-----------------------------------------------------------------------
! OBTAIN STARTING VECTORS
!-----------------------------------------------------------------------

mSel = nSel
m_Sel = mSel
if (NAC == 0) m_Sel = 0
call mma_allocate(iSel,m_Sel,label='iSel')
call mma_allocate(ExplE,m_Sel,label='ExplE')
call mma_allocate(ExplV,m_Sel,mSel,label='ExplV')
nMaxSel = nConf
if (N_ELIMINATED_GAS_MOLCAS > 0) nmaxSel = nCSF_HEXS

call CStart_CI_Util(CIVEC,LW1,TUVX,iSel,ExplE,ExplV,nMaxSel,IFINAL)

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
  if ((nSel == nConf) .or. ((N_ELIMINATED_GAS_MOLCAS > 0) .and. (nSel == nCSF_HEXS))) then
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
    !call David5(nAc,stSym,nDet,MAXJT,ITERCI,
    !call David5(nAc,stSym,nDet,ItLimit,ITERCI,CI_conv,Threshold,LW1,TUVX,iSel,ExplE,ExplV)

    call David5(nDet,ItLimit,IterCI,CI_conv,Threshold,iSel,ExplE,ExplV,LW1,TUVX)

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

return

end subroutine DavCtl
