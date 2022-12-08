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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Tra_Ctl2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)
!***********************************************************************
!                                                                      *
!     main control section for:                                        *
!     - transformation of ERIs from AO to MO basis                     *
!     - Fock matrix generation                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*), D1I(*), D1A(*), ExFac
real(kind=wp), intent(inout) :: PUVX(*)
real(kind=wp), intent(_OUT_) :: TUVX(*)
real(kind=wp), intent(out) :: FI(*), FA(*)
integer(kind=iwp), intent(in) :: IPR
logical(kind=iwp), intent(in) :: lSquare
#include "rasdim.fh"
#include "general.fh"
integer(kind=iwp) :: iAsh, iBas, iDisk, iFro, iIsh, ij_Bas_pairs, ij_Orb_pairs, ijSym, iOff, iOrb, iStack, iSym, jAsh, jBas, jFro, &
                     jIsh, jOrb, jSym, kAsh, kBas, kFro, kIsh, kl_bas_pairs, kl_Orb_pairs, kOrb, kSym, kSymMax, lAsh, lBas, lFro, &
                     lIsh, lOrb, lSym, nPUVX, off_ltMat(mxSym), off_PUVX(mxSym,mxSym,mxSym), off_sqMat(mxSym)

if (IPR > 1) then
  write(u6,*)
  write(u6,*) ' Enter transformation section'
  write(u6,*) ' ============================'
  write(u6,*)
end if

! generate offsets

iStack = 0
do iSym=1,nSym
  off_sqMat(iSym) = iStack
  iBas = nBas(iSym)
  iStack = iStack+iBas*iBas
end do

iStack = 0
do iSym=1,nSym
  off_ltMat(iSym) = iStack
  iBas = nBas(iSym)
  iStack = iStack+(iBas*iBas+iBas)/2
end do

iStack = 0
do iSym=1,nSym
  iOrb = nOrb(iSym)
  do jSym=1,nSym
    jAsh = nAsh(jSym)
    ijSym = Mul(iSym,jSym)
    do kSym=1,nSym
      kAsh = nAsh(kSym)
      lSym = Mul(ijSym,kSym)
      if (lSym <= kSym) then
        lAsh = nAsh(lSym)
        kl_Orb_pairs = kAsh*lAsh
        if (kSym == lSym) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
        off_PUVX(iSym,jSym,kSym) = iStack
        iStack = iStack+iOrb*jAsh*kl_Orb_pairs
      end if
    end do
  end do
end do
nPUVX = iStack

! Init Fock matrices
FI(1:nTot1) = Zero
FA(1:nTot1) = Zero

! start transformation section

if (IPR >= 5) then
  write(u6,*) ' Symmetry  Basis functions   total orbitals    active orbitals'
  write(u6,*) ' -------------------------------------------------------------'
end if
do iSym=1,nSym
  iBas = nBas(iSym)
  iOrb = nOrb(iSym)
  iFro = nFro(iSym)
  iIsh = nIsh(iSym)
  iAsh = nAsh(iSym)
  do jSym=1,iSym
    jBas = nBas(jSym)
    jOrb = nOrb(jSym)
    jFro = nFro(jSym)
    jIsh = nIsh(jSym)
    jAsh = nAsh(jSym)
    ijSym = Mul(iSym,jSym)
    kSymMax = nSym
    if (.not. lSquare) kSymMax = iSym
    do kSym=1,kSymMax
      kBas = nBas(kSym)
      kOrb = nOrb(kSym)
      kFro = nFro(kSym)
      kIsh = nIsh(kSym)
      kAsh = nAsh(kSym)
      lSym = Mul(ijSym,kSym)
      if (lSym <= kSym) then
        lBas = nBas(lSym)
        lOrb = nOrb(lSym)
        lFro = nFro(lSym)
        lIsh = nIsh(lSym)
        lAsh = nAsh(lSym)

        if (iBas*jBas*kBas*lBas /= 0) then
          ij_Bas_pairs = iBas*jBas
          ij_Orb_pairs = iAsh*jAsh
          if (iSym == jSym) then
            ij_Bas_pairs = (iBas*iBas+iBas)/2
            ij_Orb_pairs = (iAsh*iAsh+iAsh)/2
          end if
          kl_Bas_pairs = kBas*lBas
          kl_Orb_pairs = kAsh*lAsh
          if (kSym == lSym) then
            kl_Bas_pairs = (kBas*kBas+kBas)/2
            kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
          end if
          call TraDrv(IPR,lSquare,iSym,jSym,kSym,lSym,iBas,jBas,kBas,lBas,iOrb,jOrb,kOrb,lOrb,iFro,jFro,kFro,lFro,iIsh,jIsh,kIsh, &
                      lIsh,iAsh,jAsh,kAsh,lAsh,ij_Bas_pairs,kl_Bas_pairs,ij_Orb_pairs,kl_Orb_pairs,off_PUVX,off_sqMat,off_ltMat, &
                      mxSym,CMO,PUVX,D1I,FI,D1A,FA,ExFac)
        end if

      end if
    end do
  end do
end do
if (IPR >= 5) then
  write(u6,*) ' -------------------------------------------------------------'
end if

! Synchronize Fock matrices if running parallel:
call GADsum(FI,nTot1)
call GADsum(FA,nTot1)

! print FI and FA
if (IPR >= 10) then
  write(u6,*)
  write(u6,*) ' FI in AO-basis'
  write(u6,*) ' --------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    if (iOrb > 0) call TriPrt(' ',' ',FI(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
  write(u6,*)
  write(u6,*) ' FA in AO-basis'
  write(u6,*) ' --------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    if (iOrb > 0) call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+(iOrb*iOrb+iOrb)/2
  end do
end if

! Synchronize PUVX if running parallel:
call GAdsum(PUVX,nPUVX)

! select integrals TUVX
call Get_TUVX(PUVX,TUVX)

! save integrals on disk
iDisk = 0
call DDaFile(LUINTM,1,PUVX,nPUVX,iDisk)

return

end subroutine Tra_Ctl2
