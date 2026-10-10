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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine OLagFro1(NBSQT,nOLag,Trf,DFro,OLag)

use caspt2_global, only: FIFA_all, if_FroSCF
use caspt2_module, only: NBAS, NBAST, NDEL, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT, nOLag
real(kind=wp), intent(in) :: Trf(NBSQT)
real(kind=wp), intent(out) :: DFro(NBSQT)
real(kind=wp), intent(inout) :: OLag(nOLag)
integer(kind=iwp) :: iEps, iMO, iOrb, iSQ, iSym, jOrb, nBasI, nFroI, nIshI, nOrbI
real(kind=wp), allocatable :: EPS(:), WRK(:,:), ZPT2(:,:)

! The frozen-inactive density (DFro) is the multiplier of the condition that the Fock matrix is zero
! between the frozen and inactive orbitals.
! The frozen orbitals are determined by the Fock matrix constructed with the RDM of RASSCF.
! With if_FroSCF, the orbital energies of RASSCF are used, and Trf transforms the PT2 orbitals to the canonical ones of RASSCF.
! Otherwise, the Fock matrix of CASPT2 is used, which is diagonal in the frozen and inactive blocks.

call mma_allocate(EPS,nBasT,Label='EPS')
if (if_FroSCF) then
  ! (canonical) orbital energies of RASSCF
  call Get_dArray('RASSCF OrbE',EPS,nBasT)
else
  ! (canonical) orbital energies of PT2
  iMO = 0
  iEps = 0
  do iSym=1,nSym
    nOrbI = nBas(iSym)-nDel(iSym)
    do iOrb=1,nOrbI
      EPS(iEps+iOrb) = FIFA_all(iMO+iOrb+nOrbI*(iOrb-1))
    end do
    iMO = iMO+nOrbI*nOrbI
    iEps = iEps+nBas(iSym)
  end do
end if

DFro(:) = Zero
iMO = 1
iSQ = 0
iEps = 0
do iSym=1,nSym
  nBasI = nBas(iSym)
  nOrbI = nBas(iSym)-nDel(iSym)
  nFroI = nFro(iSym)
  if ((nOrbI > 0) .and. (nFroI > 0)) then
    nIshI = nIsh(iSym)
    call mma_allocate(ZPT2,nFroI,nIshI,Label='ZPT2')
    call mma_allocate(WRK,nFroI,nIshI,Label='WRK')
    !! Make sure that the frozen orbital derivative is zero
    !! (it does not appear in the PT2 energy)
    OLag(iMO:iMO+nOrbI*nFroI-1) = Zero
    do jOrb=1,nIshI
      do iOrb=1,nFroI
        ZPT2(iOrb,jOrb) = -Half*(OLag(iMO+iOrb-1+nOrbI*(nFroI+jOrb-1))-OLag(iMO+nFroI+jOrb-1+nOrbI*(iOrb-1)))
      end do
    end do
    if (if_FroSCF) then
      ! PT2 orbitals --> canonical orbitals of RASSCF
      call DGEMM_('N','N',nFroI,nIshI,nFroI,One,Trf(iSQ+1),nBasI,ZPT2,nFroI,Zero,WRK,nFroI)
      call DGEMM_('N','T',nFroI,nIshI,nIshI,One,WRK,nFroI,Trf(iSQ+nFroI+1+nBasI*nFroI),nBasI,Zero,ZPT2,nFroI)
    end if
    do jOrb=1,nIshI
      do iOrb=1,nFroI
        ZPT2(iOrb,jOrb) = ZPT2(iOrb,jOrb)/(EPS(iEps+iOrb)-EPS(iEps+nFroI+jOrb))
      end do
    end do
    if (if_FroSCF) then
      ! canonical orbitals of RASSCF --> PT2 orbitals
      call DGEMM_('T','N',nFroI,nIshI,nFroI,One,Trf(iSQ+1),nBasI,ZPT2,nFroI,Zero,WRK,nFroI)
      call DGEMM_('N','N',nFroI,nIshI,nIshI,One,WRK,nFroI,Trf(iSQ+nFroI+1+nBasI*nFroI),nBasI,Zero,ZPT2,nFroI)
    end if
    do jOrb=1,nIshI
      do iOrb=1,nFroI
        DFro(iMO+iOrb-1+nOrbI*(nFroI+jOrb-1)) = ZPT2(iOrb,jOrb)
        DFro(iMO+nFroI+jOrb-1+nOrbI*(iOrb-1)) = ZPT2(iOrb,jOrb)
      end do
    end do
    call mma_deallocate(ZPT2)
    call mma_deallocate(WRK)
  end if
  iMO = iMO+nOrbI*nOrbI
  iSQ = iSQ+nBasI*nBasI
  iEps = iEps+nBasI
end do
call mma_deallocate(EPS)

end subroutine OLagFro1
