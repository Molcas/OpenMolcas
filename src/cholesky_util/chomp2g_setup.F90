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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_Setup(irc,EOcc,EVir)
!
! Jonas Bostrom, Feb 2010
!
! Purpose: Do some additional setup only needed for
!          MP2-gradients or properties.

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use ChoMP2, only: AdrR1, AdrR2, ChoMP2g_Allocated, EFrozT, EOccuT, EVirtT, iAdrOff, iAoMo, iDel, iFro, iMoAo, iMoMo, iOcc, iVir, &
                  MP2D, MP2D_e, MP2D_e_full, MP2D_full, MP2W, MP2W_e, MP2W_e_full, MP2W_full, nAoMo, nDel, nFro, nFroT, nMo, &
                  nMoAo, nMoMo, nMoType, nOcc, nOccT, nOccVirt, nOrb, nVir, nVirT
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: iE, iMOType, iProdType, iS, iSym, iSymAl, iSymP, iSymQ, jMoType, lDens, lDens_e, nb

nMOType = 3
call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)

! Initialize an  offset for writing choleskyvectors to disk
iAdrOff(1:nSym,1:nMoType**2) = 0

nOccVirT = sum(nOcc(1:nSym)*nVir(1:nSym))

nMO(1:nSym,1) = nFro(1:nSym)
nMO(1:nSym,2) = nOcc(1:nSym)
nMO(1:nSym,3) = nVir(1:nSym)

do iMoType=1,nMoType
  do jMoType=1,nMoType
    iProdType = jMOType+(iMOType-1)*nMOType
    do iSym=1,nSym
      nMoMo(iSym,iProdType) = 0
      do iSymP=1,nSym
        iSymQ = Mul(iSymP,iSym)
        iMoMo(iSymq,iSymp,iProdType) = nMoMo(iSym,iProdType)
        nMoMo(iSym,iProdType) = nMoMo(iSym,iProdType)+nMO(iSymP,iMOType)*nMO(iSymQ,jMOType)
      end do
    end do
  end do
end do
do iMoType=1,nMoType
  do iSym=1,nSym
    nMoAo(iSym,iMoType) = 0
    do iSymAl=1,nSym
      iSymP = Mul(iSymAl,iSym)
      iMoAo(iSymP,iSymAl,iMoType) = nMoAo(iSym,iMoType)
      nMoAo(iSym,iMoType) = nMoAo(iSym,iMoType)+nBas(iSymAl)*nMO(iSymP,iMoType)
    end do
  end do
end do
do iMoType=1,nMoType
  do iSym=1,nSym
    nAoMo(iSym,iMoType) = 0
    do iSymQ=1,nSym
      iSymAl = Mul(iSymQ,iSym)
      iAoMo(iSymAl,iSymQ,iMoType) = nAoMo(iSym,iMoType)
      nAoMo(iSym,iMoType) = nAoMo(iSym,iMoType)+nBas(iSymAl)*nMO(iSymQ,iMoType)
    end do
  end do
end do

! Allocate MP2_density
! --------------------

lDens = sum(nOrb(1:nSym)*nOrb(1:nSym))

ChoMP2g_Allocated = .true.

call mma_allocate(MP2D_full,lDens,Label='MP2D_full')
call mma_allocate(MP2W_full,lDens,Label='MP2W_full')
MP2D_full(:) = Zero
MP2W_full(:) = Zero

iE = 0
do iSym=1,nSym
  nb = nOrb(iSym)
  iS = iE+1
  iE = iE+nb**2
  MP2D(iSym)%A(1:nb,1:nb) => MP2D_full(iS:iE)
  MP2W(iSym)%A(1:nb,1:nb) => MP2W_full(iS:iE)
end do

! Allocate extended MP2_density (with deleted orbitals)
! -----------------------------------------------------

lDens_e = sum((nOrb(1:nSym)+nDel(1:nSym))*(nOrb(1:nSym)+nDel(1:nSym)))

call mma_allocate(MP2D_e_full,lDens_e,Label='MP2D_e_full')
call mma_allocate(MP2W_e_full,lDens_e,Label='MP2W_e_full')
MP2D_e_full(:) = Zero
MP2W_e_full(:) = Zero

iE = 0
do iSym=1,nSym
  nb = nOrb(iSym)+nDel(iSym)
  iS = iE+1
  iE = iE+nb**2
  MP2D_e(iSym)%A(1:nb,1:nb) => MP2D_e_full(iS:iE)
  MP2W_e(iSym)%A(1:nb,1:nb) => MP2W_e_full(iS:iE)
end do

! Allocate adress-field for reordered R-vectors
! ---------------------------------------------
call mma_allocate(AdrR1,nSym,nSym,nOccT,Label='AdrR1')
call mma_allocate(AdrR2,nSym,nSym,nVirT,Label='AdrR2')

! Allocate a vector for the orbital energies of frozen and virtual
! frozen molecules.
call mma_allocate(EFrozT,max(1,nFroT),Label='EFrozT')
call mma_allocate(EOccuT,max(1,nOccT),Label='EOccuT')
call mma_allocate(EVirtT,max(1,nVirT),Label='EVirtT')
! Fill them with the right things
do iSym=1,nSym
  EFrozT(iFro(iSym)+1:iFro(iSym)+nFro(iSym)) = EOcc(iFro(iSym)+nOccT+1:iFro(iSym)+nOccT+nFro(iSym))
  EOccuT(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym)) = EOcc(iOcc(iSym)+1:iOcc(iSym)+nOcc(iSym))
  EVirtT(iVir(iSym)+1:iVir(iSym)+nVir(iSym)) = EVir(iVir(iSym)+iDel(iSym)+1:iVir(iSym)+iDel(iSym)+nVir(iSym))
end do

irc = 0

end subroutine ChoMP2g_Setup
