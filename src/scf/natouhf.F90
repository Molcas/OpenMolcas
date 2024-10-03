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
! Copyright (C) 2006, Per-Olof Widmark                                 *
!***********************************************************************

subroutine NatoUHF(DensA,DensB,FockA,FockB,nBT,CMO,nBB,Ovl,Nato,Eta,Eps,nnB,nSym,nBas,nOrb)
!***********************************************************************
!                                                                      *
! This routine computes the natural orbitals for a UHF wavefunction.   *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBT, nBB, nnB, nSym, nBas(nSym), nOrb(nSym)
real(kind=wp), intent(in) :: DensA(nBT), DensB(nBT), FockA(nBT), FockB(nBT), CMO(nBB), Ovl(nBT)
real(kind=wp), intent(out) :: Nato(nBB), Eta(nnB), Eps(nnB)
integer(kind=iwp) :: indx, iOffCMO, iOffEta, iOffSqr, iOffTri, iOrb, iSym, MaxSqr, MaxTri, nCMO, nEta, nSqr, nTri
real(kind=wp) :: tmp
real(kind=wp), allocatable :: Aux1(:), Aux2(:), Aux3(:), Dens(:), Fock(:), SMat(:)

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
MaxTri = 0
MaxSqr = 0
nTri = 0
nSqr = 0
nCMO = 0
nEta = 0
do iSym=1,nSym
  MaxTri = max(MaxTri,nTri_Elem(nBas(iSym)))
  MaxSqr = max(MaxSqr,nBas(iSym)**2)
  nTri = nTri+nTri_Elem(nBas(iSym))
  nSqr = nSqr+nBas(iSym)**2
  nCMO = nCMO+nBas(iSym)*nOrb(iSym)
  nEta = nEta+nOrb(iSym)
end do
!----------------------------------------------------------------------*
! Allocate arrays                                                      *
!----------------------------------------------------------------------*
call mma_allocate(Dens,nTri,Label='Dens')
call mma_allocate(Fock,nTri,Label='Fock')
call mma_allocate(SMat,MaxSqr,Label='SMat')
call mma_allocate(Aux1,MaxSqr,Label='Aux1')
call mma_allocate(Aux2,MaxSqr,Label='Aux2')
call mma_allocate(Aux3,MaxSqr,Label='Aux3')
!----------------------------------------------------------------------*
! Add up the densities                                                 *
!----------------------------------------------------------------------*
Dens(:) = DensA(1:nTri)+DensB(1:nTri)
!----------------------------------------------------------------------*
! Copy orbitals                                                        *
!----------------------------------------------------------------------*
Nato(1:nCMO) = CMO(1:nCMO)
!----------------------------------------------------------------------*
! Average Fock matrix.                                                 *
!----------------------------------------------------------------------*
Fock(:) = Half*(FockA(1:nTri)+FockB(1:nTri))
!iOffTri = 0
!do iSym=1,nSym
!  call TriPrt('natouhf: Fock A','(20f10.4)',FockA(1+iOffTri),nBas(iSym))
!  call TriPrt('natouhf: Fock B','(20f10.4)',FockB(1+iOffTri),nBas(iSym))
!  call TriPrt('natouhf: Average Fock','(20f10.4)',FockB(1+iOffTri),nBas(iSym))
!  iOffTri = iOffTri+nTri_Elem(nBas(iSym))
!end do
!----------------------------------------------------------------------*
! Form density in MO basis: C(t)SDSC                                   *
!----------------------------------------------------------------------*
iOffSqr = 0
iOffTri = 0
iOffCMO = 0
iOffEta = 0
do iSym=1,nSym
  if (iSym > 1) then
    iOffTri = iOffTri+nTri_Elem(nBas(iSym-1))
    iOffSqr = iOffSqr+nBas(iSym-1)**2
    iOffEta = iOffEta+nOrb(iSym-1)
    iOffCMO = iOffCMO+nBas(iSym-1)*nOrb(iSym-1)
  end if
  if (nBas(iSym) <= 0) cycle

  ! Compute C(t)S (=aux1), but square S first

  call Square(Ovl(iOfftri+1),Smat,1,nBas(iSym),nBas(iSym))
  !call RecPrt('natouhf: Smat','(12f12.6)',Smat,nBas(iSym),nBas(iSym))
  call DGEMM_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym), &
              One,CMO(iOffCMO+1),nBas(iSym), &
              Smat,nBas(iSym), &
              Zero,Aux1,nOrb(iSym))
  !call RecPrt('natouhf: C(t)S','(12f12.6)',Aux1,nOrb(iSym),nBas(iSym))

  ! Compute C(t)SD (=aux3), but first unfold D (=aux2)

  call Dsq(Dens(1+iOffTri),Aux2,1,nBas(iSym),nBas(iSym))
  !call RecPrt('natouhf: Density','(12f12.6)',Aux2,nBas(iSym),nBas(iSym))
  call DGEMM_('N','N',nOrb(iSym),nBas(iSym),nBas(iSym), &
              One,Aux1,nOrb(iSym), &
              Aux2,nBas(iSym), &
              Zero,Aux3,nOrb(iSym))
  !call RecPrt('natouhf: C(t)SD','(12f12.6)',Aux3,nOrb(iSym),nBas(iSym))

  ! Compute C(t)SDS (=aux1)

  call DGEMM_('N','N',nOrb(iSym),nBas(iSym),nBas(iSym), &
              One,Aux3,nOrb(iSym), &
              Smat,nBas(iSym), &
              Zero,Aux1,nOrb(iSym))
  !call RecPrt('natouhf: C(t)SDS','(12f12.6)',Aux1,nOrb(iSym),nBas(iSym))

  ! Compute C(t)SDSC (=aux2)

  call DGEMM_Tri('N','N',nOrb(iSym),nOrb(iSym),nBas(iSym), &
                 One,Aux1,nOrb(iSym), &
                 CMO(iOffCMO+1),nBas(iSym), &
                 Zero,Aux2,nOrb(iSym))
  !call TriPrt('natouhf: C(t)SDSC','(12f12.6)',Aux2,nOrb(iSym))

  ! Shift diagonal slightly

  tmp = 1.0e-6_wp
  do iOrb=1,nOrb(iSym)
    indx = nTri_Elem(iOrb)
    Aux2(indx) = Aux2(indx)+tmp
    tmp = Half*tmp
  end do
  !call TriPrt('natouhf: C(t)SDSC + shift','(12f12.6)',Aux2,nOrb(iSym))

  ! Diagonalize

  call NIdiag(Aux2,Nato(iOffCMO+1),nOrb(iSym),nBas(iSym))
  call Pickup(Aux2,Eta(iOffEta+1),nOrb(iSym))
  Eta(iOffEta+1:iOffEta+nOrb(iSym)) = -Eta(iOffEta+1:iOffEta+nOrb(iSym))
  call SortEig(Eta(iOffEta+1),Nato(iOffCMO+1),nOrb(iSym),nBas(iSym),1,.true.)
  Eta(iOffEta+1:iOffEta+nOrb(iSym)) = -Eta(iOffEta+1:iOffEta+nOrb(iSym))
end do
!----------------------------------------------------------------------*
! Compute diagonal of average Fock matrix                              *
!----------------------------------------------------------------------*
call MkEorb_Inner(Fock,nBT,Nato,nBB,Eps,nnB,nSym,nBas,nOrb)
!call PriMO('natouhf: UHF nato',.true.,.true.,-One,1.0e6_wp,nSym,nBas,nOrb,Name,Eps,Eta,Nato,3)
!----------------------------------------------------------------------*
! Deallocate arrays                                                    *
!----------------------------------------------------------------------*
call mma_deallocate(Aux3)
call mma_deallocate(Aux2)
call mma_deallocate(Aux1)
call mma_deallocate(SMat)
call mma_deallocate(Fock)
call mma_deallocate(Dens)
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine NatoUHF
