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
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This is a very preliminary routine for testing purposes. It relies   *
! on the basis set being of ANO type.                                  *
!                                                                      *
!                                                                      *
! Absolutely NOT to be used for production!!!!                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************

!#ifdef _DEBUGPRINT_
subroutine Fmod1n()

use Index_Functions, only: nTri_Elem
use GuessOrb_Global, only: Label, MxBasis, nBas, nSym, PrintMOs
use OneDat, only: sNoOri
use Molcas, only: MxSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: i, iBas, iComp, iDummy(7,8), iOff, iOpt, ipCMO(MxSym), ipFock(MxSym), ipX, ipY, iRc, iSym, iSymlb, jBas, k, &
                     kBas, Lu, nBasMax, nBasTot, nSqrTot, nTriTot, RC
real(kind=wp) :: dsum, eps, orbene(MxBasis), Sik, Sjk
real(kind=wp), allocatable :: Aux1(:), CMO(:), EVec(:), Fock(:), Hlf(:), Nrm(:), Ovl(:), SFk(:), TFk(:)
character(len=80) :: Title
character(len=8) :: Lbl

!----------------------------------------------------------------------*
! Setup various counters.                                              *
!----------------------------------------------------------------------*
nBasMax = 0
nBasTot = 0
nTriTot = 0
nSqrTot = 0
do iSym=1,nSym
  if (nBasMax < nBas(iSym)) nBasMax = nBas(iSym)
  nTriTot = nTriTot+nTri_Elem(nBas(iSym))
  nSqrTot = nSqrTot+nBas(iSym)*nBas(iSym)
  nBasTot = nBasTot+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Make symmetric orthonormal orbital basis.                            *
!----------------------------------------------------------------------*
call mma_allocate(CMO,nSqrTot)
ipCMO(1) = 1
do iSym=1,nSym-1
  ipCMO(iSym+1) = ipCMO(iSym)+nBas(iSym)*nBas(iSym)
end do
call goLowdin(CMO)
!----------------------------------------------------------------------*
! Allocate Fock matrix.                                                *
!----------------------------------------------------------------------*
call mma_allocate(Fock,nTriTot)
ipFock(1) = 1
call mma_allocate(EVec,nBasTot)
do iSym=1,nSym-1
  ipFock(iSym+1) = ipFock(iSym)+nTri_Elem(nBas(iSym))
end do
Fock(:) = Zero
!----------------------------------------------------------------------*
! Create model Fock operator.                                          *
!----------------------------------------------------------------------*
call FockOper(RC,EVec)
if (RC /= 0) then
  call mma_deallocate(EVec)
  call mma_deallocate(Fock)
  call mma_deallocate(CMO)
  return
end if
!----------------------------------------------------------------------*
! Get overlap matrix and make normalization.                           *
!----------------------------------------------------------------------*
call mma_allocate(Ovl,nTriTot+4)
call mma_allocate(Nrm,nBasTot)
iSymlb = 1
iOpt = ibset(0,sNoOri)
Lbl = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Lbl,iComp,Ovl,iSymlb)

ipX = 1
ipY = 1
do iSym=1,nSym
  call goPickup(Ovl(ipX),Nrm(ipY),nBas(iSym))
  Nrm(ipY:ipY+nBas(iSym)-1) = One/sqrt(Nrm(ipY:ipY+nBas(iSym)-1))
  ipX = ipX+nTri_Elem(nBas(iSym))
  ipY = ipY+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Build Fock matrix.                                                   *
!----------------------------------------------------------------------*
ipX = 1
ipY = 1
iOff = 0
do iSym=1,nSym
  do iBas=1,nBas(iSym)
    do jBas=1,iBas
      dsum = Zero
      do kBas=1,nBas(iSym)
        eps = EVec(iOff+kBas)
        i = max(iBas,kBas)
        k = min(iBas,kBas)
        Sik = Ovl(ipX-1+nTri_Elem(i-1)+k)*Nrm(ipY-1+iBas)*Nrm(ipY-1+kBas)
        i = max(jBas,kBas)
        k = min(jBas,kBas)
        Sjk = Ovl(ipX-1+nTri_Elem(i-1)+k)*Nrm(ipY-1+jBas)*Nrm(ipY-1+kBas)
        dsum = dsum+eps*Sik*Sjk
      end do
      Fock(ipFock(iSym)-1+nTri_Elem(iBas-1)+jBas) = dsum
    end do
  end do
# ifdef _DEBUGPRINT_
  call TriPrt('Modified atomic Fock matrix','(12f12.6)',Fock(ipFock(iSym)),nBas(iSym))
# endif
  ipX = ipX+nTri_Elem(nBas(iSym))
  ipY = ipY+nBas(iSym)
  iOff = iOff+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Scale Fock matrix.                                                   *
!----------------------------------------------------------------------*
ipX = 1
do iSym=1,nSym
# ifdef _DEBUGPRINT_
  write(u6,*) '***'
  write(u6,*) '*** Symmetry',iSym
  write(u6,*) '***'
# endif
  do iBas=1,nBas(iSym)
    do jBas=1,iBas
      Fock(ipFock(iSym)-1+nTri_Elem(iBas-1)+jBas) = Fock(ipFock(iSym)-1+nTri_Elem(iBas-1)+jBas)* &
                                                    sqrt(Ovl(ipX-1+nTri_Elem(iBas)))*sqrt(Ovl(ipX-1+nTri_Elem(jBas)))
    end do
  end do
# ifdef _DEBUGPRINT_
  call TriPrt('Scaled atomic Fock matrix','(12f12.6)',Fock(ipFock(iSym)),nBas(iSym))
# endif
  ipX = ipX+nTri_Elem(nBas(iSym))
end do
!----------------------------------------------------------------------*
! Release overlap and norm arrays.                                     *
!----------------------------------------------------------------------*
call mma_deallocate(Nrm)
call mma_deallocate(Ovl)
!----------------------------------------------------------------------*
! Diagonalize Fock matrix.                                             *
!----------------------------------------------------------------------*
call mma_allocate(SFk,nBasMax*nBasMax)
call mma_allocate(Hlf,nBasMax*nBasMax)
call mma_allocate(TFk,nTri_Elem(nBasMax))
iOff = 0
do iSym=1,nSym
# ifdef _DEBUGPRINT_
  write(u6,*) '***'
  write(u6,*) '*** Symmetry',iSym
  write(u6,*) '***'
# endif
  if (nBas(iSym) > 0) then
    call Square(Fock(ipFock(iSym)),SFk,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),One,SFk,nBas(iSym),CMO(ipCMO(iSym)),nBas(iSym),Zero,Hlf,nBas(iSym))
    call DGEMM_Tri('T','N',nBas(iSym),nBas(iSym),nBas(iSym),One,CMO(ipCMO(iSym)),nBas(iSym),Hlf,nBas(iSym),Zero,TFk,nBas(iSym))
#   ifdef _DEBUGPRINT_
    call TriPrt('Transformed Fock matrix','(12f12.6)',TFk,nBas(iSym))
#   endif
  end if

  !call Jacob(TFk,CMO(ipCMO(iSym)),nbas(iSym),nbas(iSym))
  call NIdiag(TFk,CMO(ipCMO(iSym)),nbas(iSym),nbas(iSym))
# ifdef _DEBUGPRINT_
  call TriPrt('Diagonalized atomic Fock matrix','(12f12.6)',TFk,nBas(iSym))
# endif
  call goPickup(TFk,orbene(iOff+1),nBas(iSym))
  call goSort(orbene(iOff+1),CMO(ipCMO(iSym)),nBas(iSym),nBas(iSym))
  iOff = iOff+nBas(iSym)
end do
call mma_deallocate(TFk)
call mma_deallocate(Hlf)
call mma_deallocate(SFk)
!----------------------------------------------------------------------*
! Present data.                                                        *
!----------------------------------------------------------------------*
if (PrintMOs) call PriMO('Start orbitals',.false.,.true.,Zero,1.0e6_wp,nSym,nBas,nBas,Label,orbene,orbene,CMO(ipCMO(1)),3)
call put_darray('Guessorb',CMO(ipCMO(1)),nSqrTot)
call put_darray('Guessorb energies',orbene,nBasTot)
call Put_iArray('nOrb',nBas,nSym)
call mma_allocate(Aux1,nBasTot)
Aux1(:) = Zero
Lu = 20
Title = 'Guess orbitals'
call WrVec('GSSORB',Lu,'COE',NSYM,NBAS,NBAS,CMO(ipCMO(1)),Aux1,orbene,iDummy,Title)
call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Done, deallocate the rest.                                           *
!----------------------------------------------------------------------*
call mma_deallocate(Evec)
call mma_deallocate(Fock)
call mma_deallocate(CMO)

end subroutine Fmod1n
