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

subroutine Fmod1n(StandAlone)

use GuessOrb_Global, only: Label, MxBasis, MxSym, nBas, nSym, PrintMOs

implicit none
#include "stdalloc.fh"
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
logical StandAlone
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical Debug
logical Trace
integer iSym
integer iBas
integer jBas
integer kBas
integer iOff
integer ipCMO(MxSym)
integer ipFock(MxSym)
integer ipX
integer ipY
integer iSymlb
integer iRc
integer nBasMax
integer nBasTot
integer nTriTot
integer nSqrTot
!---
real*8 orbene(MxBasis)
real*8,dimension(:),allocatable :: CMO,Fock,EVec,Ovl,Nrm
real*8,dimension(:),allocatable :: SFk,Hlf,TFk,Aux1
!character(len=4) Name(MxAtom)
!character(len=4) Label(2,MxBasis)
!---
integer i,k
real*8 Sik
real*8 Sjk
real*8 eps
real*8 sum
integer Lu
integer iDummy(7,8)
integer RC
character*80 Title
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
if (StandAlone) then
  Debug = .false.
  Trace = .false.
else
  Debug = .false.
  Trace = .false.
end if
if (Trace) write(6,*) '>>> Entering fmod1n'
!----------------------------------------------------------------------*
! Setup various counters.                                              *
!----------------------------------------------------------------------*
nBasMax = 0
nBasTot = 0
nTriTot = 0
nSqrTot = 0
do iSym=1,nSym
  if (nBasMax < nBas(iSym)) nBasMax = nBas(iSym)
  nTriTot = nTriTot+nBas(iSym)*(nBas(iSym)+1)/2
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
  ipFock(iSym+1) = ipFock(iSym)+nBas(iSym)*(nBas(iSym)+1)/2
end do
call dCopy_(nTriTot,[0.0d0],0,Fock,1)
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
call RdOne(irc,2,'Mltpl  0',1,Ovl,iSymlb)

ipX = 1
ipY = 1
do iSym=1,nSym
  call goPickup(Ovl(ipX),Nrm(ipY),nBas(iSym))
  do iBas=1,nBas(iSym)
    Nrm(ipY-1+iBas) = 1.0d0/sqrt(Nrm(ipY-1+iBas))
  end do
  ipX = ipX+nBas(iSym)*(nBas(iSym)+1)/2
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
      sum = 0.0d0
      do kBas=1,nBas(iSym)
        eps = EVec(iOff+kBas)
        i = max(iBas,kBas)
        k = min(iBas,kBas)
        Sik = Ovl(ipX-1+i*(i-1)/2+k)*Nrm(ipY-1+iBas)*Nrm(ipY-1+kBas)
        i = max(jBas,kBas)
        k = min(jBas,kBas)
        Sjk = Ovl(ipX-1+i*(i-1)/2+k)*Nrm(ipY-1+jBas)*Nrm(ipY-1+kBas)
        sum = sum+eps*Sik*Sjk
      end do
      Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas) = sum
    end do
  end do
  if (Debug) then
    call TriPrt('Modified atomic Fock matrix','(12f12.6)',Fock(ipFock(iSym)),nBas(iSym))
  end if
  ipX = ipX+nBas(iSym)*(nBas(iSym)+1)/2
  ipY = ipY+nBas(iSym)
  iOff = iOff+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Scale Fock matrix.                                                   *
!----------------------------------------------------------------------*
ipX = 1
do iSym=1,nSym
  if (Debug) then
    write(6,*) '***'
    write(6,*) '*** Symmetry',iSym
    write(6,*) '***'
  end if
  do iBas=1,nBas(iSym)
    do jBas=1,iBas
      Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas) = Fock(ipFock(iSym)-1+iBas*(iBas-1)/2+jBas)*sqrt(Ovl(ipX-1+iBas*(iBas+1)/2))* &
                                                  sqrt(Ovl(ipX-1+jBas*(jBas+1)/2))
    end do
  end do
  if (Debug) then
    call TriPrt('Scaled atomic Fock matrix','(12f12.6)',Fock(ipFock(iSym)),nBas(iSym))
  end if
  ipX = ipX+nBas(iSym)*(nBas(iSym)+1)/2
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
call mma_allocate(TFk,nBasMax*(nBasMax+1)/2)
iOff = 0
do iSym=1,nSym
  if (Debug) then
    write(6,*) '***'
    write(6,*) '*** Symmetry',iSym
    write(6,*) '***'
  end if
  if (nBas(iSym) > 0) then
    call Square(Fock(ipFock(iSym)),SFk,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nBas(iSym),1.0d0,SFk,nBas(iSym),CMO(ipCMO(iSym)),nBas(iSym),0.0d0,Hlf,nBas(iSym))
    call MxMt(CMO(ipCMO(iSym)),nBas(iSym),1,Hlf,1,nBas(iSym),TFk,nBas(iSym),nBas(iSym))
    if (Debug) then
      call TriPrt('Transformed Fock matrix','(12f12.6)',TFk,nBas(iSym))
    end if
  end if

  !call Jacob(TFk,CMO(ipCMO(iSym)),nbas(iSym),nbas(iSym))
  call NIdiag(TFk,CMO(ipCMO(iSym)),nbas(iSym),nbas(iSym),0)
  if (Debug) then
    call TriPrt('Diagonalized atomic Fock matrix','(12f12.6)',TFk,nBas(iSym))
  end if
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
if (PrintMOs) then
  call PriMO('Start orbitals',.false.,.true.,0.0d0,1.0d6,nSym,nBas,nBas,Label,orbene,orbene,CMO(ipCMO(1)),3)
end if
call put_darray('Guessorb',CMO(ipCMO(1)),nSqrTot)
call put_darray('Guessorb energies',orbene,nBasTot)
call Put_iArray('nOrb',nBas,nSym)
call mma_allocate(Aux1,nBasTot)
do iBas=1,nBasTot
  Aux1(iBas) = 0.0d0
end do
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
if (trace) write(6,*) '<<< Exiting fmod1n'

return

end subroutine Fmod1n
