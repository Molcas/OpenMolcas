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
! This routine creates a symmetric ON basis a la Lowdin                *
! Not true if orbitals are deleted.                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine goLowdin(CMO)

use Index_Functions, only: nTri_Elem
use GuessOrb_Global, only: nBas, nDel, nSym, SThr
use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CMO(*)
integer(kind=iwp) :: iBas, iComp, iOff, iOpt, iOrb, ipCMO, ipOvl(8), irc, iSym, iSymlb, jBas, kBas, nBig, nElem, npSmat, nTot, &
                     nTriTot
real(kind=wp) :: Temp
character(len=8) :: Lbl
real(kind=wp), allocatable :: Eig(:), Ovl(:), SMat(:), Vec(:)
#ifdef _DEBUGPRINT_
real(kind=wp), allocatable :: Tmp(:,:)
#endif

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
nBig = 0
nTot = 0
nTriTot = 0
do iSym=1,nSym
  nTot = nTot+nBas(iSym)
  nTriTot = nTriTot+nTri_Elem(nBas(iSym))
  if (nBig < nBas(iSym)) nBig = nBas(iSym)
end do
!----------------------------------------------------------------------*
! Get overlap matrix                                                   *
!----------------------------------------------------------------------*
npSmat = nBig*nBig
call mma_allocate(Ovl,nTriTot+4)
ipOvl(1) = 1
call mma_allocate(Smat,npSmat)
iSymlb = 1
iOpt = ibset(0,sNoOri)
Lbl = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Lbl,iComp,Ovl,iSymlb)
do iSym=1,nSym-1
  ipOvl(iSym+1) = ipOvl(iSym)+nTri_Elem(nBas(iSym))
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call mma_allocate(Vec,nBig**2)
call mma_allocate(Eig,nBig)

ipCMO = 1
do iSym=1,nSym
  nElem = nTri_Elem(nBas(iSym))
  Smat(1:nElem) = Ovl(ipOvl(iSym):ipOvl(iSym)+nElem-1)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) '***'
  write(u6,*) '*** lowdin: symmetry',iSym
  write(u6,*) '***'
  write(u6,*)
  call TriPrt('Overlap matrix','(12f18.12)',Ovl(ipOvl(iSym)),nBas(iSym))
# endif
  call unitmat(Vec,nBas(iSym))
  call NIdiag_New(Ovl(ipOvl(iSym)),Vec,nBas(iSym),nbas(iSym))

  do iBas=1,nBas(iSym)
    call VecPhase(Vec((iBas-1)*nBas(iSym)+1),nBas(iSym))
  end do

# ifdef _DEBUGPRINT_
  call RecPrt('Transformation','(12f18.12)',Vec,nBas(iSym),nBas(iSym))
# endif
  call goPickUp(Ovl(ipOvl(iSym)),Eig,nBas(iSym))
# ifdef _DEBUGPRINT_
  call RecPrt('Overlap eigenvalues before sort','(12f18.12)',Eig,1,nBas(iSym))
# endif
  Eig(1:nBas(iSym)) = -Eig(1:nBas(iSym))
  call goSort(Eig,Vec,nBas(iSym),nBas(iSym))
  Eig(1:nBas(iSym)) = -Eig(1:nBas(iSym))
# ifdef _DEBUGPRINT_
  call RecPrt('Overlap eigenvalues after sort','(12f18.12)',Eig,1,nBas(iSym))
# endif
  nDel(iSym) = 0
  do iBas=1,nBas(iSym)
    if (Eig(iBas) < SThr) nDel(iSym) = nDel(iSym)+1
  end do
  Eig(1:nBas(iSym)) = One/sqrt(Eig(1:nBas(iSym)))
  if (.false.) then
    do iBas=1,nBas(iSym)
      do jBas=1,nBas(iSym)
        Temp = Zero
        do kBas=1,nBas(iSym)
          Temp = Temp+Eig(kBas)*Vec(nBas(iSym)*(kBas-1)+iBas)*Vec(nBas(iSym)*(kBas-1)+jBas)
        end do
        CMO(ipCMO+nBas(iSym)*(jBas-1)+(iBas-1)) = Temp
      end do
    end do
  else
    nElem = nBas(iSym)**2
    CMO(ipCMO:ipCMO+nElem-1) = Vec(1:nElem)
    do iOrb=1,nBas(iSym)
      iOff = ipCMO+nBas(iSym)*(iOrb-1)
      CMO(iOff:iOff+nBas(iSym)-1) = Eig(iOrb)*CMO(iOff:iOff+nBas(iSym)-1)
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('Symmetric orbitals','(12f18.12)',CMO(ipCMO),nBas(iSym),nBas(iSym))
  call mma_allocate(Tmp,nBas(iSym),nBas(iSym))
  iSymlb = 1
  Lbl = 'Mltpl  0'
  call RdOne(irc,iOpt,Lbl,iComp,Ovl(ipOvl(1)),iSymlb)
  do iBas=1,nBas(iSym)
    do jBas=1,nBas(iSym)
      Temp = Zero
      do kBas=1,nBas(iSym)
        Temp = Temp+CMO(ipCMO+nBas(iSym)*(kBas-1)+(iBas-1))*CMO(ipCMO+nBas(iSym)*(jBas-1)+(kBas-1))
      end do
      Tmp(iBas,jBas) = Temp
    end do
  end do
  call RecPrt('Inverted overlap matrix','(12f18.12)',Tmp,nBas(iSym),nBas(iSym))
  call mma_deallocate(Tmp)
# endif
  ipCMO = ipCMO+nBas(iSym)*nBas(iSym)
end do

call mma_deallocate(Eig)
call mma_deallocate(Vec)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call mma_deallocate(Smat)
call mma_deallocate(Ovl)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*

end subroutine goLowdin
