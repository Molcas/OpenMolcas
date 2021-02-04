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

subroutine goLowdin(CMO)

use GuessOrb_Global, only: nBas, nDel, nSym, SThr
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
!----------------------------------------------------------------------*
! Dummy variables.                                                     *
!----------------------------------------------------------------------*
real(kind=wp), intent(out) :: CMO(*)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical(kind=iwp) :: Debug, Trace
integer(kind=iwp) :: nBig, nTot, nTri, nTriTot, iSym, iBas, jBas, kBas, iOrb, ipOvl(8), ipCMO, npSmat, irc, iSymlb
real(kind=wp) :: Temp, OrbPhase
real(kind=wp), allocatable :: Ovl(:), SMat(:), Vec(:), Eig(:), Tmp(:,:)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
Debug = .false.
Trace = .false.
if (Trace) write(u6,*) '>>> Entering golowdin'
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
nBig = 0
nTot = 0
nTriTot = 0
do iSym=1,nSym
  nTot = nTot+nBas(iSym)
  nTriTot = nTriTot+nBas(iSym)*(nBas(iSym)+1)/2
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
call RdOne(irc,2,'Mltpl  0',1,Ovl,iSymlb)
do iSym=1,nSym-1
  ipOvl(iSym+1) = ipOvl(iSym)+nBas(iSym)*(nBas(iSym)+1)/2
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call mma_allocate(Vec,nBig**2)
call mma_allocate(Eig,nBig)

ipCMO = 1
do iSym=1,nSym
  nTri = nBas(iSym)*(nBas(iSym)+1)/2
  call dCopy_(nTri,Ovl(ipOvl(iSym)),1,Smat,1)
  if (Debug) then
    write(u6,*)
    write(u6,*) '***'
    write(u6,*) '*** lowdin: symmetry',iSym
    write(u6,*) '***'
    write(u6,*)
    call TriPrt('Overlap matrix','(12f18.12)',Ovl(ipOvl(iSym)),nBas(iSym))
  end if
  call FZero(Vec,nBas(iSym)**2)
  call DCopy_(nBas(iSym),[One],0,Vec,nBas(iSym)+1)
  call NIdiag_New(Ovl(ipOvl(iSym)),Vec,nBas(iSym),nbas(iSym),0)

  do iBas=1,nBas(iSym)
    temp = OrbPhase(Vec((iBas-1)*nBas(iSym)+1),nBas(iSym))
  end do

  if (Debug) then
    call RecPrt('Transformation','(12f18.12)',Vec,nBas(iSym),nBas(iSym))
  end if
  call goPickUp(Ovl(ipOvl(iSym)),Eig,nBas(iSym))
  if (Debug) then
    call RecPrt('Overlap eigenvalues before sort','(12f18.12)',Eig,1,nBas(iSym))
  end if
  do iBas=1,nBas(iSym)
    Eig(iBas) = -Eig(iBas)
  end do
  call goSort(Eig,Vec,nBas(iSym),nBas(iSym))
  do iBas=1,nBas(iSym)
    Eig(iBas) = -Eig(iBas)
  end do
  if (Debug) then
    call RecPrt('Overlap eigenvalues after sort','(12f18.12)',Eig,1,nBas(iSym))
  end if
  nDel(iSym) = 0
  do iBas=1,nBas(iSym)
    if (Eig(iBas) < SThr) nDel(iSym) = nDel(iSym)+1
  end do
  do iBas=1,nBas(iSym)
    Eig(iBas) = One/sqrt(Eig(iBas))
  end do
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
    call dCopy_(nBas(iSym)*nBas(iSym),Vec,1,CMO(ipCMO),1)
    do iOrb=1,nBas(iSym)
      Temp = Eig(iOrb)
      do iBas=1,nBas(iSym)
        CMO(ipCMO-1+iBas+nBas(iSym)*(iOrb-1)) = Temp*CMO(ipCMO-1+iBas+nBas(iSym)*(iOrb-1))
      end do
    end do
  end if
  if (Debug) then
    call RecPrt('Symmetric orbitals','(12f18.12)',CMO(ipCMO),nBas(iSym),nBas(iSym))
  end if
  if (Debug) then
    call mma_allocate(Tmp,nBas(iSym),nBas(iSym))
    iSymlb = 1
    call RdOne(irc,2,'Mltpl  0',1,Ovl(ipOvl(1)),iSymlb)
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
  end if
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
if (Trace) write(u6,*) '<<< Exiting golowdin'

return

end subroutine goLowdin
