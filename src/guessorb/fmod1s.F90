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

!#define _DEBUGPRINT_
subroutine Fmod1s()

use Index_Functions, only: nTri_Elem
use GuessOrb_Global, only: Label, nBas, nSym, PrintMOs
use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iBas, iComp, iDummy(7,8), indx, iOpt, ipTmp1, ipTmp2, ipTmp3, ipTmp4, irc, iSym, iSymlb, jBas, kBas, Lu, &
                     n2Full, nBasMax, nBasTot, nSqrTot, nTriTot, RC
real(kind=wp) :: Det, dsum, eps, Sii, Sik, Sjj, Sjk, Skk
character(len=80) :: Title
character(len=8) :: Lbl
real(kind=wp), allocatable :: Aux1(:), Aux2(:), CMOs(:), DeTr(:), Edes(:), Esym(:), Evec(:), Fdes(:), Fmo(:), Fock(:), FSym(:), &
                              Ovl(:), Sdes(:), Smat(:), SmTr(:)
#ifdef _DEBUGPRINT_
character(len=32) :: Line
#endif

!----------------------------------------------------------------------*
! Should only be called if symmetry is used.                           *
!----------------------------------------------------------------------*
if (nSym == 1) call SysAbendMsg('fmod1s','internal error 001',' ')
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
n2Full = nBasTot*nBasTot
!----------------------------------------------------------------------*
! Get symmetry transformation matrix.                                  *
!----------------------------------------------------------------------*
call mma_allocate(SmTr,n2Full)
call mma_allocate(DeTr,n2Full)
call Get_dArray('SM',SmTr,n2Full)
call Minv(SmTr,DeTr,Det,nBasTot)
#ifdef _DEBUGPRINT_
call RecPrt('Symmetrization transformation matrix','(10f12.6)',SmTr,nBasTot,nBasTot)
call RecPrt('Desymmetrization transformation matrix','(10f12.6)',DeTr,nBasTot,nBasTot)
#endif
!----------------------------------------------------------------------*
! Create model Fock operator.                                          *
!----------------------------------------------------------------------*
call mma_allocate(Esym,nBasTot)
call FockOper(RC,Esym)
if (RC /= 0) then
  call mma_deallocate(Esym)
  call mma_deallocate(DeTr)
  call mma_deallocate(SmTr)
  return
end if
#ifdef _DEBUGPRINT_
call RecPrt('Diagonal of Fock','(10f12.6)',Esym,1,nBasTot)
#endif
!----------------------------------------------------------------------*
! Desymmetrize model Fock operator                                     *
!----------------------------------------------------------------------*
call mma_allocate(Edes,nBasTot)
do iBas=1,nBasTot
  do kBas=1,nBasTot
    indx = kBas+nBasTot*(iBas-1)
    if (abs(SmTr(indx)) > 1.0e-3_wp) Edes(kBas) = Esym(iBas)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Desymmetrized diagonal of Fock','(10f12.6)',Edes,1,nBasTot)
#endif
!----------------------------------------------------------------------*
! Read and square overlap matrix.                                      *
!----------------------------------------------------------------------*
call mma_allocate(Smat,n2Full)
call mma_allocate(Ovl,nTriTot)
iSymlb = 1
iOpt = ibset(0,sNoOri)
Lbl = 'Mltpl  0'
iComp = 1
call RdOne(irc,iOpt,Lbl,iComp,Ovl,iSymlb)
Smat(:) = Zero
ipTmp1 = 1
ipTmp2 = 1
do iSym=1,nSym
  call Square(Ovl(ipTmp1),SMat(ipTmp2),1,nBasTot,nBas(iSym))
  ipTmp1 = ipTmp1+nTri_Elem(nBas(iSym))
  ipTmp2 = ipTmp2+nBas(iSym)*nBasTot+nBas(iSym)
end do
#ifdef _DEBUGPRINT_
call RecPrt('Full overlap matrix','(10f12.6)',Smat,nBasTot,nBasTot)
#endif
call mma_deallocate(Ovl)
!----------------------------------------------------------------------*
! Desymmetrize overlap matrix.                                         *
!----------------------------------------------------------------------*
call mma_allocate(Sdes,n2Full)
call mma_allocate(Aux1,n2Full)
call DGEMM_('N','N',nBasTot,nBasTot,nBasTot,One,Smat,nBasTot,DeTr,nBasTot,Zero,Aux1,nBasTot)
call DGEMM_('T','N',nBasTot,nBasTot,nBasTot,One,DeTr,nBasTot,Aux1,nBasTot,Zero,Sdes,nBasTot)
#ifdef _DEBUGPRINT_
call RecPrt('Desymmetrized overlap matrix','(10f12.6)',Sdes,nBasTot,nBasTot)
#endif
call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Build Fock matrix                                                    *
!----------------------------------------------------------------------*
call mma_allocate(Fdes,n2Full)
do iBas=1,nBasTot
  Sii = Sdes(iBas+nBasTot*(iBas-1))
  do jBas=1,nBasTot
    Sjj = Sdes(jBas+nBasTot*(jBas-1))
    dsum = Zero
    do kBas=1,nBasTot
      eps = Edes(kBas)
      Skk = Sdes(kBas+nBasTot*(kBas-1))
      Sik = Sdes(iBas+nBasTot*(kBas-1))
      Sjk = Sdes(jBas+nBasTot*(kBas-1))
      dsum = dsum+eps*Sik*Sjk/sqrt(Sii*Sjj*Skk*Skk)
    end do
    Fdes(iBas+nBasTot*(jBas-1)) = dsum
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Desymmetrized Fock matrix','(10f12.6)',Fdes,nBasTot,nBasTot)
#endif
!----------------------------------------------------------------------*
! Symmetrize Fock matrix.                                              *
!----------------------------------------------------------------------*
call mma_allocate(Fsym,n2Full)
call mma_allocate(Aux1,n2Full)
call DGEMM_('N','N',nBasTot,nBasTot,nBasTot,One,Fdes,nBasTot,SmTr,nBasTot,Zero,Aux1,nBasTot)
call DGEMM_('T','N',nBasTot,nBasTot,nBasTot,One,SmTr,nBasTot,Aux1,nBasTot,Zero,Fsym,nBasTot)
#ifdef _DEBUGPRINT_
call RecPrt('Symmetrized overlap matrix','(10f12.6)',Fsym,nBasTot,nBasTot)
#endif
call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Compact the symmetrized Fock matrix.                                 *
!----------------------------------------------------------------------*
call mma_allocate(Fock,nTriTot)
ipTmp1 = 1
ipTmp2 = 1
do iSym=1,nSym
  indx = 1
  do iBas=1,nBas(iSym)
    do jBas=1,iBas
      Fock(indx) = Fsym(iBas+nBasTot*(jBas-1))
      indx = indx+1
    end do
  end do
# ifdef _DEBUGPRINT_
  write(Line,'(a,i2)') 'Fock matrix, symmetry ',iSym
  call TriPrt(Line,'(10f12.6)',Fock(ipTmp2),nBas(iSym))
# endif
  ipTmp1 = ipTmp1+nBas(iSym)*nBasTot+nBas(iSym)
  ipTmp2 = ipTmp2+nTri_Elem(nBas(iSym))
end do
!----------------------------------------------------------------------*
! Make ON basis                                                        *
!----------------------------------------------------------------------*
call mma_allocate(CMOs,nSqrTot)
call goLowdin(CMOs)
!----------------------------------------------------------------------*
! Transform Fock matrix to MO basis and diagonalize                    *
!----------------------------------------------------------------------*
call mma_allocate(Evec,nTriTot)
call mma_allocate(Fmo,nTriTot)
call mma_allocate(Aux1,nBasMax*nBasMax)
call mma_allocate(Aux2,nBasMax*nBasMax)

ipTmp1 = 1
ipTmp2 = 1
ipTmp3 = 1
ipTmp4 = 1
do iSym=1,nSym
  if (nBas(iSym) > 0) then
    call Square(Fock(ipTmp1),Aux1,1,nBas(iSym),nBas(iSym))
    call DGEMM_('N','N',nBas(iSym),nBas(iSym),nbas(iSym),One,Aux1,nBas(iSym),CMOs(ipTmp2),nBas(iSym),Zero,Aux2,nBas(iSym))
    call DGEMM_Tri('T','N',nBas(iSym),nBas(iSym),nBas(iSym),One,CMOs(ipTmp2),nBas(iSym),Aux2,nBas(iSym),Zero,Fmo(ipTmp3),nBas(iSym))
#   ifdef _DEBUGPRINT_
    write(Line,'(a,i2)') 'MO Fock matrix, symmetry ',iSym
    call TriPrt(Line,'(10f12.6)',Fmo(ipTmp3),nBas(iSym))
#   endif
    !call Jacob(Fmo(ipTmp3),CMOs(ipTmp2),nBas(iSym),nBas(iSym))
    call NIdiag(Fmo(ipTmp3),CMOs(ipTmp2),nBas(iSym),nBas(iSym))
#   ifdef _DEBUGPRINT_
    write(Line,'(a,i2)') 'Diagonal Fock matrix, symmetry ',iSym
    call TriPrt(Line,'(10f12.6)',Fmo(ipTmp3),nBas(iSym))
#   endif
    call goPickup(Fmo(ipTmp3),Evec(ipTmp4),nBas(iSym))
    call goSort(Evec(ipTmp4),CMOs(ipTmp2),nBas(iSym),nBas(iSym))
  end if
  ipTmp1 = ipTmp1+nTri_Elem(nBas(iSym))
  ipTmp2 = ipTmp2+nBas(iSym)*nBas(iSym)
  ipTmp3 = ipTmp3+nTri_Elem(nBas(iSym))
  ipTmp4 = ipTmp4+nBas(iSym)
end do
call mma_deallocate(Fmo)
call mma_deallocate(Aux2)
call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Present data.                                                        *
!----------------------------------------------------------------------*
if (PrintMOs) call PriMO('Start orbitals',.false.,.true.,Zero,1.0e6_wp,nSym,nBas,nBas,Label,Evec,Evec,CMOs,3)
call put_darray('Guessorb',CMOs,nSqrTot)
call put_darray('Guessorb energies',Evec,nBasTot)
call Put_iArray('nOrb',nBas,nSym)
call mma_allocate(Aux1,nBasTot)
Aux1(:) = 0.0d0
Lu = 20
Title = 'Guess orbitals'
call WrVec('GSSORB',Lu,'COE',nSym,nBas,nBas,CMOs,Aux1,Evec,iDummy,Title)
call mma_deallocate(Aux1)
!----------------------------------------------------------------------*
! Done, deallocate the rest.                                           *
!----------------------------------------------------------------------*
call mma_deallocate(Evec)
call mma_deallocate(CMOs)
call mma_deallocate(Fock)
call mma_deallocate(Fsym)
call mma_deallocate(Fdes)
call mma_deallocate(Sdes)
call mma_deallocate(Smat)
call mma_deallocate(Edes)
call mma_deallocate(Esym)
call mma_deallocate(DeTr)
call mma_deallocate(SmTr)

end subroutine Fmod1s
