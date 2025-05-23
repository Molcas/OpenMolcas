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
! This routine creates start orbitals from model fock matrix elements  *
! generated by SEWARD.                                                 *
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
subroutine FckByInt(iReturncode,StandAlone)

use Index_Functions, only: nTri_Elem
use GuessOrb_Global, only: GapThr, iPrFmt, Label, nBas, nDel, nSym, PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr
#ifdef _HDF5_
use GuessOrb_Global, only: wfn_energy, wfn_mocoef, wfn_occnum, wfn_orbene, wfn_tpidx
use mh5, only: mh5_put_dset
#endif
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturncode
logical(kind=iwp), intent(in) :: StandAlone
integer(kind=iwp) :: i, i1, iBas, iComp, ijL, ijS, ijT, ik, inCMO, IndType(7,8), inEps, inFck, inOvl, iOff, iOpt, ipCOk, ipEE, &
                     ipEE0, ipOk, ipOk0, ipOkk, irc, iSym, iSymlb, j1, jBas, jk, jOff, k, kBas, kOff, kSpin, Lu, nActEl, nAsh(8), &
                     nB, nBasMax, nBasTot, nC, nD, nIsh(8), nOkk, nOrb(8), nS, nSqrTot, nTriTot
real(kind=wp) :: dActEl, ei, ej, Enr_go, tmp, tmp1, tmp2, xocc
character(len=180) :: Line
character(len=80) :: Title
character(len=8) :: Lbl
real(kind=wp), allocatable :: CMO(:), Eps(:), Fck(:), Ovl(:), T1(:), T2(:), T3(:)
#ifdef _HDF5_
integer(kind=iwp) :: IndTypeT(8,7)
character, allocatable :: typestring(:)
#endif

!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
iReturncode = 0
call getenvf('MOLCAS_TEST',Line)
!----------------------------------------------------------------------*
! Do some counting                                                     *
!----------------------------------------------------------------------*
nBasTot = 0
nBasMax = 0
nTriTot = 0
nSqrTot = 0
do iSym=1,nSym
  nBasTot = nBasTot+nbas(iSym)
  nBasMax = max(nbasmax,nBas(iSym))
  nTriTot = nTriTot+nTri_Elem(nBas(iSym))
  nSqrTot = nSqrTot+nBas(iSym)*nBas(iSym)
end do
!----------------------------------------------------------------------*
! Get model Fock matrix.                                               *
!----------------------------------------------------------------------*
inFck = nTriTot+6
call mma_allocate(Fck,inFck)
iRc = -1
iSymlb = 1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Lbl = 'FckInt'
iComp = 1
call RdOne(irc,iOpt,Lbl,iComp,Fck,iSymlb)
if (iRc /= 0) then
  iReturncode = 1
  call mma_deallocate(Fck)
  write(u6,*) '***'
  write(u6,*) '*** WARNING:'
  write(u6,*) '*** Guessorb did not produce start orbitals!!!'
  write(u6,*) '***'
  return
end if
#ifdef _DEBUGPRINT_
i = 1
do iSym=1,nSym
  call TriPrt('FckInt','(12f12.6)',Fck(i),nBas(iSym))
  !call NrmClc(Fck(i),nTri_Elem(nBas(iSym)),'FckbyInt','Fck(i)')
  i = i+nTri_Elem(nBas(iSym))
end do
#endif
!----------------------------------------------------------------------*
! Make symmetric orthonormal orbital basis.                            *
!----------------------------------------------------------------------*
inCMO = nSqrTot
call mma_allocate(CMO,inCMO)
call goLowdin(CMO)
#ifdef _DEBUGPRINT_
i = 1
do iSym=1,nSym
  nB = nBas(iSym)
  call RecPrt('CMO','(12f12.6)',CMO(i),nB,nB)
  !call NrmClC(CMO(i),nB**2,'FckbyInt','CMO(i)')
  i = i+nB*nB
end do
#endif
!----------------------------------------------------------------------*
! Get overlap matrix                                                   *
!----------------------------------------------------------------------*
inOvl = nTriTot+6
call mma_allocate(Ovl,inOvl)
iSymlb = 1
Lbl = 'Mltpl  0'
call RdOne(irc,iOpt,Lbl,iComp,Ovl,iSymlb)
#ifdef _DEBUGPRINT_
i = 1
do iSym=1,nSym
  call TriPrt('Ovlp','(12f12.6)',Ovl(i),nBas(iSym))
  !call NrmClc(Ovl(i),nTri_Elem(nBas(iSym)),'FckbyInt','Ovl(i)')
  i = i+nTri_Elem(nBas(iSym))
end do
#endif
!----------------------------------------------------------------------*
! Transform: F = S eps S                                               *
!----------------------------------------------------------------------*
call mma_allocate(T1,nBasMax**2)
call mma_allocate(T2,nBasMax**2)
call mma_allocate(T3,nBasMax**2)
ijT = 1
ijS = 1
ijL = 1
do iSym=1,nSym
  nB = nBas(iSym)
  if (nB > 0) then
    call Square(Fck(ijT),T1,1,nB,nB)
    call Square(Ovl(ijT),T2,1,nB,nB)
    call DGEMM_('N','N',nB,nB,nB,One,T1,nB,T2,nB,Zero,T3,nB)
    call DGEMM_Tri('T','N',nB,nB,nB,One,T2,nB,T3,nB,Zero,Fck(ijT),nB)
#   ifdef _DEBUGPRINT_
    call TriPrt('Fock matrix with metric','(12f12.6)',Fck(ijT),nB)
    !call NrmClc(Fck(ijT),nTri_Elem(nB),'FckbyInt','Fck(ijT)')
#   endif
  end if
  ijT = ijT+nTri_Elem(nB)
  ijS = ijS+nB*nB
  ijL = ijL+nB
end do
call mma_deallocate(T3)
call mma_deallocate(T2)
call mma_deallocate(T1)
!----------------------------------------------------------------------*
! Diagonalize the model Fock matrix                                    *
!----------------------------------------------------------------------*
inEps = nBasTot
call mma_allocate(Eps,inEps)
call mma_allocate(T1,nBasMax**2)
call mma_allocate(T2,nBasMax**2)
call mma_allocate(T3,nBasMax**2)
ijT = 1
ijS = 1
ijL = 1
do iSym=1,nSym
  nB = nBas(iSym)
  nS = nBas(iSym)-nDel(iSym)
  if (nB > 0) then
    call Square(Fck(ijT),T1,1,nB,nB)
    call DGEMM_('N','N',nB,nS,nB,One,T1,nB,CMO(ijS),nB,Zero,T2,nB)
    call DGEMM_Tri('T','N',nS,nS,nB,One,CMO(ijS),nB,T2,nB,Zero,T3,nS)
#   ifdef _DEBUGPRINT_
    call TriPrt('Transformed Fock matrix','(12f12.6)',T3,nB)
    !call NrmClc(T3,nTri_Elem(nB),'FckbyInt','Transformed Fck')
#   endif
    call NIdiag(T3,CMO(ijS),nS,nB)
    call goPickup(T3,Eps(ijL),nS)
    call goSort(Eps(ijL),CMO(ijS),nS,nB)

    do i=1,nS
      call VecPhase(CMO(ijS+(i-1)*nB),nB)
    end do
  end if
  ijT = ijT+nTri_Elem(nB)
  ijS = ijS+nB*nB
  ijL = ijL+nB
end do
#ifdef _DEBUGPRINT_
i = 1
do iSym=1,nSym
  nB = nBas(iSym)
  call RecPrt('CMO','(12f12.6)',CMO(i),nB,nB)
  !call NrmClC(CMO(i),nB**2,'FckbyInt','CMO(i)')
  i = i+nB*nB
end do
#endif
call mma_deallocate(T3)
call mma_deallocate(T2)
call mma_deallocate(T1)
!----------------------------------------------------------------------*
! Diagonalize T in virtual space.                                      *
!----------------------------------------------------------------------*
iRc = -1
iSymlb = 1
Lbl = 'Kinetic'
call RdOne(irc,iOpt,Lbl,iComp,Fck,iSymlb)
if (iRc == 0) then
  call mma_allocate(T1,nBasMax**2)
  call mma_allocate(T2,nBasMax**2)
  call mma_allocate(T3,nBasMax**2)
  ijT = 1
  ijS = 1
  ijL = 1
  do iSym=1,nSym
    nB = nBas(iSym)
    nD = nDel(iSym)
    nC = 0
    do iBas=1,nB-nD
      if (Eps(ijL+iBas-1) < -1.0e-3_wp) nC = nC+1
    end do
    nS = nB-nC-nD
    if (nS > 0) then

      ! Generate standardized virtual orbitals before we proceed.
      ! The virtual orbitals generated previously are not well
      ! defined and might differ substantially with different
      ! hardware/software and compiler options. To be able to
      ! compare we will need these standardized virtual orbitals.
      ! In real production calculations this step could for all
      ! practical purposes be skipped.

      call Virt_Space(CMO(ijS),CMO(ijS+nB*nC),Ovl(ijT),nB,nC,nS)

      call Square(Fck(ijT),T1,1,nB,nB)
      call DGEMM_('N','N',nB,nS,nB,One,T1,nB,CMO(ijS+nB*nC),nB,Zero,T2,nB)

      call DGEMM_Tri('T','N',nS,nS,nB,One,CMO(ijS+nB*nC),nB,T2,nB,Zero,T3,nS)
#     ifdef _DEBUGPRINT_
      call TriPrt('Virtual space','(12f12.6)',T3,nS)
#     endif
      call NIdiag(T3,CMO(ijS+nB*nC),nS,nB)
      call goPickup(T3,Eps(ijL+nC),nS)
      call goSort(Eps(ijL+nC),CMO(ijS+nB*nC),nS,nB)
#     ifdef _DEBUGPRINT_
      call RecPrt('Eps',' ',Eps(ijL+nC),nS,1)
      call RecPrt('Virtual Orbitals',' ',CMO(ijS+nB*nC),nB,nS)
#     endif

      ! Now order degenerate orbitals. This is only important for
      ! verification runs.

      do iBas=nC+1,nB-nD-1
        ei = Eps(ijL+iBas-1)
        tmp1 = Zero
        do kBas=1,nB
          ik = ijS+(iBas-1)*nB+kBas-1
          tmp1 = tmp1+abs(CMO(ik)*real(kBas,kind=wp))
        end do
        do jBas=iBas+1,nB-nD
          ej = Eps(ijL+jBas-1)
          if (abs(ei-ej) < 1.0e-12_wp) then
            tmp2 = Zero
            do kBas=1,nB
              jk = ijS+(jBas-1)*nB+kBas-1
              tmp2 = tmp2+abs(CMO(jk)*real(kBas,kind=wp))
            end do
            if (tmp2 > tmp1) then
              tmp = tmp2
              tmp2 = tmp1
              tmp1 = tmp
              Eps(ijL+iBas-1) = ej
              Eps(ijL+jBas-1) = ei
              ei = ej
              i1 = ijS+(iBas-1)*nB
              j1 = ijS+(jBas-1)*nB
              call DSwap_(nB,CMO(i1),1,CMO(j1),1)
            end if
          end if

        end do
      end do

      ! Introduce "standard" phase.

      do iBas=1,nB
        call VecPhase(CMO(ijS+(iBas-1)*nB),nB)
      end do

#     ifdef _DEBUGPRINT_
      call RecPrt('Eps',' ',Eps(ijL+nC),nS,1)
      call RecPrt('Virtual Orbitals',' ',CMO(ijS+nB*nC),nB,nS)
#     endif
      Eps(ijL+nC:ijL+nB-nD-1) = Eps(ijL+nC:ijL+nB-nD-1)+Three
      Eps(ijL+nB-nD:ijL+nB-1) = 999.0_wp
      do iBas=1,nB-nD
        if (Eps(ijL+iBas-1) > TThr) nDel(iSym) = nDel(iSym)+1
      end do
    end if
    ijT = ijT+nTri_Elem(nB)
    ijS = ijS+nB*nB
    ijL = ijL+nB
  end do
  call mma_deallocate(T3)
  call mma_deallocate(T2)
  call mma_deallocate(T1)
  !--------------------------------------------------------------------*
  ! Print orbital space data.                                          *
  !--------------------------------------------------------------------*
  if (StandAlone) then
    write(u6,'(a,es10.3)') 'Threshold for linear dependence due to S:',SThr
    write(u6,'(a,es10.3)') 'Threshold for linear dependence due to T:',TThr
    write(u6,*)
    write(u6,'(a,8i5)') 'Total number of basis functions',(nBas(iSym),iSym=1,nSym)
    write(u6,'(a,8i5)') 'Deleted orbitals               ',(nDel(iSym),iSym=1,nSym)
    write(u6,*)
  end if
end if
!----------------------------------------------------------------------*
! Present data.                                                        *
!----------------------------------------------------------------------*
call mma_allocate(T1,nBasTot)
call mma_allocate(T2,nBasTot)
T1(:) = Zero
call GoPop(Eps,T1,T2,nBasTot,PrintEor,PrThr,GapThr)
iBas = 0
dActEl = Zero
IndType(1:5,1:nSym) = 0
do iSym=1,nSym
  IndType(6,iSym) = nBas(iSym)-nDel(iSym)
  IndType(7,iSym) = nDel(iSym)
  do kBas=1,nBas(iSym)-nDel(iSym)
    iBas = iBas+1
    if (T1(iBas) > 1.99_wp) then
      IndType(2,iSym) = IndType(2,iSym)+1
      IndType(6,iSym) = IndType(6,iSym)-1
    else if (T1(iBas) > 0.01_wp) then
      IndType(4,iSym) = IndType(4,iSym)+1
      IndType(6,iSym) = IndType(6,iSym)-1
      dActEl = dActEl+T1(iBas)
    end if
  end do
end do
nActEl = int(dActEl+Half)
if (PrintMOs) then
  call PriMO('Start orbitals (virtuals shifted)',.true.,.true.,Zero,PrThr,nSym,nBas,nBas,Label,Eps,T1,CMO,iPrFmt)
  call xflush(u6)
end if
if (PrintPop) call Charge(nSym,nBas,Label,CMO,T1,Ovl,2,.true.,.true.)
call put_darray('Guessorb',CMO,nSqrTot)
call put_darray('Guessorb energies',Eps,nBasTot)
nOrb(1:nSym) = nBas(1:nSym)-nDel(1:nSym)
call Put_iArray('nOrb',nOrb,nSym)
call Put_iArray('nDel_go',nDel,nSym)
call Put_iArray('nDel',nDel,nSym)
nIsh(1:nSym) = IndType(2,1:nSym)
call Put_iArray('nIsh',nIsh,nSym)
nAsh(1:nSym) = IndType(4,1:nSym)
call Put_iArray('nAsh',nAsh,nSym)
call Put_iScalar('nActel',nActEl)
kSpin = 1 ! always same alpha and beta orbs
call Put_iScalar('Multiplicity',kSpin)
Enr_go = Zero
ipEE0 = 1
ipOk0 = 1
do iSym=1,nSym
  do i=0,nIsh(iSym)+nAsh(iSym)-1
    ipEE = ipEE0+i
    ipOk = ipOk0+i
    Enr_go = Enr_go+T1(ipOk)*Eps(ipEE)
  end do
  ipEE0 = ipEE0+nBas(iSym)
  ipOk0 = ipOk0+nBas(iSym)
end do
call Put_dScalar('Last energy',Enr_go)
#ifdef _HDF5_
call mh5_put_dset(wfn_energy,Enr_go)
#endif
Lu = 20
Title = 'Guess orbitals'
call WrVec('GSSORB',Lu,'COEI',nSym,nBas,nBas,CMO,T1,Eps,IndType,Title)
#ifdef _HDF5_
IndTypeT(:,:) = transpose(IndType(:,:))
call mma_allocate(typestring,nBasTot)
call orb2tpstr(nSym,nBas,IndTypeT(:,1),IndTypeT(:,2),IndTypeT(:,3),IndTypeT(:,4),IndTypeT(:,5),IndTypeT(:,6),IndTypeT(:,7), &
               typestring)
call mh5_put_dset(wfn_tpidx,typestring)
call mma_deallocate(typestring)
call mh5_put_dset(wfn_mocoef,CMO)
call mh5_put_dset(wfn_occnum,T1)
call mh5_put_dset(wfn_orbene,Eps)
#endif

! Compute density matrix (re-use memory allocated in Ovl)
iOff = 1
jOff = 1
kOff = 1
do iSym=1,nSym
  ipOkk = iOff
  nOkk = nIsh(iSym)+nAsh(iSym)
  ipCOk = jOff
  if (nBas(iSym) > 0) then
    do k=0,nOkk-1
      xocc = sqrt(T1(k+ipOkk))
      CMO(ipCOk:ipCOk+nBas(iSym)-1) = xocc*CMO(ipCOk:ipCOk+nBas(iSym)-1)
      ipCOk = ipCOk+nBas(iSym)
    end do
    call DGEMM_Tri('N','T',nBas(iSym),nBas(iSym),nOkk,One,CMO(jOff),nBas(iSym),CMO(jOff),nBas(iSym),Zero,Ovl(kOff),nBas(iSym))
    iOff = iOff+nBas(iSym)
    jOff = jOff+nBas(iSym)**2
    kOff = kOff+nTri_Elem(nBas(iSym))
  end if
end do
call Fold_tMat(nSym,nBas,Ovl,Ovl)
call Put_dArray('D1ao',Ovl,nTriTot)

call mma_deallocate(T2)
call mma_deallocate(T1)
!----------------------------------------------------------------------*
! Done, deallocate the rest.                                           *
!----------------------------------------------------------------------*
call mma_deallocate(Eps)
call mma_deallocate(Ovl)
call mma_deallocate(CMO)
call mma_deallocate(Fck)

end subroutine FckByInt
