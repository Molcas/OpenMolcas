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
! Copyright (C) 2004, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine BMtrx_Internal(nsAtom,nDimBC,nIter,mAtoms,iIter,mTR,TRVec,iTabAI,iTabAtoms,iTabBonds,nBonds,nMax,iRef,nQQ,nWndw)
!***********************************************************************
!                                                                      *
!     Objective: to handle curvilinear internal coordinates.           *
!                                                                      *
!                                                                      *
!     Authors: R. Lindh, Dept. of Theoretical Chemistry                *
!              University of Lund, SWEDEN.                             *
!              2004                                                    *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Slapaf_Info, only: Analytic_Hessian, BM, BMx, BSet, Cx, dBM, Degen, dqInt, dqInt_Aux, Gx, Gx0, HSet, HWRS, iBM, idBM, iOptC, &
                       KtB, lOld, MaxItr, mB_Tot, mdB_Tot, mq, NAC, nqBM, Numerical, PrQ, qInt, Smmtrc
use Kriging_Mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Five, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nsAtom, nDimBC, nIter, mAtoms, iIter, mTR, iTabAI(2,mAtoms), nMax, iTabAtoms(2,0:nMax,nsAtom), &
                                 nBonds, iTabBonds(3,nBonds), iRef, nWndw
real(kind=wp), intent(in) :: TRVec(nDimBC,mTR)
integer(kind=iwp), intent(out) :: nQQ
#include "warnings.h"
#include "print.fh"
integer(kind=iwp) :: i, i_Dim, iAtom, iB, iDum(6), iEnd, iGhi, iGlow, iPrint, iq, iqA, iQD, iqO, iQQ, iqR, iqRF, iqT, iRout, iSt, &
                     iX, ixyz, jIter, LuIC, M, mIter, N, nB, nB_Tot, ndB_Tot, nK, nq, nqA, nqB, nqO, nqRF, nqT, NRHS, nX
real(kind=wp) :: Alpha, Dum(1), Thr_ElRed, Thr_raw, Thr_small
logical(kind=iwp) :: g12K = .false., Proc, Proc_dB, Proc_H
character(len=32) :: filnam
character(len=14) :: cDum
integer(kind=iwp), allocatable :: Ind(:,:)
real(kind=wp), allocatable :: Degen2(:), EVal(:), F_c(:), G(:), GRef(:), GxR(:), K(:), KtBt(:,:), KtBu(:), KtM(:,:), Mult(:), &
                              Proj(:), qVal(:,:), Temp2(:,:)
character(len=14), allocatable :: qLbl(:)
integer(kind=iwp), external :: isfreeunit

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iPrint = 99
#else
iRout = 128
iPrint = nPrint(iRout)+1
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Proj,nDimBC,Label='Proj')
!                                                                      *
!***********************************************************************
!                                                                      *
Thr_raw = 3.0e-2_wp
if (HWRS) then
  Thr_ElRed = Thr_raw**2
else
  Thr_ElRed = Thr_raw
end if

i = 0
do iX=1,3*nsAtom
  iAtom = (iX+2)/3
  ixyz = iX-(iAtom-1)*3
  if (Smmtrc(ixyz,iAtom)) then
    i = i+1
    Proj(i) = One/Degen(ixyz,iAtom)
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! First some words about the handling of the symmetry. In this code
! I have selected to use only the unique components of the vectors.
! Hence, some care has to be taken to account for the degeneracy of
! some of the elements. This is handled by inserting the u matrix
! at appropriate places.

nQQ = nDimBC-mTR
if (allocated(qInt)) then
  if (size(qInt,1) /= nQQ) then
    call mma_deallocate(qInt)
    call mma_deallocate(dqInt)
  end if
end if
if (.not. allocated(qInt)) then
  call mma_allocate(qInt,nQQ,MaxItr,Label=' qInt')
  call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
  qInt(:,:) = Zero
  dqInt(:,:) = Zero
end if
if (allocated(dqInt_Aux)) then
  if (size(dqInt_Aux,1) /= nQQ) call mma_deallocate(dqInt_Aux)
end if
if ((.not. allocated(dqInt_Aux)) .and. (nSet > 1)) then
  call mma_allocate(dqInt_Aux,nQQ,MaxItr,nSet-1,Label='dqInt_Aux')
  dqInt_Aux(:,:,:) = Zero
end if

call mma_allocate(Degen2,nDimBC)
i = 0
do ix=1,3*nsAtom
  iAtom = (ix+2)/3
  ixyz = ix-(iAtom-1)*3
  if (Smmtrc(ixyz,iAtom)) then
    i = i+1
    Degen2(i) = Degen(ixyz,iAtom)
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Degen2',' ',Degen2,nDimBC,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
filnam = 'INTCOR'
LuIC = isfreeunit(4)
call molcas_open(LuIC,filnam)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! nq: Number of redundant internal coordinates.

Proc = .false.     ! Flag for processing B
Proc_dB = .false.  ! Flag for processing dB

Thr_small = 30.0_wp*deg2rad
do while (Thr_small > 1.0e-6_wp)
  call Get_Curvil(nq,nqRF,nqB,nqA,nqT,nqO,nsAtom,iIter,nIter,Cx,Proc,Dum,1,cDum,iRef,Dum,Dum,LuIC,iDum,iIter,Dum,1,1,Proc_dB, &
                  iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,mB_Tot,mdB_Tot,Dum,Dum,iDum,iDum,1,1,iDum,Thr_small)

  if (nq >= nQQ) exit
  Thr_small = Thr_small-Five*deg2rad
end do

rewind(LuIC)

if (nq == 0) then
  call WarningMessage(2,' BMtrx_internal: nq == 0')
  call Quit(_RC_INTERNAL_ERROR_)
end if

iGlow = nqRF+nqB+1
iGhi = nqRF+nqB+nqA
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'nq, nqB, nqA, nqT, nqO=',nq,nqB,nqA,nqT,nqO
#endif

! Now allocate some arrays which depend on nq

if (allocated(BM)) call mma_deallocate(BM)
if (allocated(iBM)) call mma_deallocate(iBM)
if (allocated(nqBM)) call mma_deallocate(nqBM)
call mma_allocate(BM,mB_Tot,Label='BM')
call mma_allocate(iBM,mB_Tot,Label='iBM')
call mma_allocate(nqBM,nq,Label='nqBM')
mq = nq

call mma_allocate(qVal,nq,nIter,Label='qVal')
call mma_allocate(qLbl,nq,Label='qLbl')
call mma_allocate(F_c,nq,Label='F_c')
call mma_allocate(Mult,nq,Label='Mult')
call mma_allocate(Ind,3,nq,Label='Ind')
call mma_allocate(GRef,9*nqA*nIter,Label='GRef')

F_c(:) = Zero
Mult(:) = Zero
qVal(:,:) = Zero
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Process the redundant internal coordinates for the current structure.

Proc = .true.
Proc_dB = .false.

call Get_Curvil(iq,iqRF,iqR,iqA,iqT,iqO,nsAtom,iIter,nIter,Cx,Proc,qVal,nq,qLbl,iRef,F_c,Mult,LuIC,Ind,iIter,GRef,iGlow,iGHi, &
                Proc_dB,iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,nB_Tot,ndB_Tot,BM,Dum,iBM,iDum,mB_Tot,mdB_Tot,nqBM,Thr_small)
rewind(LuIC)

if (iq /= nq) then
  call WarningMessage(2,' Error in BMtrx_internal')
  write(u6,*) 'In BMtrx_internal: iq /= nq'
  write(u6,*) 'iq=',iq
  write(u6,*) 'nq=',nq
  call Abend()
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (iPrint >= 49) then
  write(u6,*) 'nq, nqB, nqA, nqT, nqO=',nq,nqB,nqA,nqT,nqO
  call RecPrt('q-values',' ',qVal,nq,nIter)
  call RecPrt('Force Constant matrix in redundant basis',' ',F_c,1,nq)
  call RecPrt('Multiplicity factors',' ',Mult,1,nq)
  call RecPrt('Cx',' ',Cx,3*nsAtom,nIter)
  call RecPrt('Gx',' ',Gx,3*nsAtom,nIter)
end if
#endif

! Notation:
! X: cartesian coordinates
! q: redundant internal coordinates
! Q: nonredundant internal coordinates
! u: the degeneracy matrix
!
! dq = B u dx

! Start processing the B matrix for the current structure to
! generate the K matrix (dQ/dq).

call mma_allocate(K,nq**2,Label='K')
call mma_allocate(KtBu,nDimBC**2,Label='KtBu')
!                                                                      *
!***********************************************************************
!                                                                      *
! Stick in the metric of the force constants in the redundant
! space.

if (HWRS) then
  ! Scale each coordinate with the force constant

  i = 1
  do iq=1,nq
    nB = nqBM(iq)
    BM(i:i+nB-1) = F_c(iq)*BM(i:i+nB-1)
    i = i+nB
  end do
# ifdef _DEBUGPRINT_
  if (iPrint >= 99) then
    i = 1
    do iq=1,nq
      nB = nqBM(iq)
      call RecPrt('fcB',' ',BM(i),1,nB)
      i = i+nB
    end do
  end if
# endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Eliminate redundancy to produce the nonredundant internal coordinates.
!
!       t        t
! dQ = K M dq = K M B u dx

i = 1
do iq=1,nq
  nB = nqBM(iq)
  BM(i:i+nB-1) = Mult(iq)*BM(i:i+nB-1)
  i = i+nB
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (BSet) then
  call mma_allocate(G,nq*nq,Label='G')
  call mma_allocate(EVal,nTri_Elem(nq),Label='EVal')
  call ElRed2(nq,nDimBC,G,EVal,K,nK,Proj,g12K,Thr_ElRed,BM,iBM,mB_Tot,nqBM)

  if (nK > nQQ) call Remove_TR(nq,nDimBC,nQQ,K,nK,TRVec,mTR,BM,iBM,nqBM,mB_Tot)
  call mma_deallocate(EVal)
  call mma_deallocate(G)

  call Put_dArray('K',K,nq*nQQ)
else
  call Get_dArray('K',K,nq*nQQ)
end if
#ifdef _DEBUGPRINT_
if (iPrint >= 99) call RecPrt('K',' ',K,nq,nQQ)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Print the nonreduntant coordinates.

if (PrQ) then
  call mma_allocate(Temp2,nq,nQQ,Label='Temp2')
  Temp2(:,:) = reshape(K(1:nq*nQQ),[nq,nQQ])
  do iq=1,nq
    Temp2(iq,:) = Temp2(iq,:)/Mult(iq)
  end do
  call PrintQ(Temp2,qLbl,nq,nQQ,LuIC,Mult)
  call mma_deallocate(Temp2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (PrQ .and. (iPrint >= 6)) then
  write(u6,*)
  write(u6,'(A)') ' ******************************************'
  write(u6,'(A)') ' * Statistics of the internal coordinates *'
  write(u6,'(A)') ' ******************************************'
  write(u6,'(A,I5)') ' Translations and Rotations:     ',nqRF
  write(u6,'(A,I5)') ' Bonds                     :     ',nqB
  write(u6,'(A,I5)') ' Angles                    :     ',nqA
  write(u6,'(A,I5)') ' Torsions                  :     ',nqT
  write(u6,'(A,I5)') ' Out-of-plane angles       :     ',nqO
  write(u6,*)
end if
if (nq < nQQ) then
  call WarningMessage(2,' Error in BMtrx_internal')
  write(u6,*) 'In BMtrx_internal: nq < nQQ'
  write(u6,*) 'nq=',nq
  write(u6,*) 'nQQ=',nQQ
  call Abend()
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     The nonredundant internal coordinates are now defined!           *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Compute the values and the gradients of the nonredundant internal
! curvilinear coordinates for all the iterations.
!
! Loop over
! 1) the last to the first iteration if normal execusion
! 2) only the last if in transformation from internal to
!    cartesian

call mma_allocate(GxR,nDimBC,Label='GxR')

if (BSet) then
  iSt = nIter
  iEnd = iSt-min(nIter,nWndw+1)+1
else
  iEnd = nIter
  iSt = nIter
end if

! Note that the loop is in reverse order.

do jIter=iSt,iEnd,-1
  Proc_H = HSet .and. (jIter == iRef) .and. (.not. lOld)
  ! iOptC(8) = constrained optimization
  Proc_dB = Proc_H .and. (Analytic_Hessian .or. Numerical .or. btest(iOptC,8))
  ! Compute and store dBQQ in the reference structure
  if (Proc_dB) then
    if (allocated(dBM)) call mma_deallocate(dBM)
    if (allocated(idBM)) call mma_deallocate(idBM)
    call mma_allocate(dBM,mdB_Tot,Label='dBM')
    call mma_allocate(idBM,mdB_Tot*2,Label='idBM')
  else
    if (.not. allocated(dBM)) call mma_allocate(dBM,1,Label='dBM')
    if (.not. allocated(idBM)) call mma_allocate(idBM,1*2,Label='idBM')
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  call Get_Curvil(iq,iqRF,iqR,iqA,iqT,iqO,nsAtom,jIter,nIter,Cx,Proc,qVal,nq,qLbl,iRef,F_c,Mult,LuIC,Ind,iIter,GRef,iGlow,iGHi, &
                  Proc_dB,iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,nB_Tot,ndB_Tot,BM,dBM,iBM,idBM,mB_Tot,mdB_Tot,nqBM, &
                  Thr_small)
  rewind(LuIC)

  if (iq /= nq) then
    write(u6,*) 'In BMtrx_internal: iq /= nq'
    write(u6,*) 'iq=',iq
    write(u6,*) 'nq=',nq
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Form the gradients of Q
  !                                       T
  ! Form the K(t)MB matrix for dQ = K M B u dx

  i = 1
  do iq=1,nq
    nB = nqBM(iq)
    BM(i:i+nB-1) = Mult(iq)*BM(i:i+nB-1)
    i = i+nB
  end do
  KtBu(1:nQQ*nDimBC) = Zero
  do iQQ=0,nQQ-1
    i = 1
    do iq=1,nq
      nB = nqBM(iq)
      do iB=0,nB-1
        i_Dim = iBM(i)
        KtBu((i_Dim-1)*nQQ+iQQ+1) = KtBu((i_Dim-1)*nQQ+iQQ+1)+K(iQQ*nq+iq)*BM(i)
        i = i+1
      end do
    end do
  end do
# ifdef _DEBUGPRINT_
  if (iPrint >= 99) then
    call RecPrt(' The BM matrix',' ',BM(:),1,size(BM))
    call RecPrt(' The K matrix',' ',K,nq,nQQ)
    call RecPrt(' The K(t)B matrix',' ',KtBu,nQQ,nDimBC)
  end if
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Store the B matrix for the last structure

  if (jIter == nIter) then
    nX = 3*nsAtom
    call mma_allocate(BMx,nX,nX,Label='BMx')
    BMx(:,:) = Zero

    ! modify from compact to full Cartesian storage.

    i_Dim = 0
    do iX=1,nX
      iAtom = (iX+2)/3
      ixyz = iX-(iAtom-1)*3
      if (Smmtrc(ixyz,iAtom)) then
        i_Dim = i_Dim+1
        do iQQ=1,nQQ
          iQD = (i_Dim-1)*nQQ+iQQ
          BMx(iX,iQQ) = KtBu(iQD)
        end do
      else
        BMx(iX,1:nQQ) = Zero
      end if
    end do

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Branch out if only values are to be computed.

  if (BSet .or. Numerical) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Form the gradient for iteration jIter for the new definition
    ! of the K matrix.
    !
    ! dq/dx dE/dq = dE/dq
    !
    ! The B-matrix is stored (3*natom x nQQ)
    ! KtBu is stored nQQ, nDimBC

    call mma_allocate(KtBt,nDimBC,nQQ,Label='KtBt')
    call TRNSPS(nQQ,nDimBC,KtBu,KtBt)

    ! Strip KtB of the degeneracy factor (full).

    do iQQ=1,nQQ
      KtBt(1:nDimBC,iQQ) = KtBt(1:nDimBC,iQQ)/Degen2(1:nDimBC)
    end do

    M = nDimBC
    N = nQQ
    NRHS = 1

    call NRed(Gx(:,:,jIter),GxR,3*nsAtom,nDimBC,Smmtrc)
    call Eq_Solver('N',M,N,NRHS,KtBt,.false.,Degen2,GxR,dqInt(:,jIter))
    if (nSet > 1) then
      call NRed(Gx0(:,:,jIter),GxR,3*nsAtom,nDimBC,Smmtrc)
      call Eq_Solver('N',M,N,NRHS,KtBt,.false.,Degen2,GxR,dqInt_Aux(:,jIter,1))
    end if
    if (nSet > 2) then
      call NRed(NAC(:,:,jIter),GxR,3*nsAtom,nDimBC,Smmtrc)
      call Eq_Solver('N',M,N,NRHS,KtBt,.false.,Degen2,GxR,dqInt_Aux(:,jIter,2))
    end if
    !call RecPrt('GxR   ',' ',GxR,nDimBC,1)
    !call RecPrt('KtB   ',' ',KtBt,nQQ,nDimBC)
    !call RecPrt('drInt ',' ',dqInt(:,jIter),nQQ,1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Save pointer to KtB to be used in backtransforming the
    ! Hessian to internal coordinates.

    if (Proc_H) then
      call mma_allocate(KtB,nDimBC,nQQ,Label='KtB')
      KtB(:,:) = KtBt(:,:)
      !call RecPrt('KtB',' ',KtB,nDimBC,nQQ)
    end if
    call mma_deallocate(KtBt)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the values
!      t
! Q = K M q

jIter = 1
mIter = nIter
if (.not. BSet) then
  jIter = nIter
  mIter = 1
end if
call mma_allocate(KtM,nQQ,nq,Label='KtM')
do iq=1,nq
  Alpha = Mult(iq)
  do iQQ=1,nQQ
    KtM(iQQ,iq) = K(iq+(iQQ-1)*nq)*Alpha
  end do
end do
call DGEMM_('N','N',nQQ,mIter,nq,One,KtM,nQQ,qVal(:,jIter),nq,Zero,qInt(:,jIter),nQQ)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (iPrint >= 49) then
  call RecPrt(' The K Matrix',' ',K,nq,nQQ)
  call RecPrt(' q-values',' ',qVal,nq,nIter)
  call RecPrt('Q-values',' ',qInt,nQQ,nIter)
  call RecPrt('Cx',' ',Cx,3*nsAtom,nIter)
end if
if (BSet .and. (iPrint >= 49)) then
  call RecPrt('Q-gradients',' ',dqInt,nQQ,nIter)
  call RecPrt('Gx',' ',Gx,3*nsAtom,nIter)
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call mma_deallocate(KtM)
if (allocated(GxR)) call mma_deallocate(GxR)
call mma_deallocate(KtBu)
call mma_deallocate(K)
call mma_deallocate(GRef)
call mma_deallocate(Ind)
call mma_deallocate(Mult)
call mma_deallocate(F_c)
call mma_deallocate(qLbl)
call mma_deallocate(qVal)

call mma_deallocate(Degen2)
call mma_deallocate(Proj)

close(LuIC)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine BMtrx_Internal
