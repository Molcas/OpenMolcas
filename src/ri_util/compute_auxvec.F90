!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Compute_AuxVec(ipVk,ipZpk,myProc,nProc,ipUk)

use pso_stuff
use Basis_Info, only: nBas, nBas_Aux
use Gateway_global, only: force_out_of_core
use RICD_Info, only: Do_RI, Cholesky
use Symmetry_Info, only: nIrrep
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use ExTerm, only: iMP2prpt, DMLT

implicit real*8(a-h,o-z)
integer ipVk(nProc), ipZpk(nProc)
integer, optional :: ipUk(nProc)
#include "stdalloc.fh"
#include "real.fh"
#include "cholesky.fh"
#include "etwas.fh"
#include "exterm.fh"
type(DSBA_Type) DSQ, ChM(5), DLT2
logical DoExchange, DoCAS, Estimate, Update
integer nIOrb(0:7), nV_l(0:7), nV_t(0:7)
integer nU_l(0:7), nU_t(0:7)
integer ipTxy(0:7,0:7,2), jp_V_k
#include "chotime.fh"
character*8 Method
real*8, allocatable :: TmpD(:), Zv(:), Qv(:), Scr(:)

!                                                                      *
!***********************************************************************
!                                                                      *
DoExchange = Exfac /= Zero

nV_ls = 0
do i=0,nIrrep-1
  nV_l(i) = NumCho(i+1) ! local # of vecs in parallel run
  nV_ls = nV_ls+nV_l(i)
  nV_t(i) = nV_l(i)
end do
call GAIGOP(nV_t,nIrrep,'+') ! total # of vecs
if (nV_t(0) == 0) then
  call WarningMessage(2,'Compute_AuxVec: no total symmetric vectors!!')
  call Abend()
end if

if (iMp2prpt == 2) then
  if (.not. present(ipUk)) then
    write(6,*) 'No ipUk input present!'
    call Abend()
  end if
  nU_ls = 0
  do i=0,nIrrep-1
    nU_l(i) = NumCho(i+1) ! local # of vecs in parallel run
    nU_ls = nU_ls+nU_l(i)
    nU_t(i) = nU_l(i)
  end do
  call GAIGOP(nU_t,nIrrep,'+') ! total # of vecs
  if (nU_t(0) == 0) then
    call WarningMessage(2,'Compute_AuxVec: no total symmetric vectors!!')
    call Abend()
  end if
end if

NChVMx = 0
nQMax = 0
do i=0,nIrrep-1
  NChVMx = max(NChVMx,nV_t(i))
  nQMax = max(nQMax,nBas_Aux(i))
  nChOrb(i,1) = 0
  nChOrb(i,2) = 0
end do
nQvMax = nQMax*NChVMx
call mma_allocate(Scr,nQMax)

DoCAS = lPSO

if (nV_ls >= 1) then ! can be = 0 in a parallel run

  jp_V_k = ipVk(myProc)
  jp_Z_p_k = ipZpk(myProc)
  jp_U_k = 1
  if (iMp2prpt == 2) then
    jp_U_k = ipUk(myProc)
  end if
  !*********************************************************************
  !                                                                    *
  !      Get (and transform) the density matrices                      *
  !                                                                    *
  !*********************************************************************

  Timings = .false.
  !Timings = .true.

  call Get_iArray('nIsh',nIOrb,nIrrep)

  if (iMp2prpt /= 2) then
    if (DoCAS .and. lSA) then
      nSA = 5
      do i=1,nSA
        call Allocate_DT(DMLT(i),nBas,nBas,nSym,aCase='TRI')
        DMLT(i)%A0(:) = D0(:,i)
      end do
      ! Refold some density matrices
      do iIrrep=0,nIrrep-1
        ij = 0
        do iBas=1,nBas(iIrrep)
          do jBas=1,iBas-1
            ij = ij+1
            DMLT(1)%SB(iIrrep+1)%A1(ij) = Two*DMLT(1)%SB(iIrrep+1)%A1(ij)
            DMLT(3)%SB(iIrrep+1)%A1(ij) = Two*DMLT(3)%SB(iIrrep+1)%A1(ij)
            DMLT(5)%SB(iIrrep+1)%A1(ij) = Two*DMLT(5)%SB(iIrrep+1)%A1(ij)
          end do
          ij = ij+1
        end do
      end do
    else
      call Allocate_DT(DMLT(1),nBas,nBas,nSym,aCase='TRI')
      call Get_D1AO_Var(DMLT(1)%A0,nDens)
    end if
  else
    call Allocate_DT(DMLT(1),nBas,nBas,nSym,aCase='TRI')
    call Get_D1AO(DMLT(1)%A0,nDens)
  end if

  if (nKdens == 2) then
    call Allocate_DT(DMLT(2),nBas,nBas,nSym,aCase='TRI')
    ! spin-density matrix
    call Get_D1SAO_Var(DMLT(2)%A0,nDens)
    call daxpy_(nDens,-One,DMLT(1)%A0,1,DMLT(2)%A0,1)
    call dscal_(nDens,-Half,DMLT(2)%A0,1) ! beta DMAT
    call daxpy_(nDens,-One,DMLT(2)%A0,1,DMLT(1)%A0,1) ! alpha DMAT
  else if ((nKdens > 4) .or. (nKdens < 1)) then
    call WarningMessage(2,'Compute_AuxVec: invalid nKdens!!')
    call Abend()
  end if
  if (iMp2prpt == 2) then
    call Allocate_DT(DLT2,nBas,nBas,nSym,aCase='TRI')
    call Get_D1AO_Var(DLT2%A0,nDens)
    call daxpy_(nDens,-One,DMLT(1)%A0,1,DLT2%A0,1)
  else
  end if
  !*********************************************************************
  !                                                                    *
  !     Compute Fr+In+Ac localized orbitals                            *
  !     using Cholesky  decomposition for PD matrices                  *
  !     using Eigenvalue decomposition for non-PD matrices (SA-CASSCF) *
  !                                                                    *
  !*********************************************************************
  !DoExchange = Exfac /= Zero

  call Get_cArray('Relax Method',Method,8)
  if (Method == 'MCPDFT ') exfac = 1.0d0
  DoExchange = Exfac /= Zero

  if (DoExchange .or. DoCAS) then
    do i=1,nKdens
      call Allocate_DT(ChM(i),nBas,nBas,nSym)
    end do
    if (lSA) then
      call mma_allocate(TmpD,nDens,Label='TmpD')
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! PD matrices

    call Allocate_DT(DSQ,nBas,nBas,nSym)
    do j=1,nKvec
      if (lSA) then
        if (j == 1) then
          call dcopy_(nDens,DMLT(1)%A0,1,TmpD,1)
        else if (j == 2) then
          call dcopy_(nDens,DMLT(3)%A0,1,TmpD,1)
        end if
        call UnFold(TmpD,nDens,DSQ%A0,size(DSQ%A0),nIrrep,nBas)
      else
        call UnFold(DMLT(j)%A0,nDens,DSQ%A0,size(DSQ%A0),nIrrep,nBas)
      end if

      do i=0,nIrrep-1
        call CD_InCore(DSQ%SB(i+1)%A2,nBas(i),ChM(j)%SB(i+1)%A2,nBas(i),nChOrb(i,j),1.0d-12,irc)
      end do
      if (irc /= 0) then
        write(6,*) 'Compute_AuxVec: failed to get Cholesky MOs !'
        call Abend()
      end if
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! non-PD matrices

    if (lSA) then

      do i=3,4

        ! Get the appropriate density matrix

        if (i == 3) then
          call dcopy_(nDens,DMLT(2)%A0,1,TmpD,1)
        else if (i == 4) then
          call dcopy_(nDens,DMLT(4)%A0,1,TmpD,1)
        end if

        ! And eigenvalue-decompose it

        iOffDSQ = 0
        do isym=0,nIrrep-1

          ChM(i)%SB(iSym+1)%A2(:,:) = Zero
          call dcopy_(nbas(isym),[One],0,ChM(i)%SB(iSym+1)%A2,nBas(isym)+1)
          call NIdiag(TmpD(1+iOffDSQ),ChM(i)%SB(iSym+1)%A2,nBas(isym),nBas(isym))

          ! First sort eigenvectors and eigenvalues

          do j=1,nBas(isym)
            irun = iOffDSQ+j*(j+1)/2
            do k=j,nBas(isym)
              jrun = iOffDSQ+k*(k+1)/2
              if (TmpD(irun) < TmpD(jrun)) then
                tmp = TmpD(irun)
                TmpD(irun) = TmpD(jrun)
                TmpD(jrun) = tmp
                do l=1,nBas(isym)
                  tmp = ChM(i)%SB(iSym+1)%A2(l,j)
                  ChM(i)%SB(iSym+1)%A2(l,j) = ChM(i)%SB(iSym+1)%A2(l,k)
                  ChM(i)%SB(iSym+1)%A2(l,k) = tmp
                end do
              end if
            end do
          end do

          Cho_thrs = 1.0d-12

          l = 0
          do j=1,nBas(isym)
            if (TmpD(iOffDSQ+j*(j+1)/2) > Cho_thrs) then
              l = l+1
              tmp = sqrt(TmpD(iOffDSQ+j*(j+1)/2))
              do k=1,nBas(isym)
                ChM(i)%SB(iSym+1)%A2(k,l) = ChM(i)%SB(iSym+1)%A2(k,j)*tmp
              end do
            end if
          end do
          npos(isym,i-2) = l

          do j=1,nBas(isym)
            if (-TmpD(iOffDSQ+j*(j+1)/2) > Cho_thrs) then
              l = l+1
              irun = (l-1)*nBas(isym)
              jrun = (j-1)*nBas(isym)
              tmp = sqrt(-TmpD(iOffDSQ+j*(j+1)/2))
              do k=1,nBas(isym)
                ChM(i)%SB(iSym+1)%A2(k,l) = ChM(i)%SB(iSym+1)%A2(k,j)*tmp
              end do
            end if
          end do
          nChOrb(isym,i) = l

          iOffDSQ = iOffDSQ+nBas(isym)*(nBas(isym)+1)/2
        end do

      end do

      ! Refold the other DM

      do iIrrep=0,nIrrep-1
        ij = 0
        do iBas=1,nBas(iIrrep)
          do jBas=1,iBas-1
            ij = ij+1
            DMLT(2)%SB(iIrrep+1)%A1(ij) = Two*DMLT(2)%SB(iIrrep+1)%A1(ij)
            DMLT(4)%SB(iIrrep+1)%A1(ij) = Two*DMLT(4)%SB(iIrrep+1)%A1(ij)
          end do
          ij = ij+1
        end do
      end do
    end if
    call Deallocate_DT(DSQ)
    if (lSA) call mma_deallocate(TmpD)
  end if
  !*********************************************************************
  !                                                                    *
  !      First contract the RI vectors with the density matrix         *
  !                                                                    *
  !*********************************************************************
  ! Pointers to the Cholesky vectors of P2
  mAO = 0
  iOff = 0
  do kIrrep=0,nIrrep-1 ! compound symmetry
    iOff2 = 0
    do jIrrep=0,nIrrep-1
      iIrrep = ieor(jIrrep,kIrrep)
      if (iIrrep < jIrrep) then
        nnAorb = nASh(iIrrep)*nAsh(jIrrep)
      else if (iIrrep == jIrrep) then
        nnAorb = nAsh(iIrrep)*(nAsh(iIrrep)+1)/2
      else
        cycle
      end if
      ipTxy(iIrrep,jIrrep,1) = 1+iOff2+iOff
      ipTxy(jIrrep,iIrrep,1) = ipTxy(iIrrep,jIrrep,1)
      if (lSA) then
        ipTxy(iIrrep,jIrrep,2) = ipTxy(iIrrep,jIrrep,1)+n_Txy
        ipTxy(jIrrep,iIrrep,2) = ipTxy(iIrrep,jIrrep,2)
      end if
      iOff2 = iOff2+nnAorb
    end do
    iOff = iOff+iOff2*nnP(kIrrep)
    mAO = mAO+nBas(kIrrep)*nASh(kIrrep)
  end do

  allocate(AOrb(nADens))
  do iADens=1,nADens
    call Allocate_DT(AOrb(iADens),nAsh,nBas,nIrrep)
  end do

  ! Reordering of the active MOs :  C(a,v) ---> C(v,a)

  iCount = 0
  do iIrrep=0,nIrrep-1

    jCount = iCount+nBas(iIrrep)*nIOrb(iIrrep)

    iCount = iCount+nBas(iIrrep)**2
    if (nBas(iIrrep)*nASh(iIrrep) == 0) cycle

    do iADens=1,nADens

      do i=1,nASh(iIrrep)
        kOff1 = 1+jCount+nBas(iIrrep)*(i-1)
        AOrb(iADens)%SB(iIrrep+1)%A2(i,1:nBas(iIrrep)) = CMO(kOff1:kOff1-1+nBas(iIrrep),iADens)
      end do

    end do ! iADens

  end do ! iIrrep

  if (nKdens == 2) then ! for Coulomb term
    call daxpy_(nDens,One,DMLT(2)%A0,1,DMLT(1)%A0,1)
  end if

  ! Add the density of the environment in a OFE calculation (Coulomb)

  call OFembed_dmat(DMlt(1)%A0,nDens)

  !nScreen = 10 ! Some default values for the screening parameters
  !dmpK = One
  Estimate = .false.
  Update = .true.
  call Cho_Get_Grad(irc,nKdens,DMLT,DLT2,ChM,Txy,n_Txy*nAdens,ipTxy,DoExchange,lSA,nChOrb,AOrb,nAsh,DoCAS,Estimate,Update, &
                    V_k(jp_V_k,1),nV_k,U_k(jp_U_k),Z_p_k(jp_Z_p_k,1),nZ_p_k,nnP,npos)

  if (irc /= 0) then
    call WarningMessage(2,'Compute_AuxVec: failed in Cho_Get_Grad !!')
    call Abend()
  end if

  if (DoCAS .or. DoExchange) then
    do i=1,nKdens
      call deallocate_DT(ChM(i))
    end do
  end if
  if (iMp2prpt == 2) then
    call deallocate_DT(DLT2)
  end if

end if ! no vectors on this node

! For parallel run: reordering of the V_k(tilde) vector from
! the "node storage" to the Q-vector storage
if (nProc > 1) then
  do i=1,size(V_K,2)
    call Reord_Vk(ipVk,nProc,myProc,nV_l,nV_t,[1],1,V_k(:,i))
  end do
end if
!***********************************************************************
!                                                                      *
!     Second step: contract with the Q-vectors to produce V_k          *
!            ~  T                                                      *
!     V = V Q                                                          *
!                                                                      *
!***********************************************************************

if (Cholesky .and. (.not. Do_RI)) then ! to cope with the calls below
  nBas_Aux(0) = nBas_Aux(0)+1
end if

call mma_maxDBLE(MemMax)

if (Force_out_of_Core) MemMax = 4*(nQvMax)/10
nQv = min(MemMax,nQvMax)
call mma_allocate(Qv,nQv,Label='Qv')

! Coulomb

do i=0,nJdens-1
  call Mult_Vk_Qv_s(V_k(ipVk(1),1+i),nV_t(0),Qv,nQv,Scr,nQMax,nBas_Aux,nV_t(0),nIrrep,'T')
  call dcopy_(nV_k,Scr,1,V_k(ipVk(1),1+i),1)
end do

! MP2

if (iMp2prpt == 2) then
  call Mult_Vk_Qv_s(U_k(ipUk(1)),nU_t(0),Qv,nQv,Scr,nQMax,nBas_Aux,nU_t(0),nIrrep,'T')
  call dcopy_(nV_k,Scr,1,U_k(ipUk(1)),1)
end if

! Active term

if (DoCAS) then ! reorder Zp_k

  call mma_allocate(Zv,nZ_p_k,Label='Zv')

  do iAvec=1,nAvec
    if (nProc > 1) call Reord_Vk(ipZpk(1),nProc,myProc,nV_l,nV_t,nnP,nIrrep,Z_p_k(:,iAVec))

    call Mult_Zp_Qv_s(Z_p_k(ipZpk(1),iAvec),nZ_p_k,Qv,nQv,Zv,nZ_p_k,nV_t,nnP,nBas_Aux,nIrrep,'T')

    call dcopy_(nZ_p_k,Zv,1,Z_p_k(ipZpk(1),iAvec),1)
  end do
  call mma_deallocate(Zv)
end if

call mma_deallocate(Qv)
call mma_deallocate(Scr)

! Exchange

if (DoExchange) then
  DoCholExch = .true.
  do iSO=1,nKvec
    call Mult_RijK_QKL(iSO,nBas_aux,nIrrep)
  end do
  if (iMp2prpt == 2) then
    call Mult_with_Q_MP2(nBas_aux,nBas,nIrrep)
  end if
end if
if (Cholesky .and. (.not. Do_RI)) nBas_Aux(0) = nBas_Aux(0)-1

return

end subroutine Compute_AuxVec
