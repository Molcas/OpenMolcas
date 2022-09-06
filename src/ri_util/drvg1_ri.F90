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
! Copyright (C) 2007, Roland Lindh                                     *
!***********************************************************************

subroutine Drvg1_RI(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
!     Object: superdriver for gradients for the RI/DF approximation    *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh, Dep. Chem. Phys., Lund University, Sweden  *
!             January '07                                              *
!                                                                      *
!***********************************************************************

use Basis_Info, only: nBas, nBas_Aux
use pso_stuff
use RICD_Info, only: Do_RI, Cholesky
use Symmetry_Info, only: nIrrep
use Para_Info, only: myRank, nProcs
use Data_Structures, only: Deallocate_DT
use ExTerm, only: iMP2prpt, LuAVector, LuBVector

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
#include "cholesky.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "exterm.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
real*8 Grad(nGrad), Temp(nGrad)
character*6 Fname
character*7 Fname2
character*8 Method
logical Found
integer nAct(0:7)
real*8, allocatable :: V_k_new(:,:), U_k_new(:)
integer, allocatable :: iZk(:), iVk(:), iUk(:)
real*8, allocatable :: DMTmp(:), Tmp(:)
integer, allocatable :: SO_ab(:), ij2(:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Compute_AuxVec(ipVk,ipZpk,myProc,nProc,ipUk)
    integer nProc, myProc
    integer ipVk(nProc), ipZpk(nProc)
    integer, optional :: ipUk(nProc)
  end subroutine Compute_AuxVec
  subroutine Effective_CD_Pairs(ij2,nij_Eff)
    integer, allocatable :: ij2(:,:)
    integer nij_Eff
  end subroutine Effective_CD_Pairs
  subroutine Drvg1_2Center_RI(Grad,Temp,nGrad,ij2,nij_Eff)
    integer nGrad, nij_Eff
    real*8 Grad(nGrad), Temp(nGrad)
    integer, allocatable :: ij2(:,:)
  end subroutine Drvg1_2Center_RI
  subroutine Drvg1_3Center_RI(Grad,Temp,nGrad,ij3,nij_Eff)
    integer nGrad, nij_Eff
    real*8 Grad(nGrad), Temp(nGrad)
    integer, allocatable :: ij3(:,:)
  end subroutine Drvg1_3Center_RI
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
DoCholExch = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 33
iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
call FZero(Temp,nGrad)
call mma_allocate(Tmp,nGrad,Label='Tmp')
!                                                                      *
!***********************************************************************
!                                                                      *
BufFrac = 0.1d0
call Cho_X_Init(irc,BufFrac)
if (irc /= 0) then
  call WarningMessage(2,' Drvg1_RI: Cho_X_Init failed')
  call Abend()
end if
!
!*******************************************
!
! Decide if it's MP2

iMp2Prpt = 0
call Get_cArray('Relax Method',Method,8)
if (Method == 'MBPT2   ') then
  call Get_iScalar('mp2prpt',iMp2Prpt)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of the Cholesky approach compute the A and Q matrices.

if (Cholesky .and. (.not. Do_RI)) then

  if (nIrrep /= 1) then
    call WarningMessage(2,'Error in Drvg1_RI')
    write(6,*) ' CD gradients with symmetry is not implemented yet!'
    call Abend()
  end if

  call Cho_X_CalculateGMat(irc)
  if (iRC /= 0) then
    call WarningMessage(2,'Error in Drvg1_RI')
    write(6,*) 'Failure during G matrix construction'
    call Abend()
  end if

  ! Now compute the Q matrix.
  !
  !     Note that, as the A matrix is
  !     computed in the full-pivoted (rows and columns) storage,
  !     also the resulting Q matrix is full-pivoted.
  !     This is necessary for the ReMap_V_k to work (see below).
  !     In the RI case, only the column pivoting of Q is
  !     preserved. One day we may want to unify the two cases.
  !
  !     (In Cholesky the Q matrix is stored as squared. In RI,
  !      it is, in general, rectangular as lin. dep. may occur
  !      among its columns).

  call ICopy(nIrrep,NumCho,1,nBas_Aux,1)
  call GAIGOP(nBas_Aux,nIrrep,'+')
  call Gen_QVec(nIrrep,nBas_Aux)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Prepare handling of two-particle density.

call PrepP()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize the number of sets of densities and auxiliary vectors
nAdens = 1
nAVec = 1
nKdens = 1
nJdens = 1

call Qpg_iScalar('SCF mode',Found)
if (Found) then
  call Get_iScalar('SCF mode',iUHF) ! either 0 or 1
else
  iUHF = 0
end if
nKdens = nKdens+iUHF
nKvec = nKdens

if (lPSO .and. lSA) then
  nJdens = 5
  nKdens = 4
  nKVec = 2
  nAdens = 2
  nAvec = 4
end if

!MGD Could be more efficient memory-wise when symmetry
!    Decompose the 2-particle active density matrix
mAO = 0
if (lPSO) then
  call Get_iArray('nAsh',nAct,nIrrep)
  n_Txy = 0
  do ijsym=0,nIrrep-1
    ntmp = 0
    do jSym=0,nIrrep-1
      isym = ieor(jSym,ijsym)
      if (iSym > jSym) then
        ntmp = ntmp+nAct(iSym)*nAct(jSym)
      else if (iSym == jSym) then
        ntmp = ntmp+nAct(iSym)*(nAct(iSym)+1)/2
      end if
    end do
    n_Txy = n_Txy+ntmp**2
    mAO = mAO+nAct(ijsym)*nBas(ijsym)
  end do
  m_Txy = nAdens
  call mma_allocate(Txy,n_Txy,nAdens,Label='Txy')
  call mma_allocate(DMdiag,nG1,nAdens,Label='DMdiag')
  call mma_allocate(DMtmp,nG1*(nG1+1)/2,Label='DMtmp')
  call iZero(nnP,nIrrep)
  call Compute_txy(G1(1,1),nG1,Txy,n_Txy,nAdens,nIrrep,DMdiag,DMtmp,nAct)
  call mma_deallocate(DMtmp)
else
  call mma_allocate(Txy,1,1,Label='Txy')
  call mma_allocate(DMdiag,1,1,Label='DMdiag')
end if
n_ij2K = 0
nZ_p_k = 0
nZ_p_l = 0
nZ_p_k_New = 0
do i=0,nIrrep-1
  iOff_ij2K(i+1) = n_ij2K
  n_ij2K = n_ij2K+nBas(i)*(nBas(i)+1)/2
  nZ_p_k = nZ_p_k+nnP(i)*nBas_Aux(i)  ! Global size
  nZ_p_l = nZ_p_l+nnP(i)*NumCho(i+1)  ! Local size
  nZ_p_k_New = nZ_p_k_New+nnP(i)*nBas(i)*(nBas(i)+1)/2
end do
if (Do_RI) nZ_p_k = nZ_p_k-nnP(0)

! Allocate the "global" Z_p_k array

if (lPSO) then
  call mma_allocate(Z_p_k,nZ_p_k,nAVec,Label='Z_p_k')
else
  nZ_p_k = 1
  call mma_allocate(Z_p_k,1,nAVec,Label='Z_p_k')
end if
Z_p_k(:,:) = Zero

! Preprocess the RI and Q vectors as follows

! Allocate memory for V_k

nV_k = nBas_Aux(0)
if (Do_RI) nV_k = nV_k-1

nAux_Tot = 0
do iIrrep=0,nIrrep-1
  nAux_Tot = nAux_Tot+nBas_Aux(iIrrep)
end do

call mma_allocate(V_K,nV_k,nJdens,Label='V_k')
V_k(:,:) = Zero
if (iMp2prpt == 2) then
  call mma_allocate(U_K,nV_k,Label='U_k')
  U_k(:) = Zero
else
  call mma_allocate(U_K,1,Label='U_k')
end if
!                ~
! 1) Compute the V_k vector
!             ~
! 2) Contract V_k and Q (transpose) vectors producing the V_k

! Note: the above two points apply to Z_p_k as well (active space)

call mma_allocate(iVk,[0,nProcs-1],Label='iVk')
call mma_allocate(iZk,[0,nProcs-1],Label='iZk')
iVk(:) = 0
iZk(:) = 0
!iVk(myRank) = NumCho(1)*nJdens
iVk(myRank) = NumCho(1)
!iZk(myRank) = nZ_p_l*nAvec  ! store the local size of Zk
iZk(myRank) = nZ_p_l         ! store the local size of Zk
call GAIGOP(iVk,nProcs,'+')
call GAIGOP(iZk,nProcs,'+')  ! distribute to all nodes

! Compute the starting position in the global sense for each node.

iStart = 1
jStart = 1
do j=0,nProcs-1  ! Loop over all nodes
  itmp = iVk(j)
  iVk(j) = iStart
  iStart = iStart+itmp

  jtmp = iZk(j)
  iZk(j) = jStart
  jStart = jStart+jtmp
end do

if (iMp2prpt == 2) then
  call mma_allocate(iUk,[0,nProcs-1],Label='iUk')
  iUk(:) = 0
  iUk(myRank) = NumCho(1)
  call GAIGOP(iUk,nProcs,'+')
  kStart = 1
  do j=0,nProcs-1
    kTmp = iUk(j)
    iUk(j) = kStart
    kStart = kStart+kTmp
  end do

  call Compute_AuxVec(iVk,iZk,myRank+1,nProcs,ipUk=iUk)

else

  call Compute_AuxVec(iVk,iZk,myRank+1,nProcs)

end if
!                                                                      *
!***********************************************************************
!                                                                      *

if (Cholesky .and. (.not. Do_RI)) then

  ! Map from Cholesky auxiliary basis to the full
  ! 1-center valence product basis.

  call mma_allocate(ij2K,n_ij2K,Label='ij2K')
  ij2K(:) = 0
  nV_k_New = nBas(0)*(nBas(0)+1)/2
  call mma_allocate(V_k_new,nV_k_New,nJdens,Label='V_k_new')
  V_k_new(:,:) = Zero

  if (iMp2prpt == 2) then
    call mma_allocate(U_k_new,nV_k_New,Label='U_k_new')
    U_k_new(:) = Zero
  end if

  call mma_allocate(SO_ab,2*nAux_Tot,Label='SO_ab')
  SO_ab(:) = 0
  iOff = 1
  do iSym=1,nSym
    call CHO_X_GET_PARDIAG(iSym,SO_ab(iOff))

    if ((iSym == 1) .and. (iMp2prpt == 2)) then
      call ReMap_U_k(U_k,nV_k,U_k_New,nV_k_New,SO_ab)
    end if
    m_ij2K = nBas(iSym-1)*(nBas(iSym-1)+1)/2
    do i=0,nJDens-1
      call ReMap_V_k(iSym,V_k(1,1+i),nV_k,V_k_new(1,1+i),nV_k_New,SO_ab(iOff),ij2K(iOff_ij2K(iSym)+1),m_ij2K)
    end do
    iOff = iOff+2*nBas_Aux(iSym-1)
  end do

  nV_k = nV_k_new

  call mma_deallocate(SO_ab)
  call mma_deallocate(V_k)
  call mma_allocate(V_k,nV_k,nJdens,Label='V_k')
  V_k(:,:) = V_k_new(:,:)
  call mma_deallocate(V_k_new)

  if (iMp2prpt == 2) then
    call mma_deallocate(U_k)
    call mma_allocate(U_k,nV_k,Label='U_k')
    U_k(:) = U_k_new(:)
    call mma_deallocate(U_k_new)
  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Get the effective list of shell-pairs in case of CD

  call Effective_CD_Pairs(ij2,nij_Eff)
else
  nij_Eff = 0
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Open C-vector-files if nSym is equal to 1

if (DoCholExch) then
  do i=1,nKvec
    do jSym=1,nSym
      iSeed = 7+jSym+(i-1)*nSym
      LuCVector(jSym,i) = IsFreeUnit(iSeed)
      if (i == 1) then
        write(Fname,'(A4,I1,I1)') 'CVEA',jSym
      else if (i == 2) then
        write(Fname,'(A4,I1,I1)') 'CVEB',jSym
      end if
      call DANAME_MF_WA(LuCVector(jSym,i),Fname)
    end do
  end do
  ! Initialize timings
  do i=1,2
    tavec(i) = 0.0d0
    tbvec(i) = 0.0d0
  end do
end if
if (imp2prpt == 2) then
  do i=1,2
    iSeed = 8+nSym
    LuAVector(i) = IsFreeUnit(iSeed)
    write(Fname2,'(A5,I1)') 'AMP2V',i
    call DaName_MF_WA(LuAVector(i),Fname2)
    iSeed = 9+nSym
    LuBVector(i) = IsFreeUnit(iSeed)
    write(Fname,'(A5,I1)') 'BMP2V',i+2
    call DaName_MF_WA(LuBVector(i),Fname)
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute contributions due to the "2-center" two-electron integrals

Case_2C = .true.
call Drvg1_2center_RI(Temp,Tmp,nGrad,ij2,nij_Eff)
call GADGOP(Tmp,nGrad,'+')
if (iPrint >= 15) call PrGrad(' RI-Two-electron contribution - 2-center term',Tmp,nGrad,ChDisp)
call DaXpY_(nGrad,One,Temp,1,Grad,1) ! Move any 1-el contr.
call dcopy_(nGrad,Tmp,1,Temp,1)
call DScal_(nGrad,-One,Temp,1)
Case_2C = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute contributions due to the "3-center" two-electron integrals

Case_3C = .true.
call Drvg1_3center_RI(Temp,Tmp,nGrad,ij2,nij_Eff)
call GADGOP(Tmp,nGrad,'+')
if (iPrint >= 15) call PrGrad(' RI-Two-electron contribution - 3-center term',Tmp,nGrad,ChDisp)
call DaXpY_(nGrad,Two,Tmp,1,Temp,1)
Case_3C = .false.
if (allocated(Txy)) call mma_deallocate(Txy)
if (allocated(DMdiag)) call mma_deallocate(DMdiag)
if (allocated(AOrb)) then
  do iADens=1,nADens
    call Deallocate_DT(AOrb(iADens))
  end do
  deallocate(AOrb)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoCholExch) then
  do i=1,nKvec
    do jSym=1,nSym
      call DaClos(luCVector(jSym,i))
    end do
  end do
end if
if (iMp2prpt == 2) then
  do i=1,2
    call DaClos(LuAVector(i))
    call DaClos(LuBVector(i))
  end do
end if

if (Cholesky .and. (.not. Do_RI)) then
  call mma_deallocate(ij2)
  call mma_deallocate(ij2K)
end if
call CloseP()
call mma_deallocate(iZk)
call mma_deallocate(iVk)

if (allocated(iUk)) call mma_deallocate(iUk)
if (allocated(Z_p_k)) call mma_deallocate(Z_p_k)
if (allocated(V_k)) call mma_deallocate(V_k)
if (allocated(U_k)) call mma_deallocate(U_k)

call Cho_X_Final(irc)
if (irc /= 0) then
  call WarningMessage(2,' Drvg1_RI: Cho_X_Final failed')
  call Abend()
end if
call mma_deallocate(Tmp)
if (iPrint >= 15) call PrGrad(' RI-Two-electron contribution - Temp',Temp,nGrad,ChDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
call SavTim(6,TCpu2-TCpu1,TWall2-TWall1)

#ifdef _CD_TIMING_
Drvg1_CPU = TCpu2-TCpu1
Drvg1_Wall = TWall2-TWall1
write(6,*) '-------------------------'
write(6,*) 'Time spent in Cho_get_grad:'
write(6,*) 'Wall/CPU',ChoGet_Wall,ChoGet_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Mult_Rijk_Qkl:'
write(6,*) 'Wall/CPU',rMult_Wall,rMult_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Prepp:'
write(6,*) 'Wall/CPU',Prepp_Wall,Prepp_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Pget_ri2:'
write(6,*) 'Wall/CPU',Pget2_Wall,Pget2_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Pget_ri3:'
write(6,*) 'Wall/CPU',Pget3_Wall,Pget3_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Drvg1_ri:'
write(6,*) 'Wall/CPU',Drvg1_Wall,Drvg1_CPU
write(6,*) '-------------------------'
Total_Dens_Wall = ChoGet_Wall+rMult_Wall+Prepp_Wall+Pget2_Wall+Pget3_Wall
Total_Dens_CPU = ChoGet_CPU+rMult_CPU+Prepp_CPU+Pget2_CPU+Pget3_CPU
Total_Der_Wall = Drvg1_Wall-Total_Dens_Wall
Total_Der_CPU = Drvg1_CPU-Total_Dens_CPU
Total_Der_Wall2 = TwoEl2_Wall+TwoEl3_Wall
Total_Der_CPU2 = TwoEl2_CPU+TwoEl3_CPU

write(6,*) 'Total Time for Density:'
write(6,*) 'Wall/CPU',Total_Dens_Wall,Total_Dens_CPU
write(6,*) '-------------------------'
write(6,*) 'Total TIme for 2-center Derivatives:'
write(6,*) 'Wall/CPU',Twoel2_Wall,Twoel2_CPU
write(6,*) '-------------------------'
write(6,*) 'Total TIme for 3-center Derivatives:'
write(6,*) 'Wall/CPU',Twoel3_Wall,Twoel3_CPU
write(6,*) '-------------------------'
write(6,*) 'Total Time for Derivatives:'
write(6,*) 'Wall/CPU',Total_Der_Wall2,Total_Der_CPU2
write(6,*) '-------------------------'
write(6,*) 'Derivative check:'
write(6,*) 'Wall/CPU',Total_Der_Wall,Total_Der_CPU
write(6,*) '-------------------------'
#endif

return

end subroutine Drvg1_RI
