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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine oit_sp(rkappa,sigma,i1,r3,p11,r4,p12,D,FA,rm1,rm2,focki)
!                          ~
! Constructs  F  = <0|[Q  ,H]|0>
!              pq       pq

use MCLR_Data, only: nDensC, nDens2, nNA, ipMat, nA, nDens, nMBA
use input_mclr, only: nSym, nAsh, nBas, nIsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp

implicit none
integer i1
real*8 rkappa(nDensC), sigma(ndensC), FA(ndens2), D(*), p12(*), p11(*), rm1(*), rm2(*), Focki(*)
real*8, allocatable :: K(:), FAtemp(:), Fock(:), Q(:), Q1(:)
integer iSym, jSpin, iS, iAsh, jAsh, ipF1, ipF2, ipFI1, IA
real*8 R1, Fact, Reco, R3, R4, Dij

isym = 1
! sign1  Kappa**t=Sign1*Kappa
! sign2  <0|[Qip,H]|0>=Aip+sign2*Api
r1 = real(i1,kind=wp)
Fact = -One ! bullshit
reco = -One !(k*a+reco*a*k)
jspin = 1 ! triplet
call mma_allocate(K,ndens2,Label='K')
call mma_allocate(FAtemp,ndens2,Label='FAtemp')
call mma_allocate(Fock,ndens2,Label='Fock')
call mma_allocate(Q,ndens2,Label='Q')
call mma_allocate(Q1,ndens2,Label='Q1')
call dcopy_(nmba,[Zero],0,rm1,1)
call dcopy_(nmba,[Zero],0,rm2,1)
call dcopy_(ndens2,[Zero],0,Focki,1)
Q(:) = Zero
Q1(:) = Zero
call Unc(rkappa,K,isym,r1)

call R2ElInt_SP(K,rm1,rm2,FockI,FAtemp,iSym,ReCo,Fact,jspin,D,FA)

call dcopy_(ndens2,[Zero],0,Fock,1)

! Q  = sum(jkl)=(pj|kl)d(ijkl)
!  pi

! <o|E(S)  E- E(S) |o>(pb|cd)
!        ab cd    ad
call CreQ_sp(Q,rm1,P11,isym)
call DSCAL_(ndens,r3,Q,1)

! <o|E(S)  E- E(S) |o>(pb|cd)
!        ab cd   ad
call CreQ_sp(Q1,rm2,p12,isym)
call daxpy_(ndens,r4,Q1,1,Q,1)

do iS=1,nSym

  !        A
  ! F  =2 F
  !  pi    pi

  call DaXpY_(nIsh(is)*nBas(is),-r1*Two,FAtemp(ipMat(is,is)),1,Fock(ipMat(is,is)),1)

  do iAsh=1,nAsh(iS)
    do jAsh=1,nAsh(is)
      Dij = D(iash+nA(is)+(jAsh+nA(is)-1)*nNA)
      ipF1 = ipMat(is,is)+(Nish(is)+iAsh-1)*nBas(is)
      ipF2 = ipMat(is,is)+Nish(is)+iAsh-1
      ipFI1 = ipMat(is,is)+(Nish(is)+jAsh-1)*nBas(is)

      !         I
      ! F  = F + F  D
      !  pa   pa  pb ab

      call DaXpY_(nBas(is),-r1*Dij,FockI(ipFI1),1,Fock(ipF1),1)
      call DaXpY_(nIsh(is),-Dij,FockI(ipFI1),1,Fock(ipF2),nbas(is))
    end do
  end do

  ! F  = F  + Q
  !  pa   pa   pa

  call DaXpY_(nAsh(is)*nBas(is),-r1,Q(nbas(is)*nish(is)+ipMat(is,is)),1,Fock(ipMat(is,is)+nBas(is)*nIsh(is)),1)
  do iA=nish(is),nish(is)+nAsh(is)-1
    call DaXpY_(nBas(is),-One,Q(nbas(is)*ia+ipMat(is,is)),1,Fock(ipMat(is,is)+iA),nbas(is))

  end do
end do

call Compress(Fock,Sigma,isym)

call mma_deallocate(Q1)
call mma_deallocate(Q)
call mma_deallocate(K)
call mma_deallocate(FAtemp)
call mma_deallocate(FOck)

end subroutine oit_sp
