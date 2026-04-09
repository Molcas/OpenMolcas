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

use MCLR_Data, only: ipMat, nA, nDens, nDensC, nMBA, nNA
use input_mclr, only: nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rkappa(nDensC), R3, p11(*), R4, p12(*), D(*), FA(nDens)
real(kind=wp), intent(out) :: sigma(nDensC)
integer(kind=iwp), intent(in) :: i1
real(kind=wp), intent(_OUT_) :: rm1(*), rm2(*), Focki(*)
integer(kind=iwp) :: IA, iAsh, ipF1, ipF2, ipFI1, iS, iSym, jAsh, jSpin
real(kind=wp) :: Dij, Fact, R1, Reco
real(kind=wp), allocatable :: FAtemp(:), Fock(:), K(:), Q(:), Q1(:)

isym = 1
! sign1  Kappa**t=Sign1*Kappa
! sign2  <0|[Qip,H]|0>=Aip+sign2*Api
r1 = real(i1,kind=wp)
Fact = -One ! bullshit
reco = -One !(k*a+reco*a*k)
jspin = 1 ! triplet
call mma_allocate(K,nDens,Label='K')
call mma_allocate(FAtemp,nDens,Label='FAtemp')
call mma_allocate(Fock,nDens,Label='Fock')
call mma_allocate(Q,nDens,Label='Q')
call mma_allocate(Q1,nDens,Label='Q1')
rm1(1:nmba) = Zero
rm2(1:nmba) = Zero
Focki(1:nDens) = Zero
Q(:) = Zero
Q1(:) = Zero
call Unc(rkappa,K,isym,r1)

call R2ElInt_SP(K,rm1,rm2,FockI,FAtemp,iSym,ReCo,Fact,jspin,D,FA)

Fock(:) = Zero

! Q  = sum(jkl)=(pj|kl)d(ijkl)
!  pi

! <o|E(S)  E- E(S) |o>(pb|cd)
!        ab cd    ad
call CreQ_sp(Q,rm1,P11,isym)
Q(:) = r3*Q(:)

! <o|E(S)  E- E(S) |o>(pb|cd)
!        ab cd   ad
call CreQ_sp(Q1,rm2,p12,isym)
Q(:) = Q(:)+r4*Q1(:)

do iS=1,nSym

  !        A
  ! F  =2 F
  !  pi    pi

  Fock(ipMat(is,is):ipMat(is,is)+nIsh(is)*nBas(is)-1) = Fock(ipMat(is,is):ipMat(is,is)+nIsh(is)*nBas(is)-1)- &
                                                        r1*Two*FAtemp(ipMat(is,is):ipMat(is,is)+nIsh(is)*nBas(is)-1)

  do iAsh=1,nAsh(iS)
    do jAsh=1,nAsh(is)
      Dij = D(iash+nA(is)+(jAsh+nA(is)-1)*nNA)
      ipF1 = ipMat(is,is)+(Nish(is)+iAsh-1)*nBas(is)
      ipF2 = ipMat(is,is)+Nish(is)+iAsh-1
      ipFI1 = ipMat(is,is)+(Nish(is)+jAsh-1)*nBas(is)

      !         I
      ! F  = F + F  D
      !  pa   pa  pb ab

      Fock(ipF1:ipF1+nBas(is)-1) = Fock(ipF1:ipF1+nBas(is)-1)-r1*Dij*FockI(ipFI1:ipFI1+nBas(is)-1)
      Fock(ipF2:ipF2+nIsh(is)*nBas(is)-1:nBas(is)) = Fock(ipF2:ipF2+nIsh(is)*nBas(is)-1:nBas(is))-Dij*FockI(ipFI1:ipFI1+nIsh(is)-1)
    end do
  end do

  ! F  = F  + Q
  !  pa   pa   pa

  ipF1 = ipMat(is,is)+nBas(is)*nIsh(is)
  ipF2 = ipMat(is,is)+nBas(is)*(nIsh(is)+nAsh(is))-1
  Fock(ipF1:ipF2) = Fock(ipF1:ipF2)-r1*Q(ipF1:ipF2)
  do iA=nIsh(is),nIsh(is)+nAsh(is)-1
    Fock(ipMat(is,is)+iA:ipMat(is,is)+iA+nBas(is)**2-1:nBas(is)) = &
      Fock(ipMat(is,is)+iA:ipMat(is,is)+iA+nBas(is)**2-1:nBas(is))-Q(ipMat(is,is)+nBas(is)*iA:ipMat(is,is)+nBas(is)*(iA-1)-1)
  end do
end do

call Compress(Fock,Sigma,isym)

call mma_deallocate(Q1)
call mma_deallocate(Q)
call mma_deallocate(K)
call mma_deallocate(FAtemp)
call mma_deallocate(FOck)

end subroutine oit_sp
