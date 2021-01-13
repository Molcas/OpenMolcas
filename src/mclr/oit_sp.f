************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine oit_sp(rkappa,sigma,i1,i2,r3,p11,r4,p12,D,FA,
     &                  rm1,rm2,focki)
*
*                              ~
*     Constructs  F  = <0|[Q  ,H]|0>
*                  pq       pq
*
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
      Real*8 rkappa(nDensC), sigma(ndensC),FA(ndens2),D(*),
     &       p12(*),p11(*),rm1(*),rm2(*),Focki(*)
      Real*8, Allocatable:: K(:), FAtemp(:), Fock(:), Q(:), Q1(:)
*
      irec(i,j)=i+nna*(j-1)
*
      isym=1
*     sign1  Kappa**t=Sign1*Kappa
*     sign2  <0|[Qip,H]|0>=Aip+sign2*Api
      r1=DBLE(i1)
      Fact=-1.0d0 ! bullshit
      reco=-1.0d0 !(k*a+reco*a*k)
      jspin=1 ! triplet
      Call mma_allocate(K,ndens2,Label='K')
      Call mma_allocate(FAtemp,ndens2,Label='FAtemp')
      Call mma_allocate(Fock,ndens2,Label='Fock')
      Call mma_allocate(Q,ndens2,Label='Q')
      Call mma_allocate(Q1,ndens2,Label='Q1')
      call dcopy_(nmba,[0.0d0],0,rm1,1)
      call dcopy_(nmba,[0.0d0],0,rm2,1)
      call dcopy_(ndens2,[0.0d0],0,Focki,1)
      Q(:)=0.0D0
      Q1(:)=0.0D0
      Call Unc(rkappa,K,isym,r1)

      Call R2ElInt_SP(K,rm1,rm2,
     %             FockI,FAtemp,
     &             nDens2,iSym,ReCo,Fact,jspin,D,FA)
*
*
*
      call dcopy_(ndens2,[0.0d0],0,Fock,1)
*
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
*                          __
*     <o|E(S)  E- E(S) |o>(pb|cd)
*            ab cd    ad
      Call CreQ_sp(Q,rm1,P11,isym)
      Call DSCAL_(ndens,r3,Q,1)

*                             __
*     <o|E(S)  E- E(S) |o>(pb|cd)
*            ab cd   ad
      Call CreQ_sp(Q1,rm2,p12,isym)
      call daxpy_(ndens,r4,Q1,1,Q,1)
*
      Do iS=1,nSym
*
*             A
*      F  =2 F
*       pi    pi
*
       Call DaXpY_(nIsh(is)*nBas(is),-r1*2.0d0,
     &            FAtemp(ipMat(is,is)),1,
     &            Fock(ipMat(is,is)),1)

       Do iAsh=1,nAsh(iS)
        Do jAsh=1,nAsh(is)
         Dij=D(irec(iash+nA(is),jAsh+nA(is)))
         ipF1= ipMat(is,is)+(Nish(is)+iAsh-1)*nBas(is)
         ipF2= ipMat(is,is)+Nish(is)+iAsh-1
         ipFI1=ipMat(is,is)+(Nish(is)+jAsh-1)*nBas(is)
*
*                I
*        F  = F + F  D
*         pa   pa  pb ab
*
          Call DaXpY_(nBas(is),-r1*Dij,
     &               FockI(ipFI1),1,
     &               Fock(ipF1),1)
          Call DaXpY_(nIsh(is),-Dij,
     &               FockI(ipFI1),1,
     &               Fock(ipF2),nbas(is))
        End Do
       End Do
*
*
*      F  = F  + Q
*       pa   pa   pa
*
       Call DaXpY_(nAsh(is)*nBas(is),-r1,
     &            Q(nbas(is)*nish(is)+ipMat(is,is)),1,
     &            Fock(ipMat(is,is)+nBas(is)*nIsh(is)),1)
       Do iA=nish(is),nish(is)+nAsh(is)-1
        Call DaXpY_(nBas(is),-1.0d0,
     &            Q(nbas(is)*ia+ipMat(is,is)),1,
     &            Fock(ipMat(is,is)+iA),nbas(is))

       End Do
       End Do
*
      Call Compress(Fock,Sigma,isym)
*
      Call mma_deallocate(Q1)
      Call mma_deallocate(Q)
      Call mma_deallocate(K)
      Call mma_deallocate(FAtemp)
      Call mma_deallocate(FOck)
*
      return
* Avoid unused argument warnings
      if (.false.) call Unused_integer(i2)
      end
