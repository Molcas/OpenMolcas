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
#include "WrkSpc.fh"
      Real*8 rkappa(nDensC), sigma(ndensC),FA(ndens2),D(*),
     &       p12(*),p11(*),rm1(*),rm2(*),Focki(*)
*
      irec(i,j)=i+nna*(j-1)
*
      isym=1
*     sign1  Kappa**t=Sign1*Kappa
*     sign2  <0|[Qip,H]|0>=Aip+sign2*Api
      r1=DBLE(i1)
      r2=DBLE(i2)
      Fact=-1.0d0 ! bullshit
      reco=-1.0d0 !(k*a+reco*a*k)
      jspin=1 ! triplet
      Call GetMem('ATemp','ALLO','REAL',ipK,ndens2)
      Call GetMem('ATemp','ALLO','REAL',ipFA,ndens2)
      Call GetMem('FTemp','ALLO','REAL',ipFock,ndens2)
      Call GetMem('QTemp','ALLO','REAL',ipq,ndens2)
      Call GetMem('QTemp','ALLO','REAL',ipq1,ndens2)
      call dcopy_(nmba,0.0d0,0,rm1,1)
      call dcopy_(nmba,0.0d0,0,rm2,1)
      call dcopy_(ndens2,0.0d0,0,Focki,1)
      call dcopy_(ndens2,0.0d0,0,Work(ipQ),1)
      call dcopy_(ndens2,0.0d0,0,Work(ipQ1),1)
      Call Unc(rkappa,Work(ipK),isym,r1)

      Call R2ElInt_SP(Work(ipK),rm1,rm2,
     %             FockI,Work(ipFA),
     &             nDens2,iSym,ReCo,Fact,jspin,D,FA)
*
*
*
      call dcopy_(ndens2,0.0d0,0,Work(ipFock),1)
*
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
*                          __
*     <o|E(S)  E- E(S) |o>(pb|cd)
*            ab cd    ad
      Call CreQ_sp(Work(ipQ),rm1,P11,isym)
      Call DSCAL_(ndens,r3,Work(ipQ),1)

*                             __
*     <o|E(S)  E- E(S) |o>(pb|cd)
*            ab cd   ad
      Call CreQ_sp(Work(ipQ1),rm2,p12,isym)
      call daxpy_(ndens,r4,Work(ipQ1),1,Work(ipq),1)
*
      Do iS=1,nSym
*
*             A
*      F  =2 F
*       pi    pi
*
       Call DaXpY_(nIsh(is)*nBas(is),-r1*2.0d0,
     &            Work(ipFA-1+ipMat(is,is)),1,
     &            Work(ipFock+ipMat(is,is)-1),1)

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
     &               Work(ipfock-1+ipF1),1)
          Call DaXpY_(nIsh(is),-Dij,
     &               FockI(ipFI1),1,
     &               Work(ipfock-1+ipF2),nbas(is))
        End Do
       End Do
*
*
*      F  = F  + Q
*       pa   pa   pa
*
       Call DaXpY_(nAsh(is)*nBas(is),-r1,
     &            Work(nbas(is)*nish(is)+ipQ+ipMat(is,is)-1),1,
     &            Work(ipFock-1+ipMat(is,is)+nBas(is)*nIsh(is)),1)
       Do iA=nish(is),nish(is)+nAsh(is)-1
        Call DaXpY_(nBas(is),-1.0d0,
     &            Work(ipQ+nbas(is)*ia+ipMat(is,is)-1),1,
     &            Work(ipFock+ipMat(is,is)-1+iA),nbas(is))

       End Do
       End Do
*
      Call Compress(Work(ipFock),Sigma,isym)
*
      Call GetMem('QTemp','Free','REAL',ipq,ndens2)
      Call GetMem('QTemp','Free','REAL',ipq1,ndens2)
      Call GetMem('ATemp','FREE','REAL',ipK,ndens2)
      Call GetMem('ATemp','FREE','REAL',ipFA,ndens2)
      Call GetMem('FTemp','FREE','REAL',ipFock,ndens2)
*
      return
      end
