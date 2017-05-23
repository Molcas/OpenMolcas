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
      SubRoutine RInt_ns(rkappa,rmo,Fock,Focki,
     &               idsym,reco,jspin,rie)
*
*                              ~
*     Constructs  F  = <0|[E  ,H]|0>
* added in rinttd ( + <0|[[E  , Kappa],H]|0> )
*                  pq       pq                 pq
*
* Some modifications (subroutines that ends with _ns) to handle
* non anti-symmetric orbital rotations.
* Instead of using Q-Q^t as the "Q" contribution to the gradient
* Q^A-Q^{Bt} is used (see the paper). Major modifications in read2.
* Instead of constructing one set of MO integrals/particle we construct
* one set of integrals used for constructing Q^A and one for Q^B
*
c
c Fock is E*d/dx(lambda)
c rkappa is d/dx(lambda)
c
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"
      Real*8 Fock(nDens2),rkappa(nDens2),
     &       Focki(ndens2),rMO(*)
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      Call GETMEM('FA','ALLO','REAL',ipFA,ndens2)
c Fact controls the sign of H(k)
      Fact=1.0d0
      One=1.0d0
      Call GetMem('MOTemp1','ALLO','REAL',ipMT1,nmba)
      Call GetMem('MOTemp2','ALLO','REAL',ipMT2,nmba)
      call dcopy_(nmba,0.0d0,0,Work(ipMT1),1)
      call dcopy_(nmba,0.0d0,0,Work(ipMT2),1)

*      Call GetMem('MOTemp2','CHEC','REAL',ipMT2,nmba)
      Call R2ElInt_ns(rkappa,work(ipMT1),work(ipMT2),
     %             focki,Work(ipFA),
     &             nDens2,idSym,ReCo,Fact,jspin)
C
      If (.false.) Then  !If (idsym.eq.2) Then
      jpCMO=1
       Do iSym=1,nSym
            Write(6,'(A,i2.2)') 'Inactive fackmatrix = ',iSym
            Call RecPrt(' ' ,' ',focki(jpCMO),nBas(iSym),nBas(iSym))
            jpCMO=jpCMO+nBas(iSym)*nBas(iSym)
       End do
       jpCMO=ipFA
       Do  iSym=1,nSym
            Write(6,'(A,i2.2)') 'Active fockmatrix     = ',iSym
            Call RecPrt(' ' ,' ',Work(jpCMO),nBas(iSym),nBas(iSym))
            jpCMO=jpCMO+nBas(iSym)*nBas(iSym)
       end do
      End If
C

*      Call GetMem('MO1emp1','CHEC','REAL',ipMT2,nmba)
*
*
*
      call dcopy_(ndens2,0.0d0,0,Fock,1)
*
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
      If (iMethod.eq.2) Then
       Call GetMem('QTemp','ALLO','REAL',ipqA,ndens2)
       Call GetMem('QTemp','ALLO','REAL',ipqB,ndens2)
       Call CreQ_td(Work(ipQB),Work(ipMT1),Work(ipG2sq),idsym)
       Call CreQ_td(Work(ipQA),Work(ipMT2),Work(ipG2sq),idsym)
      End If
*
*      Call RECPRT('IpQB',' ',Work(ipQB),nDens2,1)
*      Call RECPRT('IpQA',' ',Work(ipQA),nDens2,1)
*      Call GetMem('MO12mp1','CHEC','REAL',ipMT2,nmba)
*
      Do iS=1,nSym
       jS=iEOr(iS-1,idsym-1)+1
*
*            I    A
*      F  =2( F  + F  )
*       pi     pi   pi
*
       Call DGEADD2(2.0d0,
     &              Focki(ipMat(is,js)),nBas(is),'N',
     &              Fock(ipMat(is,js)),nBas(is),'N',
     %              Fock(ipMat(is,js)),nBas(is),
     &              nbas(is),nish(js))
       Call DGEADD2(-2.0d0,
     &             Focki(ipMat(is,js)),nBas(is),'N',
     &             Fock(ipMat(is,js)),nBas(is),'N',
     %             Fock(ipMat(is,js)),nBas(is),
     &             nish(is),nBas(js))
*      Call GetMem('MO1e2p1','CHEC','REAL',ipMT2,nmba)
       If (iMethod.eq.2) Then
c
c 61-121 is all to do with multiconf
c
       Call DGEADD2(2.0d0,
     &              Work(ipFA-1+ipMat(is,js)),nBas(is),'N',
     &              Fock(ipMat(is,js)),nBas(is),'N',
     %              Fock(ipMat(is,js)),nBas(is),
     &              nbas(is),nIsh(js))
       Call DGEADD2(-2.0d0,
     &             Work(ipFA-1+ipMat(is,js)),nBas(is),'N',
     &             Fock(ipMat(is,js)),nBas(is),'N',
     %             Fock(ipMat(is,js)),nBas(is),
     &             nish(is),nBas(js))
*      Call GetMem('MO1e241','CHEC','REAL',ipMT2,nmba)

       Do iAsh=1,nAsh(jS)
        Do jAsh=1,nAsh(js)

*                I
*        F  = F - F  D
*         ap   ap  ap ba
*
          Dij=Work(ipg1t+itri(iash+nA(js),jAsh+nA(js))-1)
          ipF= ipMat(is,js)+(Nish(js)+iAsh-1)*nBas(is)
          ipFI=ipMat(is,js)+(Nish(js)+jAsh-1)*nBas(is)
*
*                I
*        F  = F + F  D
*         pa   pa  pb ab
*
          Call DaXpY_(nBas(is),Dij,
     &               focki(ipFI),1,
     &               Fock(ipF),1)
        End Do
       End Do
       Do iAsh=1,nAsh(iS)
        Do jAsh=1,nAsh(is)
         ipF=ipMat(is,js)+nIsh(is)+jAsh-1
         ipFI=ipMat(is,js)+nIsh(is)+iAsh-1
         Dij=Work(ipg1t+itri(iash+nA(is),jAsh+nA(is))-1)

*
*                I
*        F  = F - F  D
*         pa   pa  pb ab
*
          Call DaXpY_(nBas(js),-Dij,
     &               focki(ipFI),nbas(is),
     &               Fock(ipF),nbas(is))
        End Do
       End Do
*      Call GetMem('MO12241','CHEC','REAL',ipMT2,nmba)

       Call DGEADD(Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),'N',
     &              Work(ipQB+ipMatba(is,js)-1),nBas(is),'N',
     %              Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),
     &              nBas(is),nAsh(js))
       Call DGESUB(Fock(ipMat(is,js)+nish(is)),nBas(is),'N',
     &              Work(ipQA+ipMatba(js,is)-1),nBas(js),'T',
     %              Fock(ipMat(is,js)+nish(is)),nBas(is),
     &              nash(is),nBas(js))
      End If
* Transpose ipsc2
*     Call GetMem('Temp','ALLO','REAL',ipT,nbas(is)*nbas(jS))
*     call dcopy_(nbas(is)*nbas(jS),0.0d0,0,Work(ipT),1)
*     Call DGETMO(Fock(ipmat(is,js)),nbas(is),
*    &                nbas(is),nbas(js),Work(ipT),
*    &                nbas(js))
*     call dcopy_(nBas(jS)*nBas(iS),Work(ipT),
*    &       1,Fock(ipmat(js,is)),1)
*     Call GetMem('Temp','FREE','REAL',ipT,nbas(is)*nbas(jS))
*
      End Do
*      Call GetMem('M112241','CHEC','REAL',ipMT2,nmba)
      If (imethod.eq.2) Then
      Call GetMem('QTEMP','FREE','REAL',ipQA,ndens2)
      Call GetMem('QTEMP','FREE','REAL',ipQB,ndens2)
      End If

*
      Call DSCAL_(ndens2,-2.0d0,fock,1)
      Call AddGrad(rKappa,Fock,idsym,2.0d0*(-fact))
      Call PickMO_td(Work(ipMT1),rmo,idsym)
*
      Call GetMem('MOTemp2','FREE','REAL',ipMT2,nmba)
      Call GetMem('MOTemp1','FREE','REAL',ipMT1,nmba)
      Call GETMEM('FA','FREE','REAL',ipFA,ndens2)
*
      return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(rie)
      end
