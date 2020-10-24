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
      SubRoutine RInt_ns(rkappa,rmo,Fock,Focki,idsym,reco,jspin,rie)
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
#include "real.fh"
#include "stdalloc.fh"
#include "glbbas_mclr.fh"
      Real*8 Fock(nDens2),rkappa(nDens2),
     &       Focki(ndens2),rMO(*)
      Real*8, Allocatable:: FA(:), MT1(:), MT2(:), QA(:), QB(:)
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      Call mma_allocate(FA,ndens2,Label='FA')
c Fact controls the sign of H(k)
      Fact=One
      Call mma_allocate(MT1,nmba,Label='MT1')
      Call mma_allocate(MT2,nmba,Label='MT2')
      MT1(:) = Zero
      MT2(:) = Zero

      Call R2ElInt_ns(rkappa,MT1,MT2,focki,FA,
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

      Fock(:) = Zero
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
      If (iMethod.eq.2) Then
       Call mma_allocate(qA,ndens2,Label='QA')
       Call mma_allocate(qB,ndens2,Label='QB')
       Call CreQ_td(QB,MT1,Work(ipG2sq),idsym)
       Call CreQ_td(QA,MT2,Work(ipG2sq),idsym)
      End If
*
*      Call RECPRT('QB',' ',QB,nDens2,1)
*      Call RECPRT('QA',' ',QA,nDens2,1)
*
      Do iS=1,nSym
       jS=iEOr(iS-1,idsym-1)+1
*
*            I    A
*      F  =2( F  + F  )
*       pi     pi   pi
*
       Call DGEADD2(Two,
     &              Focki(ipMat(is,js)),nBas(is),'N',
     &              Fock(ipMat(is,js)),nBas(is),'N',
     %              Fock(ipMat(is,js)),nBas(is),
     &              nbas(is),nish(js))
       Call DGEADD2(-Two,
     &             Focki(ipMat(is,js)),nBas(is),'N',
     &             Fock(ipMat(is,js)),nBas(is),'N',
     %             Fock(ipMat(is,js)),nBas(is),
     &             nish(is),nBas(js))
       If (iMethod.eq.2) Then
c
c 61-121 is all to do with multiconf
c
       Call DGEADD2(Two,
     &              FA(ipMat(is,js)),nBas(is),'N',
     &              Fock(ipMat(is,js)),nBas(is),'N',
     %              Fock(ipMat(is,js)),nBas(is),
     &              nbas(is),nIsh(js))
       Call DGEADD2(-Two,
     &             FA(ipMat(is,js)),nBas(is),'N',
     &             Fock(ipMat(is,js)),nBas(is),'N',
     %             Fock(ipMat(is,js)),nBas(is),
     &             nish(is),nBas(js))

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

       Call DGEADD(Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),'N',
     &              QB(ipMatba(is,js)),nBas(is),'N',
     %              Fock(ipMat(is,js)+nbas(is)*nish(js)),nBas(is),
     &              nBas(is),nAsh(js))
       Call DGESUB(Fock(ipMat(is,js)+nish(is)),nBas(is),'N',
     &              QA(ipMatba(js,is)),nBas(js),'T',
     %              Fock(ipMat(is,js)+nish(is)),nBas(is),
     &              nash(is),nBas(js))
      End If
* Transpose ipsc2
*     Call mma_allocate(T,nbas(is)*nbas(jS),Label='T')
*     call dcopy_(nbas(is)*nbas(jS),[0.0d0],0,T,1)
*     Call DGETMO(Fock(ipmat(is,js)),nbas(is),
*    &                nbas(is),nbas(js),T,
*    &                nbas(js))
*     call dcopy_(nBas(jS)*nBas(iS),T,1,Fock(ipmat(js,is)),1)
*     Call mma_deallocate(T)
*
      End Do
      If (imethod.eq.2) Then
         Call mma_deallocate(QA)
         Call mma_deallocate(QB)
      End If
*
      Call DSCAL_(ndens2,-Two,fock,1)
      Call AddGrad(rKappa,Fock,idsym,Two*(-fact))
      Call PickMO_td(MT1,rmo,idsym)
*
      Call mma_deallocate(MT2)
      Call mma_deallocate(MT1)
      Call mma_deallocate(FA)
*
      return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(rie)
      end
