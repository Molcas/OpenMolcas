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
      SubRoutine r2elint_ns(rKappa,rMO1,rmo2,FockI,FockA,nF,
     &                   iDSym,sign,Fact,jspin)
*
************************************************************************
*
*     Constructs the one index transformed Fock-matrixes
*     and (pj|kl).
*     rKappa : the transformation matrix
*     iDSym  : Symmetry of transformation
*     sign   : 1:real -1:complex
*     jspin  : 0:singlet 1:triplet
*
************************************************************************
*
      use Arrays, only: CMO, G1t, FAMO, FIMO
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "stdalloc.fh"
      Logical lFI,lFA,lMo
      Parameter ( One = 1.0d0 )
      Real*8 rKappa(nDens2),rMO1(nMba),rmo2(*),FockI(nDens2),
     &       FockA(nDens2)
      Real*8 rdum(1)
      Real*8, Allocatable:: T1(:), Tmp2(:), T3(:), T4(:), DIL(:),
     &                      DI(:), DIR(:), FI(:), Scr(:), FI1(:),
     &                      K1(:), DAL(:), DAR(:), DA(:), FA1(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      ndens22=ndens2
      iAM=0
      iBM=0
      Do i=1,nSym
       Do j=1,nSym
        nDens22=Max(nDens22,nBas(i)*nBas(j))
       End do
       iAM=Max(iAM,nAsh(i)+nIsh(i))
       iBM=Max(iBM,nBAs(i))
      End Do
      imem=nDens22

      Call mma_allocate(T1,imem,Label='T1')
      Call mma_allocate(Scr,imem,Label='Scr')
      Call mma_allocate(Tmp2,nDens22,Label='Tmp2')
      Call mma_allocate(T3,nDens22,Label='T3')
      Call mma_allocate(T4,nDens22,Label='T4')
      Call mma_allocate(DIL,nDens2,Label='DIL')
      Call mma_allocate(DI,nCMO,Label='DI')
      Call mma_allocate(DIR,nDens2,Label='DIR')
      Call mma_allocate(FI,ndens2,Label='FI')
      Call mma_allocate(FI1,ndens2,Label='FI1')
      Call mma_allocate(K1,ndens2,Label='K1')

      FockI(:)=0.0d0
      FockA(:)=0.0d0
      FI(:)   =0.0d0
      FI1(:)  =0.0d0
      K1(:)   =0.0d0
      DIR(:)  =0.0d0
      DIL(:)  =0.0d0
      DI(:)   =0.0d0
      lFI=.true.
      lFa=.false.
      lMo=.false.

      If (iMethod.eq.2) Then
         Call mma_allocate(DAL,nDens2,Label='DAL')
         Call mma_allocate(DAR,nDens2,Label='DAR')
         Call mma_allocate(DA,nCMO,Label='DA')
         Call mma_allocate(FA1,nDens2,Label='FA1')
         lFa=.true.
         lMo=.true.
      Else
         Call mma_allocate(DAL,1,Label='DAL')
         Call mma_allocate(DAR,1,Label='DAR')
         Call mma_allocate(DA,1,Label='DA')
         Call mma_allocate(FA1,1,Label='FA1')
      End If
      FA1(:)=0.0D0
      DAL(:)=0.0D0
      DAR(:)=0.0D0
      DA(:) =0.0D0
c THIS IS THE ENTIRE DENSITY FOR MULTICONF
      Do iS=1,nSym
         Do iB=1,nIsh(iS)
            ip=ipCM(iS)+(ib-1)*nBas(is)+ib-1
            DI(ip)=2.0d0
         End Do
      End Do
      If (iMethod.eq.2) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
          iA=nA(is)+ib
          jA=nA(is)+jb
          ip2=itri(iA,jA)
          DA(ip)=G1t(ip2)
         End Do
        End Do
       End Do
      End If

*
*     Construct {kappa,(ij|kl)} and the fock matrix contribution
*     from one index tranformation of contracted indexes
*
*      If (iDsym.eq.2) Then
*           Call RecPrt('FI',' ',FI,ndens2,1)
*           Call RecPrt('DI',' ',DI,ndens2,1)
*           Stop 10
*      End If
      FacR=Fact
      Call Read2_ns(rmo1,rmo2,
     &           FockI,FockA,
     &           T1,imem,Scr,Tmp2,
     &           T3,T4,
     &           DIR,DIL,
     &           DI,
     &           DAR,DAL,
     &           DA,
     &           rkappa,idsym,Sign,Facr,jSpin,
     &           lFA,lfi,lMo,
     &           CMO)
*      If (iDsym.eq.2) Then
*           Call RecPrt('FI',' ',FI,ndens2,1)
*           Call RecPrt('DI',' ',DI,ndens2,1)
*           Stop 10
*      End If
      Do iS=1,nsym
      js=ieor(idsym-1,is-1)+1
      If (nbas(is)*nbas(js).ne.0)
     &Call DGETMO(rkappa(ipmat(is,js)),nbas(is),nbas(is),nbas(js),
     &            K1(ipmat(js,is)),nbas(js))
      End Do
      Call DSCAL_(ndens2,-1.0d0,K1,1)
      DIR(:)=0.0d0
      DIL(:)=0.0d0
      If (imethod.eq.2) Then
         DAR(:)=0.0d0
         DAL(:)=0.0d0
      End If
      Call Read2_ns(rdum,rdum,
     &           FI1,FA1,
     &           T1,imem,Scr,Tmp2,
     &           T3,T4,
     &           DIR,DIL,
     &           DI,
     &           DAR,DAL,
     &           DA,
     &           K1,idsym,Sign,Facr,jSpin,
     &           lFA,lfi,.false.,
     &           CMO)
*      If (iDsym.eq.2)
*     &      Call RecPrt('FI1',' ',FI1,ndens2,1)
*
*      Calculate contribution from uncontracted indexes.
*
      Do iS=1,nSym
       jS=iEOr(iS-1,iDSym-1)+1
       If (nBas(iS)*nBas(jS).ne.0) Then
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            FIMO(ipCM(iS)),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            FIMO(ipCM(jS)),nBas(jS),
     &            One,FockI(ipMat(iS,jS)),nBas(is))
       If (iMethod.eq.2) Then
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            FAMO(ipCM(iS)),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            FAMO(ipCM(jS)),nBas(jS),
     &            One,FockA(ipMat(iS,jS)),nBas(is))
       End If
       End If
      End Do
*
      Call mma_deallocate(DA)
      Call mma_deallocate(DAR)
      Call mma_deallocate(DAL)
      Call mma_deallocate(FA1)
      Call mma_deallocate(FI)
      Call mma_deallocate(FI1)
      Call mma_deallocate(K1)
      Call mma_deallocate(DIR)
      Call mma_deallocate(DI)
      Call mma_deallocate(DIL)
      Call mma_deallocate(T4)
      Call mma_deallocate(T3)
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Scr)
      Call mma_deallocate(T1)

*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nF)
      End
