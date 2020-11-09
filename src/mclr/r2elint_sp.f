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
      SubRoutine r2elint_sp(rKappa,rMO1,rmo2,FockI,FockA,nF,
     &                   iDSym,sign,Fact,jspin,D,FA)
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
      use Arrays, only: CMO, FIMO
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "spin.fh"
      Logical lFI,lFA,lMo
      Parameter ( One = 1.0d0 )
      Real*8 rKappa(nDens2),rMO1(nMba),rmo2(*),FockI(nDens2),
     &       FockA(nDens2),D(*),FA(*)
      Real*8, Allocatable:: T1(:), Tmp2(:), T3(:), T4(:), DIL(:),
     &                      DI(:), DIR(:), FI(:),
     &                      FA2(:), DAR(:), DA(:), DAL(:)
*
      irec(i,j)=i+nna*(j-1)
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
      Call mma_allocate(Tmp2,nDens22,Label='Tmp2')
      Call mma_allocate(T3,nDens22,Label='T3')
      Call mma_allocate(T4,nDens22,Label='T4')
      Call mma_allocate(DIL,nDens2,Label='DIL')
      Call mma_allocate(DI,nCMO,Label='DI')
      Call mma_allocate(DIR,nDens2,Label='DIR')
      Call mma_allocate(FI,ndens2,Label='FI')

      FockI(:)=0.0d0
      FockA(:)=0.0d0
      FI(:)  =0.0D0
      DIL(:)  =0.0D0
      DIR(:)  =0.0D0

      lFI=.true.
      lFa=.false.
      lMo=.false.
      If (iMethod.eq.2) Then
         Call mma_allocate(DAL,nDens2,Label='DAL')
         Call mma_allocate(DAR,nDens2,Label='DAR')
         Call mma_allocate(DA,nCMO,Label='DA')
         Call mma_allocate(FA2,nDens2,Label='FA2')
         lFa=.true.
         lMo=.true.
      Else
         Call mma_allocate(DAL,1,Label='DAL')
         Call mma_allocate(DAR,1,Label='DAR')
         Call mma_allocate(DA,1,Label='DA')
         Call mma_allocate(FA2,1,Label='FA2')
      End If

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
            DA(ip)=D(irec(iA,jA))
         End Do
        End Do
       End Do
      End If

*
*     Construct {kappa,(ij|kl)} and the fock matrix contribution
*     from one index tranformation of contracted indexes
*
      FacR=Fact
      Call Read2_2(rmo1,rmo2,
     &           FI,FA2,
     &           T1,imem,Tmp2,
     &           T3,T4,
     &           DIR,DIL,
     &           DI,
     &           DAR,DAL,
     &           DA,
     &           rkappa,idsym,Sign,Facr,jSpin,
     &           lFA,lfi,lMo,
     &           CMO)
*
*      Calculate contribution from uncontracted indexes.
*
      Do iS=1,nSym
       jS=iEOr(iS-1,iDSym-1)+1
       If (nBas(iS)*nBas(jS).ne.0) Then
       Call DGEMM_('T','N',
     &             nBas(iS),nBas(jS),nBas(iS),
     &             1.0d0,CMO(ipCM(iS)),nBas(iS),
     &             FI(ipMat(iS,jS)),nBas(iS),
     &             0.0d0,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            FIMO(ipCM(iS)),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockI(ipMat(iS,jS)),nBas(iS))
       Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            FIMO(ipCM(jS)),nBas(jS),
     &            One,FockI(ipMat(iS,jS)),nBas(is))
       If (iMethod.eq.2) Then
         Call DGEMM_('T','N',
     &               nBas(iS),nBas(jS),nBas(iS),
     &               1.0d0,CMO(ipCM(iS)),nBas(iS),
     &               FA2(ipMat(iS,jS)),nBas(iS),
     &               0.0d0,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(iS),Sign*Facr,
     &            FA(ipCM(iS)),nBas(is),
     &            rkappa(ipMat(iS,jS)),nBas(iS),
     &            One,FockA(ipMat(iS,jS)),nBas(iS))
         Call DGEMM_('N','N',nBas(iS),nBas(jS),nBas(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nBas(is),
     &            FA(ipCM(jS)),nBas(jS),
     &            One,FockA(ipMat(iS,jS)),nBas(is))
       End If
       End If
      End Do

      Call mma_deallocate(FA2)
      Call mma_deallocate(DA)
      Call mma_deallocate(DAR)
      Call mma_deallocate(DAL)
      Call mma_deallocate(FI)
      Call mma_deallocate(DIR)
      Call mma_deallocate(DI)
      Call mma_deallocate(DIL)
      Call mma_deallocate(T4)
      Call mma_deallocate(T3)
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(T1)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nF)
      End
