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
*                                                                      *
************************************************************************
*                                                                      *
      use Arrays, only: CMO, G1t, FAMO, FIMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One, Two
      use MCLR_Data, only: nDens2, nMBA, ipCM, ipMat, nA, nCMO
      use input_mclr, only: nSym,nAsh,nIsh,nBas,iMethod
      Implicit None
      Real*8 rKappa(nDens2),rMO1(nMBA),rMO2(nMBA),FockI(nDens2),
     &       FockA(nDens2)
      Integer nF,iDSym,jSpin
      Real*8 sign,Fact

      Logical lFI,lFA,lMo
      Real*8 rdum(1)
      Real*8, Allocatable:: T1(:), Tmp2(:), T3(:), T4(:), DIL(:),
     &                      DI(:), DIR(:), FI(:), FI1(:),
     &                      K1(:), DAL(:), DAR(:), DA(:), FA1(:)
      Integer nDens22, iAM, iBM, iMem, iS, iB, ip, jB, iA, jA, ip2, jS
      Real*8 FacR
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      SubRoutine Read2_ns(rMO1,rMO2,FockI,FockA,
     &                  Temp1,nTemp,Temp2,Temp3,Temp4,
     &                  DI13,DI24,DI,
     &                  DA13,DA24,DA,
     &                  rkappa,idsym,
     &                  Signa,Fact,jSpin,lfat,lfit,lMOt,CMO)
      use MCLR_Data, only: nMBA, nDens2, nCMO
      Implicit None
      real*8 rmo1(nMBA),rmo2(nMBA),FockI(nDens2),FockA(nDens2)
      Integer nTemp
      real*8 Temp1(ntemp),Temp2(nDens2),Temp3(nDens2),Temp4(nDens2),
     &       DI13(nDens2),DI24(nDens2),DI(nCMO),
     &       DA13(nDens2),DA24(nDens2),DA(nCMO),
     &       rkappa(nDens2)
      Integer iDSym
      real*8 Signa,Fact
      Integer jSpin
      Logical lFAt,lFIT,lmot
      real*8 CMO(nCMO)

      End SubRoutine Read2_ns
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      Integer i,j,itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
************************************************************************
*

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
      Call mma_allocate(FI1,ndens2,Label='FI1')
      Call mma_allocate(K1,ndens2,Label='K1')

      FockI(:)=Zero
      FockA(:)=Zero
      FI(:)   =Zero
      FI1(:)  =Zero
      K1(:)   =Zero
      DIR(:)  =Zero
      DIL(:)  =Zero
      DI(:)   =Zero
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
      FA1(:)=Zero
      DAL(:)=Zero
      DAR(:)=Zero
      DA(:) =Zero
c THIS IS THE ENTIRE DENSITY FOR MULTICONF
      Do iS=1,nSym
         Do iB=1,nIsh(iS)
            ip=ipCM(iS)+(ib-1)*nBas(is)+ib-1
            DI(ip)=Two
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
     &           T1,imem,Tmp2,
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
      Call DSCAL_(ndens2,-One,K1,1)
      DIR(:)=Zero
      DIL(:)=Zero
      If (imethod.eq.2) Then
         DAR(:)=Zero
         DAL(:)=Zero
      End If
      Call Read2_ns(rdum,rdum,
     &           FI1,FA1,
     &           T1,imem,Tmp2,
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
      Call mma_deallocate(T1)

*                                                                      *
************************************************************************
*                                                                      *
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nF)
      End SubRoutine r2elint_ns
