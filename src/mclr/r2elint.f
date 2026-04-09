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
      SubRoutine r2elint(rKappa,rMO1,rmo2,FockI,FockA,nF,
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
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One, Two
      use MCLR_Data, only: nDens2, nMBA, ipCM, ipMat, nA, nCMO
      use input_mclr, only: nSym,nAsh,nIsh,nBas,nOrb,iMethod,CasInt
      Implicit None
      Real*8 rKappa(nDens2),rMO1(nMba),rmo2(*),FockI(nDens2),
     &       FockA(nDens2)
      Integer nF,iDSym,jSpin
      Real*8 sign,Fact

      Logical lFI,lFA,lMo
      Real*8, Allocatable:: T1(:), Tmp2(:), T3(:), T4(:), DIL(:),
     &                      DI(:), DIR(:), FI(:),
     &                      DAL(:), DAR(:), DA(:), FA(:)
      Integer nDens22, iAM, iBM, iMem, iS, iB, ip, jB, iA, ip2, jS, jA
      Real*8 FacR

*
      integer i,j,itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
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

      FockI(:)=Zero
      FockA(:)=Zero
      FI(:)   =Zero
      DI(:)   =Zero
      DIL(:)  =Zero
      DIR(:)  =Zero
      lFI=.true.
      lFa=.false.
      lMo=.false.
      If (iMethod.eq.2) Then
         Call mma_allocate(DAL,nDens2,Label='DAL')
         Call mma_allocate(DAR,nDens2,Label='DAR')
         Call mma_allocate(DA,nCMO,Label='DA')
         Call mma_allocate(FA,nDens2,Label='FA')
         lFa=.true.
         lMo=.true.
      Else
         Call mma_allocate(DAL,1,Label='DAL')
         Call mma_allocate(DAR,1,Label='DAR')
         Call mma_allocate(DA,1,Label='DA')
         Call mma_allocate(FA,1,Label='FA')
      End If
      FA(:) =Zero
      DA(:) =Zero
      DAL(:)=Zero
      DAR(:)=Zero

      Do iS=1,nSym
         Do iB=1,nIsh(iS)
            ip=ipCM(iS)+(ib-1)*nOrb(is)+ib-1
            DI(ip)=Two
         End Do
      End Do
      If (iMethod.eq.2) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1
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
      FacR=Fact
      Call Read2_2(rmo1,rmo2,
     &         FockI,FockA,
     &         T1,imem,
     &         Tmp2,
     &         T3,T4,nDens22,
     &         DIR,DIL,
     &         DI,
     &         DAR,DAL,
     &         DA,
     &         rkappa,idsym,Sign,Facr,jSpin,
     &         lFA,lfi,lMo,
     &         CMO)
*
*       call recprt('1 FockI ','',FockI,nDens2,1)
*       call recprt('1 FockA ','',FockA,nDens2,1)
*
*      Calculate contribution from uncontracted indexes.
*
      Do iS=1,nSym
       jS=iEOr(iS-1,iDSym-1)+1
       If (nOrb(iS)*nOrb(jS).ne.0) Then
       If (.not.CASINT)
     & Call DGEMM_('T','N',
     &             nOrb(iS),nOrb(jS),nBas(iS),
     &             1.0d0,CMO(ipCM(iS)),nBas(iS),
     &             FI(ipMat(iS,jS)),nBas(iS),
     &             Zero,FockI(ipMat(iS,jS)),nOrb(iS))
       Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sign*Facr,
     &            FIMO(ipCM(iS)),nOrb(is),
     &            rkappa(ipMat(iS,jS)),nOrb(iS),
     &            One,FockI(ipMat(iS,jS)),nOrb(iS))
       Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nOrb(is),
     &            FIMO(ipCM(jS)),nOrb(jS),
     &            One,FockI(ipMat(iS,jS)),nOrb(is))
       If (iMethod.eq.2) Then
       If (.not.CASINT)
     &   Call DGEMM_('T','N',
     &               nOrb(iS),nOrb(jS),nBas(iS),
     &               1.0d0,CMO(ipCM(iS)),nBas(iS),
     &               FA(ipMat(iS,jS)),nBas(iS),
     &               Zero,FockA(ipMat(iS,jS)),nOrb(iS))
         Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Sign*Facr,
     &            FAMO(ipCM(iS)),nOrb(is),
     &            rkappa(ipMat(iS,jS)),nOrb(iS),
     &            One,FockA(ipMat(iS,jS)),nOrb(iS))
         Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Facr,
     &            rkappa(ipMat(iS,jS)),nOrb(is),
     &            FAMO(ipCM(jS)),nOrb(jS),
     &            One,FockA(ipMat(iS,jS)),nOrb(is))
       End If
       End If
      End Do

!       call recprt('2 FockI ','',FockI,nDens2,1)
!       call recprt('2 FockA ','',FockA,nDens2,1)

      Call mma_deallocate(DA)
      Call mma_deallocate(DAR)
      Call mma_deallocate(DAL)
      Call mma_deallocate(FA)
      Call mma_deallocate(FI)
      Call mma_deallocate(DIR)
      Call mma_deallocate(DI)
      Call mma_deallocate(DIL)
      Call mma_deallocate(T4)
      Call mma_deallocate(T3)
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(T1)

c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(nF)

      End SubRoutine r2elint
