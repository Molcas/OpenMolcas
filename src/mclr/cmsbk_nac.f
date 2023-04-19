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
* Copyright (C) 2021, Paul B Calio                                     *
************************************************************************
* ****************************************************************
* history:                                                       *
* Based on cmsbk.f from Jie J. Bao &&                            *
* Additional work from  rhs_nac.f                                *
* ****************************************************************
      Subroutine Calcbk_CMSNAC(bk,R,nTri,GDMat,zX)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "sa.fh"

******Output
      Real*8,DIMENSION(nDens2)::bk
******Input
      Real*8,DIMENSION(nRoots**2)::R
      INTEGER nTri
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
      Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
******Auxiliaries
      Real*8,DIMENSION(:),Allocatable::FOccMO,P2MOt
      INTEGER nP2,nG1
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      ng1=itri(ntash,ntash)
      nP2=itri(ng1,ng1)
      CALL mma_allocate(FOccMO,nDens2)
      CALL mma_allocate(P2MOt,nP2)
      CALL FZero(bk,nDens2)
      CALL GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,nP2)
      CALL GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)
      CALL GetPDFTFock_NAC(bk)
      CALL PutCMSFockOcc(FOccMO,nTri)

      CALL mma_deallocate(FOccMO)
      CALL mma_deallocate(P2MOt)

      RETURN
      end subroutine
******************************************************

******************************************************
      Subroutine GetPDFTFock_NAC(bk)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "sa.fh"
******Output
      Real*8,DIMENSION(nDens2)::bk
******Input
******Auxiliaries
      Real*8,DIMENSION(:),Allocatable::T,FT99,bktmp
      INTEGER IS,JS
      CALL mma_allocate(FT99,nDens2)
      CALL mma_allocate(bktmp,nDens2)
      CALL mma_allocate(T,nDens2)
      CALL Get_DArray('FxyMS           ',FT99 ,nDens2)
      CALL dcopy_(nDens2,FT99,1,T,1)

      DO IS=1,nSym
         jS=iEOR(iS-1,0)+1
         If (nBas(is)*nBas(jS).ne.0) then
           Call DGeSub(T(ipMat(iS,jS)),nBas(iS),'N',
     &                 T(ipMat(jS,iS)),nBas(jS),'T',
     &                 bktmp(ipMat(iS,jS)),nBas(iS),
     &                 nBas(iS),nBas(jS))
         End If
      END DO
      CALL daxpy_(nDens2,-2.0d0,bktmp,1,bk,1)
      CALL mma_deallocate(T)
      CALL mma_deallocate(FT99)
      CALL mma_deallocate(bktmp)
      RETURN
      end subroutine
******************************************************
******************************************************
      Subroutine GetWFFock_NAC(FOccMO,bk,R,nTri,P2MOt,NG2)
******Partially readpated from rhs_sa.f
      use stdalloc, only : mma_allocate, mma_deallocate
      use ipPage, only: W
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "sa.fh"
#include "real.fh"
******Input
      Real*8,DIMENSION(nRoots**2)::R
      INTEGER nTri,NG2
******Output
      Real*8,DIMENSION(nDens2)::FOccMO
      Real*8,DIMENSION(nDens2)::bk
      Real*8,DIMENSION(nG2)::P2MOt
******Auxiliaries
      Real*8,DIMENSION(:),Allocatable::FinCI
*     FinCI: CI Vectors in final CMS state basis
      Real*8,DIMENSION(1)::rdum
      Real*8,DIMENSION(:),Allocatable::Fock,T,G1r,G2r,G2rt,
     & CIL,CIR,G1q,G2q,G1qs,G2qs,G1m
      Real*8,DIMENSION(:),Allocatable::DMatAO,D5,D6
      INTEGER I,K,NCSFs
      Real*8 Fact
      INTEGER iB,jB,kB,lB,iDkl,iRijkl

      INTEGER IJ, KL, IJKL, IJ2, KL2
      Real*8 factor
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
       ng1=itri(ntash,ntash)
       ng2=itri(ng1,ng1)

       Call mma_allocate(FinCI,nconf1*nroots,Label='FinCI')
       Call mma_allocate(Fock,ndens2,Label='Fock')
       Call mma_allocate(T,ndens2,Label='T')
       Call mma_allocate(G1q,ng1,Label='G1q')
       Call mma_allocate(G1m,ng1,Label='G1m')
       Call mma_allocate(G2q,ng2,Label='G2q')
       Call mma_allocate(G1r,ntash**2,Label='G1r')
       Call mma_allocate(G2r,itri(ntash**2,ntash**2),Label='G2r')
       Call mma_allocate(G2rt,itri(ntash**2,ntash**2),Label='G2rt')
*******Rotate CI vectors back to those for reference states
       NCSFs=NCSF(state_sym)
       CALL DGEMM_('n','n',NCSFS,nRoots,nRoots,1.0d0,W(ipCI)%Vec,
     &             NCSFs,R,nRoots,0.0d0,FinCI,nCSFs)
       nConfL=Max(ncsf(state_sym),nint(xispsm(state_sym,1)))
       nConfR=Max(ncsf(state_sym),nint(xispsm(state_sym,1)))

       Call mma_allocate(CIL,nConfL)
       Call mma_allocate(CIR,nConfR)

       I=NACstates(1)
       J=NACstates(2)
       Call CSF2SD(FinCI(1+(J-1)*NCSFs),CIL,state_sym)
       Call CSF2SD(FinCI(1+(I-1)*NCSFs),CIR,state_sym)
       Call Densi2(2,G1r,G2rt,CIL,CIR,0,0,0,ntash**2,
     &              itri(ntash**2,ntash**2))

*Copied from rhs_nac.f
       ij=0
       Do iB=0,nnA-1
         Do jB=0,iB-1
           ij=ij+1
           G1q(ij)=(G1r(1+iB*ntAsh+jB)+
     &                    G1r(1+jB*ntAsh+iB))*Half
*          Note that the order of subtraction depends on how the matrix
*          will be used when contracting with derivative integrals
*          This is found to give the correct results:
           G1m(ij)=(G1r(1+jB*ntAsh+iB)-
     &                  G1r(1+iB*ntAsh+jB))*Half
         End Do
         ij=ij+1
         G1q(ij)=G1r(1+iB*ntAsh+iB)
         G1m(ij)=Zero
       End Do

       Do iB=1,ntash
        Do jB=1,ntash
        G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
        End Do
       End Do

       Do iB=1,ntAsh**2
         jB=itri(iB,iB)
         G2rt(jB)=Half*G2rt(jB)
       End Do
       Do iB=0,ntAsh-1
         Do jB=0,iB-1
           ij=iB*(iB+1)/2+jB
           Do kB=0,ntAsh-1
             Do lB=0,kB
               kl=kB*(kB+1)/2+lB
               If (ij.ge.kl) Then
                 factor=Quart
                 If (ij.eq.kl) factor=Half
                 ijkl=ij*(ij+1)/2+kl
                 ij2=iB*ntAsh+jB
                 kl2=kB*ntAsh+lB
                 G2q(1+ijkl)=factor*G2rt(1+ij2*(ij2+1)/2+kl2)
                 ij2=Max(jB*ntAsh+iB,lB*ntAsh+kB)
                 kl2=Min(jB*ntAsh+iB,lB*ntAsh+kB)
                 G2q(1+ijkl)=G2q(1+ijkl)+
     &                           factor*G2rt(1+ij2*(ij2+1)/2+kl2)
                 If (kB.ne.lB) Then
                   ij2=iB*ntAsh+jB
                   kl2=lB*ntAsh+kB
                   G2q(1+ijkl)=G2q(1+ijkl)+
     &                             factor*G2rt(1+ij2*(ij2+1)/2+kl2)
                   If (ij.ne.kl) Then
                     ij2=Max(jB*ntAsh+iB,kB*ntAsh+lB)
                     kl2=Min(jB*ntAsh+iB,kB*ntAsh+lB)
                     G2q(1+ijkl)=G2q(1+ijkl)+
     &                              factor*G2rt(1+ij2*(ij2+1)/2+kl2)
                   End If
                 End If
               End If
             End Do
           End Do
         End Do
         ij=iB*(iB+1)/2+iB
         Do kB=0,ntAsh-1
           Do lB=0,kB
             kl=kB*(kB+1)/2+lB
             If (ij.ge.kl) Then
               factor=Half
               If (ij.eq.kl) factor=One
               ijkl=ij*(ij+1)/2+kl
               ij2=iB*ntAsh+iB
               kl2=kB*ntAsh+lB
               G2q(1+ijkl)=factor*G2rt(1+ij2*(ij2+1)/2+kl2)
               If (kB.ne.lB) Then
                 kl2=lB*ntAsh+kB
                 G2q(1+ijkl)=G2q(1+ijkl)+
     &                           factor*G2rt(1+ij2*(ij2+1)/2+kl2)
               End If
             End If
           End Do
         End Do
       End Do

      Do iB=1,ntAsh
        Do jB=1,ntAsh
          ij=iTri(iB,jB)
          ij2=ntAsh*(iB-1)+jB
          Do kB=1,ntAsh
            Do lB=1,ntAsh
              kl=iTri(kB,lB)
              kl2=ntAsh*(kB-1)+lB
              factor=One
              If (ij.ge.kl .and. kB.eq.lB) factor=Two
              If (ij.lt.kl .and. iB.eq.jB) factor=Two
              ijkl=iTri(ij,kl)
              ijkl2=iTri(ij2,kl2)
              G2r(ijkl2)=factor*G2q(ijkl)
            End Do
          End Do
        End Do
      End Do

       Call FockGen(0.0d0,G1r,G2r,FOccMO,bk,1)

*******D1MOt: CMS-PDFT 1RDM for computing 1-electron gradient
       Call Put_DArray('D1MOt           ',G1q,ng1)

       iRC=0
       LuDens=20
       Call DaName(LuDens,'MCLRDENS')
       Call dDaFile(LuDens,1,G1m,ng1,iRC)
       Call DaClos(LuDens)

       Do iB=1,ntash
        Do jB=1,ntash
         iDij=iTri(ib,jB)
         iRij=jb+(ib-1)*ntash
         Do kB=1,ntash
          Do lB=1,ntash
           iDkl=iTri(kB,lB)
           iRkl=lb+(kb-1)*ntash
           fact=One
           if(iDij.ge.iDkl .and. kB.eq.lB) fact=0.5d0
           if(iDij.lt.iDkl .and. iB.eq.jB) fact=0.5d0
           iijkl=itri(iDij,iDkl)
           iRijkl=itri(iRij,iRkl)
           G2q(iijkl)=Fact*G2r(iRijkl)
          End Do
         End Do
        End Do
       End Do

       Call Get_dArray_chk('P2MOt',P2MOt,ng2)
       Call DaXpY_(ng2,1.0d0,G2q,1,P2MOt,1)

*******Done with the info from CMS final state

*******Doing some computation for computing non-active-active 2RDM in
*******integral_util/prepp.f
       Call mma_allocate(D5,nTri)
       Call mma_allocate(D6,nTri)
*******D5: Used in ptrans_sa when isym==jsym (PDFT parts cancel WF
*******    parts for intermediate states)
*******D6: Used in ptrans_sa when isym.ne.jsym (sum of inactive parts of
*******intermediate-state 1RDMs cancels that of the final state)
       Call mma_allocate(DMatAO,nTri)
       CALL Get_DArray('MSPDFTD6        ',D6,nTri)
       CALL GetDMatAO(G1q,DMatAO,ng1,nTri)
       CALL DaXpY_(nTri,1.0d0,DMatAO,1,D6,1)
       CALL DCopy_(nTri,DMatAO,1,D5,1)
       CALL Put_DArray('MSPDFTD5        ',D5,nTri)
       CALL Put_DArray('MSPDFTD6        ',D6,nTri)
       Call mma_deallocate(D5)
       Call mma_deallocate(D6)
       Call mma_deallocate(DMatAO)
*******Beginning of the info for CMS intermediate states
       jdisk=itoc(3)
       Call mma_allocate(G1qs,ng1*nRoots)
       Call mma_allocate(G2qs,ng2*nRoots)
       DO K=1,nRoots
        Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)
        Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
        Call dDaFile(LUJOB ,2,G2q,Ng2,jDisk)
        Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
        Call dcopy_(ng1,G1q,1,G1qs((K-1)*ng1+1),1)
        Call dcopy_(ng2,G2q,1,G2qs((K-1)*ng2+1),1)
        Call mma_allocate(DMatAO,ntri)
        CALL GetDMatAO(G1q,DMatAO,ng1,nTri)
        Call mma_deallocate(DMatAO)

        Do iB=1,ntash
         Do jB=1,ntash
         G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
         End Do
        End Do

        Do iB=1,ntash
         Do jB=1,ntash
          iDij=iTri(ib,jB)
          iRij=jb+(ib-1)*ntash
          Do kB=1,ntash
           Do lB=1,ntash
            iDkl=iTri(kB,lB)
            iRkl=lb+(kb-1)*ntash
            fact=One
            if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
            if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
            iijkl=itri(iDij,iDkl)
            iRijkl=itri(iRij,iRkl)
            G2r(iRijkl)=Fact*G2q(iijkl)
           End Do
          End Do
         End Do
        End Do

        Call FockGen(0.0d0,G1r,G2r,T,Fock,1)
        CALL Daxpy_(nDens2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),
     &    Fock,1,bk,1)
        CALL Daxpy_(nDens2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),
     &    T,1,FOccMO,1)
        Call DaXpY_(ng2,-R((I-1)*nRoots+K)*R((J-1)*nRoots+K),
     &    G2q,1,P2MOt,1)
       END DO
       Call Put_DArray('D1INTER         ',G1qs,ng1*nRoots)
       Call Put_DArray('P2INTER         ',G2qs,ng2*nRoots)
       Call mma_deallocate(G1qs)
       Call mma_deallocate(G2qs)
       Call mma_deallocate(Fock)
       Call mma_deallocate(T)
       Call mma_deallocate(G1r)
       Call mma_deallocate(G1m)
       Call mma_deallocate(G2r)
       Call mma_deallocate(G2rt)
       Call mma_deallocate(G1q)
       Call mma_deallocate(G2q)
       Call mma_deallocate(CIL)
       Call mma_deallocate(CIR)
       Call mma_deallocate(FinCI)
       RETURN
       End Subroutine
