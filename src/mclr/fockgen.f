************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       SubRoutine FockGen(d_0,rDens1,rdens2,Fock,FockOut,idSym)
************************************************************************
*                                                                      *
*   Constructs active Fock matrix and Q matrix                         *
*                                                                      *
*   Input: rkappa: Rotation matrix                                     *
*          idsym : symmetry of perturbation                            *
*                                                                      *
*                                                                      *
*   Output:MO     : MO integrals                                       *
*          Fock   : Fock matrix (one index transformed integrals)      *
*          MOtilde: MO (one index transformed integrals)               *
*                                                                      *
************************************************************************
      use Arrays, only: CMO, FIMO
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "standard_iounits.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
      Real*8 Fock(nDens2),FockOut(*), rDens2(*),rDens1(nna,nna)
      Real*8, Allocatable:: MO(:), Scr(:), G2x(:), Scr1(:,:), Ash(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(nDens2,[Zero],0,Fock,1)
*
      n1=0
      Do iS = 1, nSym
        n1=Max(n1,nBas(iS))
      End Do
      n2=n1**2

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
      end if

      If (newCho) Go to 15
      Call mma_allocate(MO,n2,Label='MO')
      Call mma_allocate(Scr,n2,Label='Scr')
*
      Do ipS=1,nSym
         Do kS=1,nSym
            Do iS=1,nSym
               jS=iEor(iEor(ipS-1,kS-1),iS-1)+1
*                                                                      *
************************************************************************
*                                                                      *
*              Coulomb term: F  =2(pk|ji)d
*                             kp          ij
*                                                                      *
************************************************************************
*                                                                      *
               If (iEOr(ipS-1,kS-1)+1.eq.iDsym .and.
     &             nBas(ipS)*nIsh(kS).gt.0           ) Then
                  Do iA=1,nAsh(iS)
                     iAA=iA+nIsh(iS)
                     Do jA=1,nAsh(jS)
                        jAA=jA+nIsh(jS)
*
                        Call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)
*
                        rD=rDens1(iA+nA(iS),jA+nA(jS))*Two
                        Call DaXpY_(nBas(ipS)*nIsh(kS),rd,
     &                              MO,1,Fock(ipMat(ipS,Ks)),1)
*
                     End Do
                  End Do
               End If
*                                                                      *
************************************************************************
*                                                                      *
*              Exchange term: F = -(pk|ji)d
*                              pl          kj
*                                                                      *
************************************************************************
*                                                                      *
               If (iEOr(ipS-1,iS-1)+1.eq.iDsym .and.
     &             nBas(ipS).gt.0                   ) Then
                  Do iA = 1, nIsh(iS)
                     ipF=ipMat(ipS,iS)+nBas(ipS)*(iA-1)
                     Do jA=1,nAsh(jS)
                        jAA=jA+nIsh(jS)
*
                        Call Coul(ipS,kS,iS,jS,iA,jAA,MO,Scr)
*
                        ipM=1+nIsh(kS)*nBas(ipS)
                        Do kA=1,nAsh(kS)
*
                           rd=rDens1(kA+nA(kS),jA+nA(jS))
                           Call DaXpY_(nBas(ipS),-rd,
     &                                MO(ipM),1,Fock(ipF),1)
                           ipM = ipM + nBas(ipS)
                        End Do
*
                     End Do
                  End Do
               End If
*                                                                      *
************************************************************************
*                                                                      *
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call CreQADD(Fock,rdens2,idsym,MO,Scr,n2)
      Call mma_deallocate(Scr)
      Call mma_deallocate(MO)
*
************************************************************************
*       new Cholesky code                                              *
************************************************************************
 15   Continue
      if (newCho) Then
        nVB=0
        nG2=0
        Do iSym=1,nSym
          nVB = nVB + nAsh(iSym)*nOrb(iSym)
          nAG2=0
          Do jSym=1,nSym
            kSym=iEOr(jsym-1,isym-1)+1
            nAG2=nAg2+nAsh(jSym)*nAsh(kSym)
          End Do
          nG2=nG2+nAG2**2
        End Do

*
**      Unfold 2-DM
*
        Call mma_allocate(G2x,nG2,Label='G2x')
        ipGx=0
        Do ijS=1,nSym
          Do iS=1,nSym
            jS=iEOR(is-1,ijS-1)+1
            Do kS=1,nSym
              lS=iEOR(kS-1,ijS-1)+1
              Do kAsh=1,nAsh(ks)
                Do lAsh=1,nAsh(ls)
c                 ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
                  ikl=nna*(lAsh+nA(lS)-1)+kAsh+nA(kS)
                  Do iAsh=1,nAsh(is)
                    Do jAsh=1,nAsh(js)
c                     iij =itri(iAsh+nA(is),jAsh+nA(jS))
                      iij=nna*(jAsh+nA(jS)-1)+iAsh+nA(iS)
                      ipGx=ipGx+1
                      G2x(ipGx)=rdens2(itri(iij,ikl))
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End Do
*
**      Get active CMO
*
        Call mma_Allocate(Ash,nVB,Label='Ash')
        ioff=0
        ioff1=0
        Do iS=1,nSym
          ioff2 = ioff + nOrb(iS)*nIsh(iS)
          Do iB=1,nAsh(iS)
            ioff3=ioff2+nOrb(iS)*(iB-1)
            call dcopy_(nOrb(iS),CMO(1+ioff3),1,Ash(ioff1+iB),nAsh(iS))
          End Do
          ioff=ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
          ioff1=ioff1+nAsh(iS)*nOrb(iS)
        End Do

*
        Call mma_allocate(Scr1,n2,2,Label='Scr1')
        Scr1(:,:)=Zero
        ipDA  =ip_of_work(rdens1(1,1))
        ipFock=ip_of_work(Fock(1))
        ipCMO =ip_of_work(CMO(1))
        ipG2x =ip_of_work(G2x(1))
        ipScr1=ip_of_work(Scr1(1,1))
        ipScr2=ip_of_work(Scr1(1,2))
        ipAsh =ip_of_work(Ash(1))
*
        Call cho_fock_mclr(ipDA,ipG2x,ipScr1,ipScr2,ipFock,
     &                    [ipAsh],ipCMO,nIsh,nAsh,LuAChoVec)
*
        Call mma_deallocate(Scr1)
        Call mma_deallocate(Ash)
        Call mma_deallocate(G2x)
      EndIf
*
************************************************************************
*       Common part                                                    *
************************************************************************
*
      Do iS=1,nSym
         If (nBas(iS).gt.0) Then
            jS=iEOr(is-1,iDSym-1)+1
            Do iA=1,nAsh(is)
               Do jA=1,nAsh(js)
                  rd=rDens1(iA+nA(iS),jA+nA(js))
                  ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
                  ip2=nBas(iS)*(nIsh(js)+jA-1) +ipmat(is,js)
                 Call DaXpY_(nBas(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
               End Do
            End Do
         End If
      End Do

*
      If (iDsym.eq.1) Then
         Do iS=1,nSym
            If (nBas(iS)*nIsh(iS).gt.0)
     &         Call DaXpY_(nBas(iS)*nIsh(is),Two*d_0,
     &                    FIMO(ipMat(is,is)),1,
     &                    Fock(ipMat(is,is)),1)
         End Do
      End If
*
      Do iS=1,nSym
         jS=iEOR(iS-1,idSym-1)+1
         If (nBas(is)*nBas(jS).ne.0)
     &      Call DGeSub(Fock(ipMat(iS,jS)),nBas(iS),'N',
     &                  Fock(ipMat(jS,iS)),nBas(jS),'T',
     &                  FockOut(ipMat(iS,jS)),nBas(iS),
     &                  nBas(iS),nBas(jS))
      End Do

*
*
      Call DScal_(nDens2,Two,FockOut,1)
      If (idSym.eq.1) Call AddGrad2(FockOut,idSym,d_0)

      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
