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
      SubRoutine RInt_Generic(rkappa,rmos,rmoa,Fock,Q,Focki,Focka,
     &                        idsym,reco,jspin)
      use Arrays, only: CMO_Inv, CMO, G1t, FAMO
*
*                              ~
*     Constructs  F  = <0|[E  ,H]|0> ( + <0|[[E  , Kappa],H]|0> )
*                  pq       pq                 pq
*
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "glbbas_mclr.fh"
#include "standard_iounits.fh"
      Real*8 Fock(nDens2),focka(nDens2),rkappa(nDens2),
     &       Focki(ndens2),Q(ndens2),rMOs(*),rmoa(*)
      Integer ipAsh(2)
      Logical Fake_CMO2,DoAct
      Real*8, Allocatable:: MT1(:), MT2(:), MT3(:), QTemp(:), DI(:),
     &                      DLT(:), Dens2(:), DA(:), G2x(:),
     &                      CoulExch(:,:), CVa(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function
*
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Focki=',DDot_(nDens2,Focki,1,Focki,1)
      Write (LuWr,*) 'Focka=',DDot_(nDens2,Focka,1,Focka,1)
#endif
*
      Fact=-One
      call dcopy_(ndens2,[0.0d0],0,Fock,1)
*
      If (.not.newCho) Then
        Call mma_allocate(MT1,nmba,Label='MT1')
        Call mma_allocate(MT2,nmba,Label='MT2')
        MT1(:)=Zero
        MT2(:)=Zero

        Call R2ElInt(rkappa,MT1,MT2,focki,focka,
     &               nDens2,idSym,ReCo,Fact,jspin)
#ifdef _DEBUGPRINT_
        Write (LuWr,*) 'MT1=',DDot_(nmba,MT1,1,MT1,1)
        Write (LuWr,*) 'MT2=',DDot_(nmba,MT2,1,MT2,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
*
*     Q  = sum(jkl)=(pj|kl)d(ijkl)
*      pi
*
        If (iMethod.eq.2) Then
           Call mma_allocate(QTemp,ndens2,Label='QTemp')
           Call CreQ(Q,MT1,Work(ipG2tpp),idsym)
           Call CreQ(QTemp,MT2,Work(ipG2tmm),idsym)
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'Q=',DDot_(nDens2,Q,1,Q,1)
           Write (LuWr,*) 'QTemp=',DDot_(nDens2,QTemp,1,QTemp,1)
#endif
           call daxpy_(ndens2,One,QTemp,1,Q,1)
           Call mma_deallocate(QTemp)
        End If
*************************************************************************
*                                                                       *
*        Cholesky code                                                  *
*                                                                       *
*************************************************************************
      Else
        Fake_CMO2=.false.
        DoAct=.true.
*
        ipkappa  =ip_of_work(rkappa(1))
        ipFockI  =ip_of_work(FockI(1))
        ipFockA  =ip_of_work(FockA(1))
        ipMO1    =ip_of_work(rmos(1))
        ipQ      =ip_of_work(Q(1))
*
**      Form inactive density
*
        Call mma_Allocate(DI,nDens2,Label='DI')
        DI(:)=Zero

        Do iS=1,nsym
         Do iB=1,nIsh(iS)
          ip=ipCM(iS)+(ib-1)*nOrb(is)+ib-1
          DI(ip)=Two
         End Do
        End Do
*
**      Form AO 1-index transform inactive density
*
        Call mma_allocate(Dens2,nDens2,Label='Dens2')
        Call mma_allocate(DLT,nDens2,Label='DLT')
        Do iS=1,nSym
          If (nOrb(iS).ne.0) Then
            Do jS=1,nSym
              If (iEOr(iS-1,jS-1)+1.eq.idsym.and.nOrb(jS).ne.0) Then
                   Call DGEMM_('N','N',
     &                         nOrb(iS),nOrb(jS),nOrb(jS),
     &                         One,rkappa(ipMat(is,js)),nOrb(iS),
     &                             DI(ipCM(js)),nOrb(jS),
     &                         Zero,Dens2(ipMat(iS,jS)),nOrb(iS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         One,Dens2(ipMat(iS,jS)),nOrb(iS),
     &                             CMO(ipCM(is)),nOrb(iS),
     &                         Zero,DLT(ipMat(jS,iS)),nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),
     &                         One,DLT(ipMat(jS,iS)),nOrb(iS),
     &                             CMO(ipCM(js)),nOrb(jS),
     &                         Zero,Dens2(ipMat(iS,jS)),nOrb(jS))

                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         One,DI(ipCM(js)),
     &                         nOrb(iS),CMO(ipCM(is)),
     &                         nOrb(iS),Zero,
     &                         DLT(ipMat(jS,iS)),nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),One,
     &                         DLT(ipMat(jS,iS)),nOrb(iS),
     &                         CMO(ipCM(js)),nOrb(jS),Zero,
     &                         DI(ipCM(js)),nOrb(jS))
              EndIf
            End Do
          EndIf
        End Do
        call Fold_Mat(nSym,nOrb,Dens2,DLT)
        Call DScal_(ndens2,ReCo,DLT,1)
*
**      Form active density and MO coefficients
*
        nVB=0
        If (iMethod.eq.2) Then
          na2=0
          nAct=0
          nG2=0
          Do iSym=1,nSym
            nVB = nVB + nAsh(iSym)*nOrb(iSym)
            na2 = na2+nAsh(iSym)**2
            nact=nact+nAsh(iSym)
            nAG2=0
            Do jSym=1,nSym
              kSym=iEOr(jsym-1,isym-1)+1
              nAG2=nAg2+nAsh(jSym)*nAsh(kSym)
            End Do
            nG2=nG2+nAG2**2
          End Do
          nAtri=nact*(nact+1)/2
          nAtri=nAtri*(nAtri+1)/2
          Call mma_allocate(CVa,nVB,2,Label='CVa')
          Call mma_allocate(DA,na2,Label='DA')
*
          ioff=0
          ioff1=1
          ioff4=1
          ioffD=0
          Do iS=1,nSym
            ioff2 = ioff + nOrb(iS)*nIsh(iS)
            ioff5 = ioff4+ nOrb(iS)*nIsh(iS)
            Do iB=1,nAsh(iS)
              ioff3=ioff2+nOrb(iS)*(iB-1)
              call dcopy_(nOrb(iS),CMO(1+ioff3),1,
     &                  CVa(ioff1-1+iB,1),nAsh(iS))
             Do jB=1,nAsh(iS)
              ip=ioffD+ib+(jB-1)*nAsh(is)
              iA=nA(is)+ib
              jA=nA(is)+jb
              ip2=itri(iA,jA)
              DA(ip)=G1t(ip2)
             End Do
            End Do
*MGD to check
            Call DGEMM_('T','T',nAsh(iS),nOrb(iS),nOrb(iS),1.0d0,
     &                  rkappa(ioff5),nOrb(iS),
     &                  CMO(1+ioff),nOrb(iS),
     &                  0.0d0,CVa(ioff1,2),nAsh(iS))
            ioff=ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
            ioff1=ioff1+nAsh(iS)*nOrb(iS)
            ioffD=ioffD+nAsh(iS)**2
            ioff4=ioff4+nOrb(iS)**2
          End Do
*
**      Expand 2-body density matrix
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
                    ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
                    Do iAsh=1,nAsh(is)
                      Do jAsh=1,nAsh(js)
                        iij =itri(iAsh+nA(is),jAsh+nA(jS))
                        ipG=ipG2+itri(iij,ikl)-1
                        ipGx=ipGx+1
                        G2x(ipGx)=Work(ipG)
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        End If
*
**      Allocate temp arrays and zero Fock matrices
*
        Call mma_allocate(CoulExch,nDens2,4,Label='CoulExch')
        CoulExch(:,:)=Zero

        call dcopy_(ndens2,[0.0d0],0,FockA,1)
        call dcopy_(ndens2,[0.0d0],0,FockI,1)
        call dcopy_(nATri,[0.0d0],0,Work(ipMO1),1)
*
**      Compute the whole thing
*
        iread=2 ! Asks to read the half-transformed Cho vectors
        ip_CMO_Inv = ip_of_work(CMO_Inv(1))
        ipCMO      = ip_of_work(CMO(1))
        ipDI       = ip_of_Work(DI(1))
        ipDLT      = ip_of_Work(DLT(1))
        ipDA       = ip_of_Work(DA(1))
        ipG2x      = ip_of_Work(G2x(1))
        ipJI       = ip_of_Work(CoulExch(1,1))
        ipKI       = ip_of_Work(CoulExch(1,2))
        ipJA       = ip_of_Work(CoulExch(1,3))
        ipKA       = ip_of_Work(CoulExch(1,4))
        ipAsh(1)   = ip_of_Work(CVa(1,1))
        ipAsh(2)   = ip_of_Work(CVa(1,2))
        Call CHO_LK_MCLR(ipDLT,ipDI,ipDA,ipG2x,ipkappa,
     &                   ipJI,ipKI,ipJA,ipKA,ipFockI,ipFockA,
     &                   ipMO1,ipQ,ipAsh,ipCMO,ip_CMO_inv,
     &                   nIsh, nAsh,nIsh,DoAct,Fake_CMO2,
     &                   LuAChoVec,LuIChoVec,iread)
*
**      Calculate contribution from uncontracted indexes
*
        Do iS=1,nSym
          jS=iEOr(iS-1,iDSym-1)+1
          If (nOrb(iS)*nOrb(jS).ne.0) Then
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Reco*Fact,
     &                  Work(ipFIMO+ipCM(iS)-1),nOrb(is),
     &                  rkappa(ipMat(iS,jS)),nOrb(iS),
     &                  One,FockI(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,
     &                  rkappa(ipMat(iS,jS)),nOrb(is),
     &                  Work(ipFIMO+ipCM(jS)-1),nOrb(jS),
     &                  One,FockI(ipMat(iS,jS)),nOrb(is))
            If (iMethod.eq.2) Then
               Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),
     &                     Reco*Fact,
     &                     FAMO(ipCM(iS)),nOrb(is),
     &                     rkappa(ipMat(iS,jS)),nOrb(iS),
     &                     One,FockA(ipMat(iS,jS)),nOrb(iS))
               Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,
     &                     rkappa(ipMat(iS,jS)),nOrb(is),
     &                     FAMO(ipCM(jS)),nOrb(jS),
     &                     One,FockA(ipMat(iS,jS)),nOrb(is))
            End If
          End If
        End Do
*
**      Deallocate
*
        Call mma_deallocate(Dens2)
        Call mma_deallocate(CoulExch)
        If (iMethod.eq.2) Then
          Call mma_deallocate(Dens2)
          Call mma_deallocate(CVa)
          Call mma_deallocate(Dens2)
        EndIf
        Call mma_deallocate(DLT)
        Call mma_deallocate(DI)
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*
      Do iS=1,nSym
         jS=iEOr(iS-1,idsym-1)+1
*
*                I    A
*        F  =2( F  + F  )
*         pi     pi   pi
*
         If (nIsh(iS)*nOrb(jS).gt.0) Then
            Call DaXpY_(nIsh(iS)*nOrb(jS),2.0d0,
     &                 Focki(ipMat(js,is)),1,
     &                 Fock(ipMat(js,is)),1)
            If (iMethod.eq.2)  Call DaXpY_(nIsh(iS)*nOrb(jS),2.0d0,
     &                                    Focka(ipMat(js,is)),1,
     &                                    Fock(ipMat(js,is)),1)
         End If
*
         If (nOrb(iS).gt.0) Then
            Do iAsh=1,nAsh(jS)
               Do jAsh=1,nAsh(js)
                  ipF=ipMat(js,is)+nIsh(js)+jAsh-1
                  ipFI=ipMat(is,js)+(nIsh(js)+iAsh-1)*nOrb(is)
                  Dij=G1t(itri(iash+nA(js),jAsh+nA(js)))

*                I
*        F  = F - F  D
*         ap   ap  ap ba
*
                  ipF= ipMat(is,js)+(Nish(js)+iAsh-1)*nOrb(is)
                  ipFI=ipMat(is,js)+(Nish(js)+jAsh-1)*nOrb(is)
*
*                I
*        F  = F + F  D
*         pa   pa  pb ab
*
                  Call DaXpY_(nOrb(is),Dij,
     &                       Focki(ipFI),1,
     &                       Fock(ipF),1)
               End Do
            End Do
         End If
*
*
*      F  = F  + Q
*       pa   pa   pa
*
         If (nAsh(iS)*nOrb(jS).gt.0)
     &      Call DaXpY_(nAsh(is)*nOrb(js),1.0d0,
     &                 Q(ipMatba(js,is)),1,
     &                 Fock(ipMat(js,is)+nOrb(js)*nIsh(is)),1)
*
*
*      F  = F  - Q
*       ap   ap   ap
*
      End Do
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif
*
      Call DYAX(ndens2,2.0d0,Fock,1,Focka,1)
      Do iS=1,nSym
        js=iEOR(is-1,idsym-1)+1
        If (nOrb(is)*nOrb(js).ne.0)
     &  Call DGESUB(Focka(ipMat(is,js)),nOrb(is),'N',
     &              Focka(ipMat(js,is)),nOrb(js),'T',
     %              Fock(ipMat(is,js)),nOrb(is),
     &              nOrb(is),nOrb(js))
      End Do
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Fock=',DDot_(nDens2,Fock,1,Fock,1)
#endif
*
      Call AddGrad(rKappa,Fock,idsym,2.0d0*fact)
      If (.not.newCho) Then
        Call mma_allocate(MT3,nmba,Label='MT3')
        Call DZAXPY(nmba,1.0d0,MT1,1,
     &             MT2,1,MT3,1)
        Call PickMO_MCLR(MT3,rmos,idsym)
*
        If (ispop.ne.0) Then
           Call DZAXPY(nmba,-1.0d0,MT1,1,MT2,1,MT3,1)
           Call PickMO_MCLR(MT3,rmoa,idsym)
        End If
        Call mma_deallocate(MT3)
        Call mma_deallocate(MT2)
        Call mma_deallocate(MT1)
      EndIf
#ifdef _DEBUGPRINT_
      Write (LuWr,*) 'Exit RInt_Generic'
#endif
      Return
      End
