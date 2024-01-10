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
      use Arrays, only: W_CMO_Inv=>CMO_Inv, W_CMO=>CMO, G1t, G2t, FAMO,
     &                  FIMO
*
*                              ~
*     Constructs  F  = <0|[E  ,H]|0> ( + <0|[[E  , Kappa],H]|0> )
*                  pq       pq                 pq
*
      use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
      Implicit Real*8(a-h,o-z)

#include "real.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "standard_iounits.fh"
      Real*8 Fock(nDens2),focka(nDens2),rkappa(nDens2),
     &       Focki(ndens2),Q(ndens2),rMOs(*),rmoa(*)
      Logical Fake_CMO2,DoAct
      Real*8, Allocatable:: MT1(:), MT2(:), MT3(:), QTemp(:),
     &                      Dens2(:),  G2x(:)
      Type (DSBA_Type) CVa(2), DLT(1), DI, DA, Kappa, JI(1), KI, JA, KA,
     &                 FkI, FkA, QVec, CMO, CMO_Inv
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
*                                                                      *
************************************************************************
*                                                                      *
      Select Case (NewCho)
*                                                                      *
************************************************************************
*                                                                      *
      Case (.FALSE.)   ! Cho-MO
*                                                                      *
************************************************************************
*                                                                      *
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
           Call CreQ(Q,MT1,G2t,idsym)
           Call CreQ(QTemp,MT2,G2t,idsym)
#ifdef _DEBUGPRINT_
           Write (LuWr,*) 'Q=',DDot_(nDens2,Q,1,Q,1)
           Write (LuWr,*) 'QTemp=',DDot_(nDens2,QTemp,1,QTemp,1)
#endif
           call daxpy_(ndens2,One,QTemp,1,Q,1)
           Call mma_deallocate(QTemp)
        End If
*                                                                      *
************************************************************************
*                                                                      *
      Case (.TRUE.)   ! Cho-Fock
*                                                                      *
************************************************************************
        Fake_CMO2=.false.
        DoAct=.true.
*
**      Form inactive density
*
        Call Allocate_DT(DI,nOrb,nOrb,nSym)
        DI%A0(:)=Zero

        Do iS=1,nsym
         Do iB=1,nIsh(iS)
          DI%SB(iS)%A2(ib,ib)=Two
         End Do
        End Do
*
**      Form AO 1-index transform inactive density
*
        Call mma_allocate(Dens2,nDens2,Label='Dens2')
        Dens2(:)=Zero
        Call Allocate_DT(DLT(1),nOrb,nOrb,nSym) ! Note SQ format
        DLT(1)%A0(:)=Zero

        If (idSym/=1) Then
           Write (6,*) 'idSym/=1, idSym=',idsym
           Call Abend()
        End If

        Do iS=1,nSym
          If (nOrb(iS).ne.0) Then
            Do jS=1,nSym
              If (iEOr(iS-1,jS-1)+1.eq.idsym.and.nOrb(jS).ne.0) Then
                   Call DGEMM_('N','N',
     &                         nOrb(iS),nOrb(jS),nOrb(jS),
     &                         One,rkappa(ipMat(is,js)),nOrb(iS),
     &                             DI%SB(js)%A2,nOrb(jS),
     &                         Zero,Dens2(ipMat(iS,jS)),nOrb(iS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         One,Dens2(ipMat(iS,jS)),nOrb(iS),
     &                             W_CMO(ipCM(is)),nOrb(iS),
     &                         Zero,DLT(1)%SB(iS)%A2,nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),
     &                         One,DLT(1)%SB(iS)%A2,nOrb(iS),
     &                             W_CMO(ipCM(js)),nOrb(jS),
     &                         Zero,Dens2(ipMat(iS,jS)),nOrb(jS))

                   Call DGEMM_('T','T',nOrb(jS),nOrb(iS),nOrb(iS),
     &                         One,DI%SB(js)%A2,nOrb(iS),
     &                             W_CMO(ipCM(is)),nOrb(iS),
     &                         Zero,DLT(1)%SB(iS)%A2,nOrb(jS))
                   Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nOrb(iS),
     &                         One, DLT(1)%SB(iS)%A2,nOrb(iS),
     &                              W_CMO(ipCM(js)),nOrb(jS),
     &                         Zero,DI%SB(js)%A2,nOrb(jS))
              EndIf
            End Do
          EndIf
        End Do
        ! Release and allocate again in LT format
        Call Deallocate_DT(DLT(1))
        Call Allocate_DT(DLT(1),nOrb,nOrb,nSym,aCase='TRI')

        call Fold_Mat(nSym,nOrb,Dens2,DLT(1)%A0)
        DLT(1)%A0(:) = ReCo * DLT(1)%A0(:)
*
**      Form active density and MO coefficients
*
        If (iMethod.eq.2) Then
          na2=0
          nAct=0
          nG2=0
          Do iSym=1,nSym
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

          Call Allocate_DT(CVa(1),nAsh,nOrb,nSym)
          Call Allocate_DT(CVa(2),nAsh,nOrb,nSym)
          Call Allocate_DT(DA,nAsh,nAsh,nSym)
*
          ioff=0
          ioff4=1
          Do iS=1,nSym
            ioff2 = ioff + nOrb(iS)*nIsh(iS)
            ioff5 = ioff4+ nOrb(iS)*nIsh(iS)
            Do iB=1,nAsh(iS)
              ioff3=ioff2+nOrb(iS)*(iB-1)
              CVa(1)%SB(iS)%A2(iB,:) =
     &           W_CMO(ioff3+1:ioff3+nOrb(iS))
             Do jB=1,nAsh(iS)
              ip2=itri(nA(is)+ib,nA(is)+jb)
              DA%SB(iS)%A2(iB,jB)=G1t(ip2)
             End Do
            End Do
*MGD to check
            Call DGEMM_('T','T',nAsh(iS),nOrb(iS),nOrb(iS),
     &                  1.0d0,rkappa(ioff5),nOrb(iS),
     &                        W_CMO(1+ioff),nOrb(iS),
     &                  0.0d0,CVa(2)%SB(iS)%A2,nAsh(iS))
            ioff=ioff+(nIsh(iS)+nAsh(iS))*nOrb(iS)
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
                        ipGx=ipGx+1
                        G2x(ipGx)=G2t(itri(iij,ikl))
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
        call dcopy_(nATri,[0.0d0],0,rMOs,1)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
        Call RecPrt('DLT',' ',DLT(1)%A0,1,SIZE(DLT(1)%A0))
        Call RecPrt('DI ',' ',DI%A0 ,1,SIZE(DI%A0 ))
        Call RecPrt('DA ',' ',DA%A0 ,1,SIZE(DA%A0 ))
        Call RecPrt('G2x',' ',G2x,1,SIZE(G2x))
        Call RecPrt('rKappa ',' ',rKappa ,1,SIZE(rKappa))
#endif
*
**      Compute the whole thing
*
        Call Allocate_DT(Kappa,nBas,nBas,nSym,Ref=rKappa)
        Call Allocate_DT(JI(1),nBas,nBas,nSym,aCase='TRI')
        JI(1)%A0(:)=Zero
        Call Allocate_DT(KI,nBas,nBas,nSym)
        KI%A0(:)=Zero
        Call Allocate_DT(JA,nBas,nBas,nSym)
        JA%A0(:)=Zero
        Call Allocate_DT(KA,nBas,nBas,nSym)
        KA%A0(:)=Zero
        Call Allocate_DT(FkI,nBas,nBas,nSym,Ref=FockI)
        FkI%A0(:)=Zero
        Call Allocate_DT(FkA,nBas,nBas,nSym,Ref=FockA)
        FkA%A0(:)=Zero
        Call Allocate_DT(QVec,nBas,nAsh,nSym,Ref=Q)
        Call Allocate_DT(CMO,nBas,nAsh,nSym,Ref=W_CMO)
        Call Allocate_DT(CMO_Inv,nBas,nAsh,nSym,Ref=W_CMO_Inv)
        iread=2 ! Asks to read the half-transformed Cho vectors

        Call CHO_LK_MCLR(DLT,DI,DA,G2x,Kappa,JI,KI,JA,KA,FkI,FkA,
     &                   rMOs,QVec,CVa,CMO,CMO_inv,
     &                   nIsh, nAsh,DoAct,Fake_CMO2,
     &                   LuAChoVec,LuIChoVec,iread)

        Call Deallocate_DT(CMO_Inv)
        Call Deallocate_DT(CMO)
        Call Deallocate_DT(QVec)
        Call Deallocate_DT(FkA)
        Call Deallocate_DT(FkI)
        Call Deallocate_DT(KA)
        Call Deallocate_DT(JA)
        Call Deallocate_DT(KI)
        Call Deallocate_DT(JI(1))
        Call Deallocate_DT(Kappa)
        Call GADSum(FockI,nDens2)
        Call GADSum(FockA,nDens2)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      nas=0
      Do iSym = 1, nSym
        Write (6,*) 'iSym=',iSym
*       Call RecPrt('FIMO ',' ', FIMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
*       Call RecPrt('FAMO ',' ', FAMO(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
*       Call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
        nas = nas + nAsh(iSym)
      End Do
      nAtri=nas*(nas+1)/2
      nAtri=nAtri*(nAtri+1)/2
*     Call RecPrt('MO1',' ',rMOs,1,nAtri)
#endif
*
**      Calculate contribution from uncontracted indexes
*
        Do iS=1,nSym
          jS=iEOr(iS-1,iDSym-1)+1
          If (nOrb(iS)*nOrb(jS).ne.0) Then
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(iS),Reco*Fact,
     &                  FIMO(ipCM(iS)),nOrb(is),
     &                  rkappa(ipMat(iS,jS)),nOrb(iS),
     &                  One,FockI(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('N','N',nOrb(iS),nOrb(jS),nOrb(jS),Fact,
     &                  rkappa(ipMat(iS,jS)),nOrb(is),
     &                  FIMO(ipCM(jS)),nOrb(jS),
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

        If (iMethod.eq.2) Then
          Call mma_deallocate(G2x)
          Call Deallocate_DT(CVa(2))
          Call Deallocate_DT(CVa(1))
          Call deallocate_DT(DA)
        EndIf

        Call deallocate_DT(DLT(1))
        Call deallocate_DT(DI)

        Call GADSum(    Q,nDens2)
        Call GADSum( rMOs,nAtri)
*                                                                      *
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
*                                                                      *
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      If (NewCho) Then
      nas=0
      Do iSym = 1, nSym
        Write (6,*) 'iSym=',iSym
        Call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
        nas = nas + nAsh(iSym)
      End Do
      nAtri=nas*(nas+1)/2
      nAtri=nAtri*(nAtri+1)/2
      Call RecPrt('MO1',' ',rMOs,1,nAtri)
      Call abend()
      End If
#endif

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
