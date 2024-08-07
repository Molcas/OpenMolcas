!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2011, Thomas Bondo Pedersen                            *
!               2011, Roland Lindh                                     *
!***********************************************************************

#include "compiler_features.h"
#ifdef _IN_MODULE_

      Subroutine OneEl_IJ(iS,jS,iPrint,Do_PGamma,                       &
     &                    xZeta,xZI,xKappa,xPCoor,                      &
     &                    Kernel,KrnlMm,Label,lOper,nComp,CoorO,        &
     &                    nOrdOp,iChO,                                  &
     &                    iStabO,nStabO,nIC,                            &
     &                    PtChrg,nGrid,iAddPot,SOInt,l_SOInt,           &
     &                    Final,nFinal,Scrtch,nScrtch,                  &
     &                    ScrSph,nScrSph,Kern,nKern)
!
!     Thomas Bondo Pedersen and Roland Lindh, February 2011.
!
!     Purpose: compute symmetry adapted one-electron integrals for
!              shell doublet iS, jS.
!
      use Real_Spherical, only: ipSph, rSph
      use iSD_data, only: iSD
      use Basis_Info, only: DBSC, Shells, MolWgh
      use Center_Info, only: DC
      use Sizes_of_Seward, only:S
      use Gateway_Info, only: FNMC
      use Symmetry_Info, only: ChOper, nIrrep
      use Constants, only: Zero, One
      use rmat, only: RMat_Type_Integrals
      use define_af, only: AngTp
      use property_label, only: PLabel
      Implicit None
      Procedure(int_kernel) :: Kernel
      Procedure(int_mem) :: KrnlMm
#include "Molcas.fh"
      Integer iS, jS, iPrint, nComp, nOrdOp, nStabO, nIC, nGrid,        &
     &        iAddPot, l_SOInt
      Logical Do_PGamma
      Integer nFinal, nScrtch, nScrSph
      Real*8, target :: Final(nFinal)
      Real*8 Scrtch(nScrtch), ScrSph(nScrSph), SOInt(l_SOInt)
      Real*8 xZeta(*),xZI(*),xKappa(*),xPCoor(*)
      Real*8 CoorO(3,nComp), PtChrg(nGrid)
      Integer lOper(nComp), iChO(nComp), iStabO(0:7)


      Real*8 A(3), B(3), RB(3)
      Character(LEN=8) Label
      Character(LEN=LENIN) dbas
      Integer nOp(2), iDCRR(0:7), iDCRT(0:7), iStabM(0:7)
      Integer nKern
      Real*8, Target:: Kern(nKern)
      Integer, external:: MemSO1

      Real*8      Coord(3*MxAtom)
#ifdef _GEN1INT_
      Logical NATEST, DO_TRAN
#endif
      Integer i, iCmp, iBas, iAO, iShell, jCmp, jBas, jAO, jShell, nSO, &
     &        iComp, iSmLbl, iShll, iAng, iPrim, mdci, iCnttp,          &
     &        iCnt, jShll, jAng, jPrim, mdcj, jCnttp, lFinal, ii,       &
     &        Lmbdr, nDCRR, nStabM, lDCRR, l_Coord, nAtoms, iSOBlk, iiC,&
     &        mSO, iIrrep, lA0, lA1, lB0, lB1, MemAux, MemBux, MemCux,  &
     &        MemKer, MemKrn, lScrtch, lScrSph, iuv, LambdT, kk, ipFnl, &
     &        nij, nijab, iab, ipX, ipY, ipZ, jj, iAtom, nDCRT,         &
     &        nOrder, NrOpr, jCnt
      Real*8 Fact, xFactor, xMass
      Integer :: iTwoj(0:7)=[1,2,4,8,16,32,64,128]
!
!     Statement functions
      Integer ixyz, nElem
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Label(1:3).eq.'MAG') Then
      iCmp   = iSD( 2,iS)
      iBas   = iSD( 3,iS)
      iAO    = iSD( 7,iS)
      iShell = iSD(11,iS)
      jCmp   = iSD( 2,jS)
      jBas   = iSD( 3,jS)
      jAO    = iSD( 7,jS)
      jShell = iSD(11,jS)
      nSO=0
      B(:)=Zero
      Do iComp = 1, nComp
         iSmLbl=lOper(iComp)
         nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
      End Do
      If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
      If (nSO.lt.1) Return
      If (l_SOInt.lt.nSO*iBas*jBas) Then
         Call WarningMessage(2,                                         &
     &                       'OneEl_IJ: insufficient SOInt dimension!')
         Call Abend()
      End If
      Call dCopy_(nSO*iBas*jBas,[Zero],0,SOInt,1)
      iShll  = iSD( 0,iS)
      iAng   = iSD( 1,iS)
      iPrim  = iSD( 5,iS)
      mdci   = iSD(10,iS)
      iCnttp = iSD(13,iS)
      iCnt   = iSD(14,iS)
      A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
      dbas= dc(mdci)%LblCnt(1:LENIN)
      Call UpCase(dbas)
      jShll  = iSD( 0,jS)
      jAng   = iSD( 1,jS)
      jPrim  = iSD( 5,jS)
      mdcj   = iSD(10,jS)
      jCnttp = iSD(13,jS)
      jCnt   = iSD(14,jS)
      B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
       if (iPrint.ge.19) Then
        write(6,*) "interacted Ato.Fun "
          Write (6,'(A,A,A,A,A)')' ***** (',AngTp(iAng),',',AngTp(jAng),&
     &                           ') *****'
       endif
      lFinal = nIC*S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
      If (lFinal.gt.nFinal) Then
         Call WarningMessage(2,'lFinal.gt.nFinal')
         Call Abend()
      End If
      Call dCopy_(lFinal,[Zero],0,Final,1)
      Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,                     &
     &               dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
      Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,                         &
     &           dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)
      Call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
      If (iPrint.ge.19) Then
         Write (6,*)
         Write (6,*) ' g      =',nIrrep
         Write (6,*) ' u      =',dc(mdci)%nStab
         Write (6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii = 0,    &
     &                    dc(mdci)%nStab-1)
         Write (6,*) ' v      =',dc(mdcj)%nStab
         Write (6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii = 0,    &
     &                    dc(mdcj)%nStab-1)
         Write (6,*) ' LambdaR=**',LmbdR
         Write (6,*) ' r      =',nDCRR
         Write (6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii = 0, nDCRR-1)
         Write (6,*) ' m      =',nStabM
         Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii = 0, nStabM-1)
      End If
       nOp(1) = NrOpr(0)
      If (nDCRR.ge.1) Then
         Do lDCRR = 0, nDCRR-1
            Call OA(iDCRR(lDCRR),B,RB)
            nOp(2) = NrOpr(iDCRR(lDCRR))
               If (iPrint.ge.49) Then
               Write (6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),         &
     &                                        (RB(i),i=1,3)
               Endif

        Call Get_nAtoms_All(nAtoms)
        l_Coord=3*nAtoms
        Call Get_dArray('Bfn Coordinates',Coord,l_Coord)
        if (nAtoms.eq.2) then
         nAtoms=nAtoms+1
        endif
#ifdef _GEN1INT_
        NATEST=(nAtoms.eq.2)
        If (label(4:5).eq.'PX') Then
            Do_Tran = .TRUE.
        Else
            Do_Tran = .FALSE.
        End If
        read(Label(6:),'(I3)') iatom
        call test_f90mod_sgto_mag(iShell,jShell,iCmp,jCmp,              &
       &                          iPrim,jPrim,iAng,jAng,                &
       &                          iPrim,jPrim,mdci,mdcj,                &
       &                          Shells(iShll)%Exp,Shells(jShll)%Exp,  &
       &                          Shells(iShll)%Cff_p(1,1,2),           &
       &                          Shells(jShll)%Cff_p(1,1,2),           &
       &                          nAtoms,                               &
       &                          Coord,nComp,Final,.TRUE.,             &
       &                          iatom,Do_Tran)
#else
         Call WarningMessage(2,'OneEl_IJ: NO Gen1int interface '//      &
     &                       'available!')
         Call Abend()
#endif

            iSOBlk = 1
            iIC = 1
            Do iComp = 1, nComp
               iSmLbl=lOper(iComp)
               mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
               If (mSO.eq.0) Then
                  Do iIrrep = 0, nIrrep-1
                     If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0)         &
     &                  iIC = iIC + 1
                  End Do
               Else
                 !write(6,*) "Symmetry adapt component"
                  Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,               &
     &                        iShell,jShell,iShll,jShll,                &
     &                        iAO,jAO,Final,                            &
     &                        iBas,jBas,nIC,iIC,SOInt(iSOBlk),mSO,nOp)
                  iSOBlk = iSOBlk + mSO*iBas*jBas
               End If
            End Do
         End Do
      End If

         else  !  MAG Integrals

!                                                                      *
!***********************************************************************
!                                                                      *
!     Check memory for SO integrals that will be generated by
!     this batch of AO integrals. Init SOInt.
!
      iCmp   = iSD( 2,iS)
      iBas   = iSD( 3,iS)
      iAO    = iSD( 7,iS)
      iShell = iSD(11,iS)
      jCmp   = iSD( 2,jS)
      jBas   = iSD( 3,jS)
      jAO    = iSD( 7,jS)
      jShell = iSD(11,jS)
      nSO=0
      Do iComp = 1, nComp
         iSmLbl=lOper(iComp)
         nSO=nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
      End Do
      If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
      If (nSO.lt.1) Return
      If (l_SOInt.lt.nSO*iBas*jBas) Then
         Call WarningMessage(2,                                         &
     &                       'OneEl_IJ: insufficient SOInt dimension!')
         Call Abend()
      End If
      Call dCopy_(nSO*iBas*jBas,[Zero],0,SOInt,1)
!
!---- Shell info
!
      iShll  = iSD( 0,iS)
      iAng   = iSD( 1,iS)
      iPrim  = iSD( 5,iS)
      mdci   = iSD(10,iS)
      iCnttp = iSD(13,iS)
      iCnt   = iSD(14,iS)
      A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
      dbas= dc(mdci)%LblCnt(1:LENIN)
      Call UpCase(dbas)

      jShll  = iSD( 0,jS)
      jAng   = iSD( 1,jS)
      jPrim  = iSD( 5,jS)
      mdcj   = iSD(10,jS)
      jCnttp = iSD(13,jS)
      jCnt   = iSD(14,jS)
      B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
!
!---- Identify if shell doublet should be computed with special
!     R-Matrix code.
!
      If (iCnttp.eq.jCnttp .and. mdcj.eq.mdci .and. dbas.eq.'DBAS' )    &
     &   Then
         RMat_type_integrals=.True.
         If (Do_PGamma) Then
            Call PGamma
            Do_PGamma=.False.
         End If
      Else
         RMat_type_integrals=.False.
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
      If (iPrint.ge.19) Then
      write(6,*) "interacted Ato.Fun "
         Write (6,'(A,A,A,A,A)')' ***** (',AngTp(iAng),',',AngTp(jAng), &
     &                          ') *****'
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Call kernel routine to get memory requirement.
!
      Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
!
!     Special additional allocation for PCC integrals
!
      If (PLabel.ne.' ') Then
         la0=iAng
         lb0=jAng
         MemAux= 1 + 3*nElem(la0)*nElem(lb0+1)*nIC
         la1=la0
         lb1=lb0+1
         MemBux= 1 + 3*nElem(la1+1)*nElem(lb1)*nIC
         If (la1.ne.0) MemBux=MemBux+3*nElem(la1-1)*nElem(lb1)*nIC
         If (lb0.ne.0) Then
            lb1=lb0-1
            MemAux=MemAux+3*nElem(la0)*nElem(lb0-1)*nIC
            MemCux=1+3*nElem(la1+1)*nElem(lb1)*nIC
            If (la1.ne.0) MemCux=MemCux+3*nElem(la1-1)*nElem(lb1)*nIC
         Else
            MemCux=0
         End If
         MemAux = MemAux + Max(MemBux,MemCux)
         MemKer = MemKer + MemAux
      End If
!
      MemKrn=MemKer*iPrim*jPrim
      If (MemKrn.gt.nKern) Then
         Call WarningMessage(2,'MemKrn.gt.nKern')
         Write (6,*)'nOrdOp,iAng,jAng=',nOrdOp,iAng,jAng
         Write (6,*) 'MemKrn=',MemKrn
         Write (6,*) 'nKern=',nKern
         Call Abend()
      End If
      Call dCopy_(MemKrn,[Zero],0,Kern,1)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Allocate memory for the final integrals all in the
!     primitive basis.
!
      lFinal = nIC*iPrim*jPrim*nElem(iAng)*nElem(jAng)
      If (lFinal.gt.nFinal) Then
         Call WarningMessage(2,'lFinal.gt.nFinal')
         Call Abend()
      End If
      Call dCopy_(lFinal,[Zero],0,Final,1)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Scratch area for contraction step
!
      lScrtch =  Max(iPrim,jPrim) * Max(iBas,jBas) * nIC*nElem(iAng)*   &
     &           nElem(jAng)
      If (lScrtch.gt.nScrtch) Then
         Call WarningMessage(2,'lScrtch.gt.nScrtch')
         Call Abend()
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Scratch area for the transformation to spherical gaussians
!
      lScrSph=nIC*iBas*jBas*nElem(iAng)*nElem(jAng)
      If (lScrSph.gt.nScrSph) Then
         Call WarningMessage(2,'lScrSph.gt.nScrSph')
         Call Abend()
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     At this point we can compute Zeta.
!     This is now computed in the ij or ji order.
!
      Call ZXia(xZeta,xZI,iPrim,jPrim,Shells(iShll)%Exp,                &
     &                                Shells(jShll)%Exp)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Find the DCR for A and B
!
      Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,                     &
     &               dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
!
!     Find the stabilizer for A and B
!
      Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,                         &
     &           dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)
!
      Call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
!
#ifdef _DEBUGPRINT_
      If (iPrint.ge.19) Then
         Write (6,*)
         Write (6,*) ' g      =',nIrrep
         Write (6,*) ' u      =',dc(mdci)%nStab
         Write (6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii = 0,    &
     &                    dc(mdci)%nStab-1)
         Write (6,*) ' v      =',dc(mdcj)%nStab
         Write (6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii = 0,    &
     &                    dc(mdcj)%nStab-1)
         Write (6,*) ' LambdaR=**',Lmbd
         Write (6,*) ' r      =',nDCRR
         Write (6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii = 0, nDCRR-1)
         Write (6,*) ' m      =',nStabM
         Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii = 0, nStabM-1)
      End If
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute normalization factor
!
      iuv = dc(mdci)%nStab*dc(mdcj)%nStab
      If (MolWgh.eq.1) Then
         Fact = DBLE(nStabO) / DBLE(LambdT)
      Else If (MolWgh.eq.0) Then
         Fact = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LambdT)
      Else
         Fact = Sqrt(DBLE(iuv))*DBLE(nStabO)/DBLE(nirrep*LambdT)
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Loops over symmetry operations acting on the basis.
!
       nOp(1) = NrOpr(0)
      If (nDCRR.ge.1) Then
         Do lDCRR = 0, nDCRR-1
            Call OA(iDCRR(lDCRR),B,RB)
            nOp(2) = NrOpr(iDCRR(lDCRR))
             If (iPrint.ge.49) Then
               Write (6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),         &
     &                                        (RB(i),i=1,3)
            End If
!
!           Compute kappa and P.
!
            Call Setup1(Shells(iShll)%Exp,iPrim,                        &
     &                  Shells(jShll)%Exp,jPrim,A,RB,xKappa,xPCoor,xZI)
!
!           Compute primitive integrals. Result is ordered ij,ab.
!
            Call Kernel(Shells(iShll)%Exp,iPrim,                        &
     &                  Shells(jShll)%Exp,jPrim,                        &
     &                  xZeta,xZI,xKappa,xPCoor,                        &
     &                  Final,iPrim*jPrim,nIC,nComp,                    &
     &                  iAng,jAng,A,RB,nOrder,Kern,                     &
     &                  MemKer,CoorO,nOrdOp,lOper,iChO,iStabM,          &
     &                  nStabM,PtChrg,nGrid,iAddPot)
                If (iPrint.ge.49) Then
                     Call RecPrt(' Primitive Integrals',' ',Final,iPrim*&
     &                           jPrim,nElem(iAng)*nElem(jAng)*nIC)
              End If
!
!           Transform from primitive to contracted basis functions.
!
             If (iPrint.ge.99) Then
                  Call RecPrt(' Left side contraction',' ',             &
     &                        Shells(iShll)%pCff,iPrim,iBas)
                  Call RecPrt(' Right side contraction',' ',            &
     &                        Shells(jShll)%pCff,jPrim,jBas)
!
              End If
!
!           Transform i,jabx to jabx,I
            kk=nElem(iAng)*nElem(jAng)
            Call DGEMM_('T','N',jPrim*kk*nIC,iBas,iPrim,                &
     &                  1.0d0,Final,iPrim,                              &
     &                        Shells(iShll)%pCff,iPrim,                 &
     &                  0.0d0,Scrtch,jPrim*kk*nIC)
!           Transform j,abxI to abxI,J
            Call DGEMM_('T','N',kk*nIC*iBas,jBas,jPrim,                 &
     &                  1.0d0,Scrtch,jPrim,                             &
     &                        Shells(jShll)%pCff,jPrim,                 &
     &                  0.0d0,ScrSph,kk*nIC*iBas)
!
            If (iPrint.ge.99) Then
               Call RecPrt(' Contracted integrals in cartesians',' ',   &
     &                     ScrSph,kk*nIC,iBas*jBas)
            End If
!
!           Transform to spherical gaussians if needed.
!

            If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
!              Result comes back as xIJAB or xIJAb
               Call CarSph(ScrSph,kk,iBas*jBas*nIC,                     &
     &                     Final,lScrSph,                               &
     &                     RSph(ipSph(iAng)),                           &
     &                     iAng,Shells(iShll)%Transf,                   &
     &                          Shells(iShll)%Prjct,                    &
     &                     RSph(ipSph(jAng)),                           &
     &                     jAng,Shells(jShll)%Transf,                   &
     &                          Shells(jShll)%Prjct,                    &
     &                     Scrtch,iCmp*jCmp)
               Call DGeTmO(Scrtch,nIC,nIC,iBas*jBas*iCmp*jCmp,Final,    &
     &                     iBas*jBas*iCmp*jCmp)
            Else
!              Transpose abx,IJ back to IJ,abx
              Call DGeTmO(ScrSph,kk*nIC,kk*nIC,iBas*jBas,Final,         &
     &                    iBas*jBas)
            End If
            If (iPrint.ge.99) Then
             Call RecPrt(' Contracted integrals in Sphericals',' ',     &
     &                   Final,iBas*jBas,iCmp*jCmp*nIC)
            End If

!---------- Tweak here for special cases
!
            ipFnl=1
            If (Label.eq.'P_matrix') Then
               nij=iBas*jBas
               nijab=nij*iCmp*jCmp
               Do iab = 1, iCmp*jCmp
                  ipx = ipFnl+(iab-1)*nij
                  ipy = ipx + nijab
                  ipz = ipy + nijab
                  call dcopy_(nij,xPCoor(1     ),1,Final(ipx),1)
                  call dcopy_(nij,xPCoor(1+nij ),1,Final(ipy),1)
                  call dcopy_(nij,xPCoor(1+2*nij),1,Final(ipz),1)
               End Do
            Else If (Label.eq.'FMMCnX') Then
               Do jj = 1, iBas*jBas*iCmp*jCmp*nIC
                  Final(ipFnl+jj-1) = (A(1) + RB(1))/2.0d0
               End do
            Else If (Label.eq.'FMMCnY') Then
               Do jj = 1, iBas*jBas*iCmp*jCmp*nIC
                  Final(ipFnl+jj-1) = (A(2) + RB(2))/2.0d0
               End do
            Else If (Label.eq.'FMMCnZ') Then
               Do jj = 1, iBas*jBas*iCmp*jCmp*nIC
                  Final(ipFnl+jj-1) = (A(3) + RB(3))/2.0d0
               End do
            Else If (Label.eq.'Kinetic') Then
!
!              multiply with 1/m, where m is the mass of an electron
!              or muon.
!
               xfactor=One/dbsc(iCnttp)%fMass
!
!              Add the Finite Nuclear Mass Correction if activated
!
               If (FNMC .AND. (                                         &
     &             A(1).eq.RB(1) .AND.                                  &
     &             A(2).eq.RB(2) .AND.                                  &
     &             A(3).eq.RB(3)) .AND.                                 &
     &             dbsc(iCnttp)%Charge.ne.Zero) Then
                     iAtom=dbsc(iCnttp)%AtmNr
!                    Get the atom mass in au (me=1)
                     xMass=dbsc(iCnttp)%CntMass
!                    Substract the electron mass to get the nuclear
!                    mass.
                     xMass=xMass-DBLE(iAtom)
!                    Write (*,*) 'xMass=',xMass
                     xfactor=xfactor+One/xMass
               End If
!              Write (*,*) 'xfactor=',xfactor
               Call DScal_(iBas*jBas*iCmp*jCmp,xfactor,Final,1)
            End If
!
!           At this point accumulate the batch of integrals onto the
!           final symmetry adapted integrals.
!
             If (iPrint.ge.99) Then
               Call RecPrt (' Accumulated SO integrals, so far...',' ', &
     &                      SOInt,iBas*jBas,nSO)
             End If

!---------- Symmetry adapt component by component
!
            iSOBlk = 1
            iIC = 1
            Do iComp = 1, nComp
               iSmLbl=lOper(iComp)
               mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
               If (mSO.eq.0) Then
                  Do iIrrep = 0, nIrrep-1
                     If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0)         &
     &                  iIC = iIC + 1
                  End Do
               Else
                  Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell, &
     &                        iShll,jShll,iAO,jAO,Final,iBas,jBas,nIC,  &
     &                        iIC,SOInt(iSOBlk),mSO,nOp)
                  iSOBlk = iSOBlk + mSO*iBas*jBas
               End If
            End Do
         End Do
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!     Multiply with factors due to projection operators
!
      If (Fact.ne.One) Call dScal_(nSO*iBas*jBas,Fact,SOInt,1)
      If (iPrint.ge.99) Then
         Write (6,*) ' Scaling SO''s', Fact
         Call RecPrt(' Accumulated SO integrals',' ',SOInt,iBas*jBas,   &
     &               nSO)
      End If
!                                                                      *
!***********************************************************************
!                                                                      *

      End If
!
      Return
      End Subroutine OneEl_IJ

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
      dummy_empty_procedure(OneEl_IJ)

#endif
