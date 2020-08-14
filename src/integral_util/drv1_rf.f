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
* Copyright (C) 1990-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv1_RF(FactOp,nOpr,FD,nFD,CCoor,lOper,Cavxyz,lMax)
************************************************************************
*                                                                      *
* Object: to compute the local multipole moment, desymmetrize the 1st  *
*         order density matrix and accumulate contributions to the     *
*         global multipole expansion.                                  *
*                                                                      *
* Called from: RctFld                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              ZXia                                                    *
*              SetUp1                                                  *
*              MltInt                                                  *
*              DGeMV    (ESSL)                                         *
*              RecPrt                                                  *
*              DCopy    (ESSL)                                         *
*              DGEMM_   (ESSL)                                         *
*              CarSph                                                  *
*              DGeTMO   (ESSL)                                         *
*              DaXpY    (ESSL)                                         *
*              SOGthr                                                  *
*              DesymD                                                  *
*              DScal    (ESSL)                                         *
*              TriPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Modified for Hermite-Gauss quadrature November '90       *
*             Modified for Rys quadrature November '90                 *
*             Modified for multipole moments November '90              *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for general kernel routines January  91         *
*             Modified for nonsymmetrical operators February  91       *
*             Modified for gradients October  91                       *
*             Modified for reaction field calculations July  92        *
*             Modified loop structure  99                              *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), C(3), FD(nFD), FactOp(nOpr), CCoor(3,nOpr),
     &       RB(3), TRB(3), TA(3), Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Character ChOper(0:7)*3
      Integer   lOper(nOpr), iStabO(0:7),
     &          iDCRR(0:7), iDCRT(0:7), iStabM(0:7), nOp(3)
      Logical AeqB
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
      Real*8, Allocatable:: Zeta(:), ZI(:), Kappa(:), PCoor(:,:)
      Real*8, Allocatable:: Kern(:), Fnl(:), Scr1(:), Scr2(:),
     &                      DAO(:), DSOpr(:), DSO(:)
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
      Call qEnter('Drv1_RF')
*
      iIrrep = 0
*
*     Auxiliary memory allocation.
*
      Call mma_allocate(Zeta,m2Max,Label='Zeta')
      Call mma_allocate(ZI,m2Max,Label='ZI')
      Call mma_allocate(Kappa,m2Max,Label='Kappa')
      Call mma_allocate(PCoor,m2Max,3,Label='PCoor')
*                                                                      *
************************************************************************
*                                                                      *
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         If (AuxShell(iShll)) Go To 100
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            If (nSO.eq.0) Go To 131
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &        ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*
*           Call kernel routine to get memory requirement.
*
            nOrdOp = lMax
            Call RFMem(nOrder,MemKer,iAng,jAng,nOrdOp)
*           Write (*,*)nOrder,MemKer,iAng,jAng,nOrdOp
            MemKrn=MemKer*m2Max
            Call mma_allocate(Kern,MemKrn,Label='Kern')
*
*           Allocate memory for the final integrals, all in the
*           primitive basis.
*
            nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
            lFinal = MaxPrm(iAng) * MaxPrm(jAng)
     &             * nElem(iAng)*nElem(jAng)
     &             * nComp
            Call mma_allocate(Fnl,lFinal,Label='Fnl')
*
*           Scratch area for contraction step
*
            nScr1 =  MaxPrm(iAng)*MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr1,nScr1,Label='Scr1')
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=MaxPrm(iAng)*MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr2,nScr2,Label='Scr2')
*
            nDAO =iPrim*jPrim*nElem(iAng)*nElem(jAng)
            Call mma_allocate(DAO,nDAO,Label='DAO')
*
*           At this point we can compute Zeta.
*
            Call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,
     &                                    Shells(jShll)%Exp)
*
            AeqB = iS.eq.jS
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,jStab(0,mdci),
     &               nStab(mdci),jStab(0,mdcj),
     &               nStab(mdcj),iDCRR,nDCRR)
            If (iPrint.ge.49) Write (6,'(10A)')
     &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
*-----------Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                 iStabM,nStabM)
*
*           Allocate memory for the elements of the Fock or 1st order
*           denisty matrix which are associated with the current shell
*           pair.
*
            Call mma_allocate(DSOpr,nSO*iPrim*jPrim,Label='DSOpr')
            Call mma_allocate(DSO,nSO*iPrim*jPrim,Label='DSO')
*
*           Gather the elements from 1st order density / Fock matrix.
*
            Call SOGthr(DSO,iBas,jBas,nSO,FD,
     &                  n2Tri(iSmLbl),iSmLbl,
     &                  iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
*
*           Project the Fock/1st order density matrix in AO
*           basis on to the primitive basis.
*
            If (iPrint.ge.99) Then
               Call RecPrt(' Left side contraction',' ',
     &                     Shells(iShll)%pCff,iPrim,iBas)
               Call RecPrt(' Right side contraction',' ',
     &                     Shells(jShll)%pCff,jPrim,jBas)
            End If
*
*           Transform IJ,AB to J,ABi
            Call DGEMM_('T','T',
     &                  jBas*nSO,iPrim,iBas,
     &                  1.0d0,DSO,iBas,
     &                        Shells(iShll)%pCff,iPrim,
     &                  0.0d0,DSOpr,jBas*nSO)
*           Transform J,ABi to AB,ij
            Call DGEMM_('T','T',
     &                  nSO*iPrim,jPrim,jBas,
     &                  1.0d0,DSOpr,jBas,
     &                        Shells(jShll)%pCff,jPrim,
     &                  0.0d0,DSO,nSO*iPrim)
*           Transpose to ij,AB
            Call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOpr,
     &                  iPrim*jPrim)
            Call mma_deallocate(DSO)
*
            If (iPrint.ge.99) Call
     &         RecPrt(' Decontracted 1st order density/Fock matrix',
     &                ' ',DSOpr,iPrim*jPrim,nSO)
*
*           Loops over symmetry operations.
*
            Do lDCRR = 0, nDCRR-1
               RB(1)  = DBLE(iPhase(1,iDCRR(lDCRR)))*B(1)
               RB(2)  = DBLE(iPhase(2,iDCRR(lDCRR)))*B(2)
               RB(3)  = DBLE(iPhase(3,iDCRR(lDCRR)))*B(3)
*
*--------------Loop over operators
*
               Do 5000 iOpr = 1, nOpr
                  If (FactOp(iOpr).eq.Zero) Go To 5000
                  call dcopy_(3,Ccoor(1,iOpr),1,C,1)
*
*-----------------Generate stabilizer of the operator.
*
                  Call SOS(iStabO,nStabO,lOper(iOpr))
*
*-----------------Find the DCR for M and S
*
                  Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &                     iStabO,nStabO,iDCRT,nDCRT)
                  If (iPrint.ge.49) Then
                     Write (6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),
     &                     i=0,nStabM-1),')'
                     Write (6,'(10A)') ' {O}=(',(ChOper(iStabO(i)),
     &                     i=0,nStabO-1),')'
                     Write (6,'(10A)') ' {T}=(',(ChOper(iDCRT(i)),
     &                     i=0,nDCRT-1),')'
                  End If
*
*-----------------Compute normalization factor due the DCR symmetrization
*                 of the two basis functions and the operator.
*
                  iuv = nStab(mdci)*nStab(mdcj)
                  FactNd = DBLE(iuv*nStabO) / DBLE(nIrrep**2*LmbdT)
                  If (MolWgh.eq.1) Then
                     FactNd = FactNd * DBLE(nIrrep)**2 / DBLE(iuv)
                  Else If (MolWgh.eq.2) Then
                     FactNd = Sqrt(DBLE(iuv))*DBLE(nStabO) /
     &                        DBLE(nIrrep*LmbdT)
                  End If
                  FactNd = FactNd * FactOp(iOpr)
*
                  Do lDCRT = 0, nDCRT-1
                     nOp(1) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
                     nOp(2) = NrOpr(iEor(iDCRT(lDCRT),
     &                             iDCRR(lDCRR)),iOper,nIrrep)
                     nOp(3) = NrOpr(0,iOper,nIrrep)

                     TA(1) =  DBLE(iPhase(1,iDCRT(lDCRT)))*A(1)
                     TA(2) =  DBLE(iPhase(2,iDCRT(lDCRT)))*A(2)
                     TA(3) =  DBLE(iPhase(3,iDCRT(lDCRT)))*A(3)
                     TRB(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*RB(1)
                     TRB(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*RB(2)
                     TRB(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*RB(3)
                     If (iPrint.ge.49) Then
                        Write (6,'(A,/,3(3F6.2,2X))')
     &                  ' *** Centers A, B, C ***',
     &                  ( TA(i),i=1,3),
     &                  (TRB(i),i=1,3),
     &                  (C(i),i=1,3)
                        Write (6,*) ' nOp=',nOp
                     End If
*
*--------------------Desymmetrize the matrix with which we will
*                    contracte the trace.
*
                     Call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                           iShell,jShell,iShll,jShll,
     &                           DAO,iPrim,jPrim,
     &                           DSOpr,nSO,nOp,FactNd)
*
*--------------------Project the spherical harmonic space onto the
*                    cartesian space.
*
                     kk = nElem(iAng)*nElem(jAng)
                     If (Shells(iShll)%Transf.or.
     &                   Shells(jShll)%Transf) Then
*
*-----------------------ij,AB --> AB,ij
                        Call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,
     &                              iCmp*jCmp,Scr1,iCmp*jCmp)
*-----------------------AB,ij --> ij,ab
                        Call SphCar(Scr1,iCmp*jCmp,iPrim*jPrim,
     &                              Scr2,nScr2,
     &                              RSph(ipSph(iAng)),
     &                              iAng,Shells(iShll)%Transf,
     &                                   Prjct(iShll),
     &                              RSph(ipSph(jAng)),
     &                              jAng,Shells(jShll)%Transf,
     &                                   Prjct(jShll),
     &                              DAO,kk)
                     End If
                     If (iPrint.ge.99) Call RecPrt(
     &                        ' Decontracted FD in the cartesian space',
     &                        ' ',DAO,iPrim*jPrim,kk)
*
*--------------------Compute kappa and P.
*
                     Call Setup1(Shells(iShll)%Exp,iPrim,
     &                           Shells(jShll)%Exp,jPrim,
     &                           TA,TRB,Kappa,PCoor,ZI)
*
*
*--------------------Compute primitive multipole moments.
*
                     Call RFInt(Shells(iShll)%Exp,iPrim,
     &                          Shells(jShll)%Exp,jPrim,
     &                          Zeta,ZI,Kappa,Pcoor,
     &                          Fnl,iPrim*jPrim,nComp,
     &                          iAng,jAng,TA,TRB,nOrder,Kern,
     &                          MemKer,C,nOrdOp)
                     If (iPrint.ge.49) Call RecPrt(' Final Integrals',
     &                                 ' ',Fnl,nDAO,nComp)
*
*--------------------Trace with 1st order density matrix and accumulate
*                    to the multipole expansion around center Q.
*
                     If (iPrint.ge.49) Call RecPrt(
     &                        ' Decontracted FD in the cartesian space',
     &                        ' ',DAO,nDAO,1)
                     If (iPrint.ge.49) Call RecPrt('Cavxyz',' ',
     &                                             Cavxyz,1,nComp)
                     Call dGeMV_('T',nDAO,nComp,
     &                         -One,Fnl,nDAO,
     &                              DAO,1,
     &                          One,Cavxyz,1)
                     If (iPrint.ge.49) Call RecPrt('Cavxyz',' ',
     &                                             Cavxyz,1,nComp)
*
                  End Do
 5000          Continue
            End Do
*
            Call mma_deallocate(DSOpr)
            Call mma_deallocate(DAO)
            Call mma_deallocate(Scr2)
            Call mma_deallocate(Scr1)
            Call mma_deallocate(Fnl)
            Call mma_deallocate(Kern)
 131        Continue
         End Do
      End Do
 100  Continue
*
      Call mma_deallocate(PCoor)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(ZI)
      Call mma_deallocate(Zeta)
*
      Call qExit('Drv1_RF')
      Return
      End
