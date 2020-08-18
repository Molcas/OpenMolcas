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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine OneEl_g(Kernel,KrnlMm,Grad,nGrad,DiffOp,CCoor,FD,
     &                        nFD,lOper,nComp,nOrdOp,Label)
************************************************************************
*                                                                      *
* Object: to compute gradients of the one electron integrals.          *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
* Called from: Drvh1                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              ZXia                                                    *
*              SetUp1                                                  *
*              Kernel                                                  *
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
*             Modified for general kernel routines January '91         *
*             Modified for nonsymmetrical operators February '91       *
*             Modified for gradients October '91                       *
************************************************************************
      use Real_Spherical
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "nac.fh"
#include "setup.fh"
CNIKO      Real*8 A(3), B(3), Ccoor(3,nComp), FD(nFD),
      Real*8 A(3), B(3), Ccoor(*), FD(nFD),
     &       RB(3), Grad(nGrad)
      Character ChOper(0:7)*3, Label*80
      Integer iDCRR(0:7), iDCRT(0:7), iStabM(0:7),
     &          IndGrd(3,2), nOp(2), iStabO(0:7), lOper(nComp)
      Logical AeqB, EQ, DiffOp, IfGrad(3,3)
      Logical FreeiSD
      Real*8, Allocatable:: Zeta(:), ZI(:), Kappa(:), PCoor(:,:)
      Real*8, Allocatable:: Krnl(:), Final(:), Scr1(:), Scr2(:)
      Real*8, Allocatable:: DAO(:), DSOpr(:), DSO(:)
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
      call dcopy_(nGrad,[Zero],0,Grad,1)
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
      If (.Not.Allocated(iSD)) Then
         Call Set_Basis_Mode('Valence')
         Call Setup_iSD()
         FreeiSD=.True.
      Else
         FreeiSD=.False.
      End If
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      nTasks = nSkal*(nSkal+1)/2
      iS = 0
      jS = 0
      Do ijS = 1, nTasks
         jS = jS + 1
         If (jS.gt.iS) Then
            iS = jS
            jS = 1
         End If
*
C     Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iCff   = iSD( 4,iS)
         iPrim  = iSD( 5,iS)
         iExp   = iSD( 6,iS)
         iAO    = iSD( 7,iS)
         ixyz   = iSD( 8,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         call dcopy_(3,Work(ixyz),1,A,1)

C        Do jS = 1, iS
            jShll  = iSD( 0,jS)
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jCff   = iSD( 4,jS)
            jPrim  = iSD( 5,jS)
            jExp   = iSD( 6,jS)
            jAO    = iSD( 7,jS)
            jxyz   = iSD( 8,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            call dcopy_(3,Work(jxyz),1,B,1)
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            If (nSO.eq.0) Go To 131
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,jStab(0,mdci),nStab(mdci),
     &                                  jStab(0,mdcj),nStab(mdcj),
     &                                                iDCRR,nDCRR)
*
*           For the CSF contribution we cannot skip the A,A case
*
            If (.Not.isCSF) Then
               If (.Not.DiffOp .and. nDCRR.eq.1 .and. EQ(A,B)) Go To 131
            End If
            If (iPrint.ge.49) Write (6,'(10A)')
     &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &         ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*
*           Call kernel routine to get memory requirement.
*
            Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
            MemKrn=MemKer*m2Max
            Call mma_allocate(Krnl,MemKrn,Label='Krnl')
*
*           Allocate memory for the final integrals, all in the
*           primitive basis.
*
            lFinal = 6 * MaxPrm(iAng) * MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng) * nComp
            Call mma_allocate(Final,lFinal,Label='Final')
*
*           Scratch area for contraction step
*
            nScr1 =  MaxPrm(iAng)*MaxPrm(jAng) * nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr1,nScr1,Label='Scr1')
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=MaxPrm(iAng)*MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr2,nScr2,Label='Scr2')
*
            Call mma_allocate(DAO,iPrim*jPrim*nElem(iAng)*nElem(jAng),
     &                        Label='DAO')
*
*           At this point we can compute Zeta.
*
            Call ZXia(Zeta,ZI,iPrim,jPrim,Work(iExp),Work(jExp))
*
            Do iCar = 0, 2
               IndGrd(iCar+1,1) = iSD(iCar+16,iS)
               IfGrad(iCar+1,1) = iSD(iCar+16,iS).ne.0
            End Do
*
            AeqB = iS.eq.jS
*
            Do iCar = 0, 2
               IndGrd(iCar+1,2) = iSD(iCar+16,jS)
               IfGrad(iCar+1,2) = iSD(iCar+16,jS).ne.0
            End Do
*
*-----------Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                             iStabM,nStabM)
*
*           Allocate memory for the elements of the Fock or 1st order
*           density matrix which are associated with the current shell
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
     &                     Work(iCff),iPrim,iBas)
               Call RecPrt(' Right side contraction',' ',
     &                     Work(jCff),jPrim,jBas)
            End If
*
*           Transform IJ,AB to J,ABi
            Call DGEMM_('T','T',
     &                  jBas*nSO,iPrim,iBas,
     &                  1.0d0,DSO,iBas,
     &                  Work(iCff),iPrim,
     &                  0.0d0,DSOpr,jBas*nSO)
*           Transform J,ABi to AB,ij
            Call DGEMM_('T','T',
     &                  nSO*iPrim,jPrim,jBas,
     &                  1.0d0,DSOpr,jBas,
     &                  Work(jCff),jPrim,
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
            nOp(1) = NrOpr(0,iOper,nIrrep)
c VV: gcc bug: one has to use this if!
          if(nDCRR.ge.1) then
            Do 140 lDCRR = 0, nDCRR-1
               RB(1)  = DBLE(iPhase(1,iDCRR(lDCRR)))*B(1)
               RB(2)  = DBLE(iPhase(2,iDCRR(lDCRR)))*B(2)
               RB(3)  = DBLE(iPhase(3,iDCRR(lDCRR)))*B(3)
               nOp(2) = NrOpr(iDCRR(lDCRR),iOper,nIrrep)
*
*              For the CSF contribution we cannot skip the A,A case
*              (lack of translational invariance is taken care of by CmbnS1)
*
               If (.Not.isCSF) Then
                  If (EQ(A,RB).and. .Not.DiffOp) Go To 140
               End If
               If (.Not.DiffOp) Then
*--------------Use the translational invariance to reduce the set of
*              gradients to compute
                  Do iCar = 1, 3
                     If (IfGrad(iCar,1).and.IfGrad(iCar,2))  Then
                        IfGrad(iCar,2) = .False.
                        IndGrd(iCar,2) = -IndGrd(iCar,2)
                     End If
                  End Do
               End If
*
               If (iPrint.ge.49) Then
                  Write (6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),
     &                  i=0,nStabM-1),')'
               End If
*
               llOper = lOper(1)
               Do iComp = 2, nComp
                  llOper = iOr(llOper,lOper(iComp))
               End Do
               Call SOS(iStabO,nStabO,llOper)
               Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,
     &                  iDCRT,nDCRT)
*
*--------------Compute normalization factor due the DCR symmetrization
*              of the two basis functions and the operator.
*
               iuv = nStab(mdci)*nStab(mdcj)
               FactNd = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
               If (MolWgh.eq.1) Then
                FactNd = FactNd * DBLE(nIrrep)**2 / DBLE(iuv)
               Else If (MolWgh.eq.2) Then
                FactNd = sqrt(DBLE(iuv))*DBLE(nStabO)/DBLE(nIrrep*LmbdT)
               End If
*
               If (iPrint.ge.49) Then
                  Write (6,'(A,/,2(3F6.2,2X))')
     &                  ' *** Centers A, RB ***',
     &                  ( A(i),i=1,3), (RB(i),i=1,3)
               End If
*
*--------------Desymmetrize the matrix with which we will
*              contracte the trace.
*
               Call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                     iShell,jShell,iShll,jShll,
     &                     DAO,iPrim,jPrim,
     &                     DSOpr,nSO,nOp,FactNd)
*
*--------------Project the spherical harmonic space onto the
*              cartesian space.
*
               kk = nElem(iAng)*nElem(jAng)
               If (Transf(iShll).or.Transf(jShll)) Then
*
*-----------------ij,AB --> AB,ij
                  Call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,
     &                        iCmp*jCmp,Scr1,iCmp*jCmp)
*-----------------AB,ij --> ij,ab
                  Call SphCar(Scr1,iCmp*jCmp,iPrim*jPrim,
     &                        Scr2,nScr2,
     &                        RSph(ipSph(iAng)),
     &                        iAng,Transf(iShll),Prjct(iShll),
     &                        RSph(ipSph(jAng)),
     &                        jAng,Transf(jShll),Prjct(jShll),
     &                        DAO,kk)
               End If
               If (iPrint.ge.99) Call RecPrt(
     &                  ' Decontracted FD in the cartesian space',
     &                  ' ',DAO,iPrim*jPrim,kk)
*
*--------------Compute kappa and P.
*
               Call Setup1(Work(iExp),iPrim,Work(jExp),jPrim,
     &                     A,RB,Kappa,PCoor,ZI)
*
*--------------Compute gradients of the primitive integrals and
*              trace the result.
*
               Call Kernel(Work(iExp),iPrim,Work(jExp),jPrim,
     &                     Zeta,ZI,Kappa,Pcoor,
     &                     Final,iPrim*jPrim,
     &                     iAng,jAng,A,RB,nOrder,Krnl,
     &                     MemKer,Ccoor,nOrdOp,Grad,nGrad,
     &                     IfGrad,IndGrd,DAO,
     &                     mdci,mdcj,nOp,lOper,nComp,
     &                     iStabM,nStabM)
               If (iPrint.ge.49) Call PrGrad(' In Oneel',
     &             Grad,nGrad,lIrrep,ChDisp,5)
*
 140        Continue
          endif
            Call mma_deallocate(DSOpr)
            Call mma_deallocate(DAO)
            Call mma_deallocate(Scr2)
            Call mma_deallocate(Scr1)
            Call mma_deallocate(Final)
            Call mma_deallocate(Krnl)
 131        Continue
C        End Do
C     End Do
      End Do
*
      If (FreeiSD) Call Free_iSD()
      Call mma_deallocate(PCoor)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(ZI)
      Call mma_deallocate(Zeta)
*
      If (iPrint.ge.15) Call PrGrad(Label,Grad,nGrad,lIrrep,ChDisp,5)
*
      Return
      End
