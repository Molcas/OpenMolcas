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
* Copyright (C) 1990-1992,1999, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv2_RF(llOper,Ccoor,nOrdOp,Fldxyz,lMax,h0,nh0)
************************************************************************
*                                                                      *
* Object: to compute the one-electron integrals. The method employed at*
*         this point is not necessarily the fastest. However, the total*
*         time for the computation of integrals will depend on the time*
*         spent in computing the two-electron integrals.               *
*         The memory at this point is assumed to be large enough to do *
*         the computation in core.                                     *
*         The data is structured with respect to four indices, two (my *
*         ny or i j) refer to primitives or basis functions and two (a *
*         b) refer to the components of the cartesian or spherical     *
*         harmonic gaussians.                                          *
*                                                                      *
* Called from: Drv1El                                                  *
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
*              SOSctt                                                  *
*              SymAd1                                                  *
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
*             Modified to reaction field calculations July  92         *
*             Modified loop structure April  99                        *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Center_Info
      use Temporary_parameters, only: PrPrt
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), Ccoor(3),
     &       Fldxyz((lMax+1)*(lMax+2)*(lMax+3)/6), h0(nh0)
      Character ChOper(0:7)*3
      Integer   nOp(2), iStabO(0:7),
     &          iDCRR(0:7), iDCRT(0:7), iStabM(0:7)
      Logical AeqB
      Real*8, Allocatable:: Zeta(:), ZI(:), Kappa(:), PCoor(:,:)
      Real*8, Allocatable:: Kern(:), Fnl(:,:), Scr1(:), Scr2(:),
     &                      SO_Int(:)
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 212
      iPrint = nPrint(iRout)
      Call qEnter('Drv2_RF')
      If (iPrint.ge.19) Then
         Write (6,*) ' In Drv2_RF: llOper'
         Write (6,'(1X,8I5)') llOper
         Write (6,*) ' In Drv2_RF: n2Tri'
         Write (6,'(1X,8I5)')  n2Tri(llOper)
      End If
*
*
      Call SOS(iStabO,nStabO,llOper)
*
*     Auxiliary memory allocation.
*
      Call mma_allocate(Zeta,S%m2Max,Label='Zeta')
      Call mma_allocate(ZI,S%m2Max,Label='ZI')
      Call mma_allocate(Kappa,S%m2Max,Label='Kappa')
      Call mma_allocate(PCoor,S%m2Max,3,Label='PCoor')
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
         If (Shells(iShll)%Aux) Go To 100
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
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
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            nSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (iPrint.ge.29) Write (6,*) ' nSO=',nSO
            If (nSO.eq.0) Go To 131
*
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &        ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*
*           Call kernel routine to get memory requirement. Observe, however
*           that kernels which will use the HRR will allocate that
*           memory internally.
*
            Call RFMem(nOrder,MemKer,iAng,jAng,nOrdOp)
*           Write(*,*)nOrder,MemKer,iAng,jAng,nOrdOp
            MemKrn=MemKer*S%m2Max
            Call mma_allocate(Kern,MemKrn,Label='Kern')
*
*           Allocate memory for the final integrals all in the
*           primitive basis.
            nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
            lFinal = S%MaxPrm(iAng) * S%MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
            Call mma_allocate(Fnl,lFinal,nComp+1,Label='Fnl')
*
*           Scratch area for contraction step
*
            nScr1 =  Max(S%MaxPrm(iAng),S%MaxPrm(jAng)) *
     &               Max(S%MaxBas(iAng),S%MaxBas(jAng)) *
     &               nComp*nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr1,nScr1,Label='Scr1')
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=nComp*S%MaxBas(iAng)*S%MaxBas(jAng)
     &           *nElem(iAng)*nElem(jAng)
            Call mma_allocate(Scr2,nScr2,Label='Scr2')
*
*           At this point we can compute Zeta.
*           This is now computed in the ij or ji order.
*
            Call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,
     &                                    Shells(jShll)%Exp)
*
            AeqB = iS.eq.jS
*
*           Allocate memory for SO integrals that will be generated by
*           this batch of AO integrals.
*
            Call mma_allocate(SO_Int,nSO*iBas*jBas,Label='SO')
            SO_Int(:)=Zero
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
*
*           Find the stabilizer for A and B
*
            Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,
     &                 dc(mdcj)%iStab,dc(mdcj)%nStab,
     &                 iStabM,nStabM)
*
*           Find the DCR for M and S
*
            Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
            If (iPrint.ge.19) Then
               Write (6,*)
               Write (6,*) ' g      =',nIrrep
               Write (6,*) ' u      =',dc(mdci)%nStab
               Write (6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),
     &               ii = 0, dc(mdci)%nStab-1)
               Write (6,*) ' v      =',dc(mdcj)%nStab
               Write (6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),
     &               ii = 0, dc(mdcj)%nStab-1)
               Write (6,*) ' LambdaR=',LmbdR
               Write (6,*) ' r      =',nDCRR
               Write (6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),
     &               ii = 0, nDCRR-1)
               Write (6,*) ' m      =',nStabM
               Write (6,'(9A)') '(M)=',(ChOper(iStabM(ii)),
     &               ii = 0, nStabM-1)
               Write (6,*) ' s      =',nStabO
               Write (6,'(9A)') '(S)=',(ChOper(iStabO(ii)),
     &               ii = 0, nStabO-1)
               Write (6,*) ' LambdaT=',LmbdT
               Write (6,*) ' t      =',nDCRT
               Write (6,'(9A)') '(R)=',(ChOper(iDCRT(ii)),
     &               ii = 0, nDCRT-1)
            End If
*
*           Compute normalization factor
*
            iuv = dc(mdci)%nStab*dc(mdcj)%nStab
            Fact = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
            If (MolWgh.eq.1) Then
               Fact = Fact * DBLE(nIrrep)**2 / DBLE(iuv)
            Else If (MolWgh.eq.2) Then
               Fact = Sqrt(DBLE(iuv))*DBLE(nStabO)/DBLE(nIrrep*LmbdT)
            End If
*
*           Loops over symmetry operations.
*
            Do 139 lDCRT = 0, nDCRT-1
            Call OA(iDCRT(lDCRT),dbsc(iCnttp)%Coor(1:3,iCnt),A)
            nOp(1) = NrOpr(iDCRT(lDCRT))
            if(jbas.lt.-99999) write(6,*) 'nDCRR=',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
             iDCRRT=iEor(iDCRR(lDCRR),iDCRT(lDCRT))
             Call OA(iDCRRT,dbsc(jCnttp)%Coor(1:3,jCnt),B)
             nOp(2) = NrOpr(iEor(iDCRT(lDCRT),iDCRR(lDCRR)))
             If (iPrint.ge.49) Write (6,'(A,3(3F6.2,2X))')
     &             '***** Centers A, B, & C. *****',
     &             (A(i),i=1,3),(B(i),i=1,3),(Ccoor(i),i=1,3)
*
*            Compute kappa and P.
*
             Call Setup1(Shells(iShll)%Exp,iPrim,
     &                   Shells(jShll)%Exp,jPrim,
     &                   A,B,Kappa,PCoor,ZI)
*
*            Compute primitive integrals. Result is ordered ij,ab.
*
             Call RFInt(Shells(iShll)%Exp,iPrim,
     &                  Shells(jShll)%Exp,jPrim,
     &                   Zeta,ZI,
     &                   Kappa,Pcoor,
     &                   Fnl,iPrim*jPrim,nComp,
     &                   iAng,jAng,A,B,nOrder,Kern,
     &                   MemKer,Ccoor,lMax)
             If (iPrint.ge.49)
     &          Call RecPrt(' Primitive Integrals',' ',
     &                      Fnl,iPrim*jPrim*
     &                      nElem(iAng)*nElem(jAng),nComp)
*
*-----------Accumulate contributions due to interaction between the
*           electric field and the multipole moments.
*
            nFnc=iPrim*jPrim*nElem(iAng)*nElem(jAng)
            call dcopy_(nFnc,[Zero],0,Fnl(1,nComp+1),1)
            Call DNaXpY(nComp,nFnc,Fldxyz,1,
     &                  Fnl,1,nFnc,
     &                  Fnl(1,nComp+1),1,0)
            If (iPrint.ge.99) Call RecPrt(' Solvation integrals',' ',
     &                        Fnl(1,nComp+1),iPrim*jPrim,
     &                        nElem(iAng)*nElem(jAng))
*
*
*------------Transform from primitive to contracted basis functions.
*            Order of transformation is fixed. It has been shown through
*            testing that the index order ij,ab will give a performance
*            that is up to 20% faster than the ab,ij index order.
*
             If (iPrint.ge.99) Then
                Call RecPrt(' Left side contraction',' ',
     &                      Shells(iShll)%pCff,iPrim,iBas)
                Call RecPrt(' Right side contraction',' ',
     &                      Shells(jShll)%pCff,jPrim,jBas)
             End If
*
*            Transform ij,x,ab to j,xabI
             kk=nElem(iAng)*nElem(jAng)
             Call DGEMM_('T','N',
     &                   jPrim*kk,iBas,iPrim,
     &                   1.0d0,Fnl(1,nComp+1),iPrim,
     &                         Shells(iShll)%pCff,iPrim,
     &                   0.0d0,Scr1,jPrim*kk)
*            Transform j,xabI to xab,IJ
             Call DGEMM_('T','N',
     &                   kk*iBas,jBas,jPrim,
     &                   1.0d0,Scr1,jPrim,
     &                         Shells(jShll)%pCff,jPrim,
     &                   0.0d0,Fnl(1,nComp+1),kk*iBas)
*
             If (iPrint.ge.99) Call
     &          RecPrt(' Contracted integrals in cartesians',' ',
     &                     Fnl(1,nComp+1),kk,iBas*jBas)
*
*            Transform to spherical gaussians if needed.
*
             If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
*
*             Result comes back as IJxAB or IJxAb
              call dcopy_(kk*iBas*jBas,Fnl(1,nComp+1),1,
     &                                Scr2,1)
              Call CarSph(Scr2,kk,iBas*jBas,
     &                    Fnl(1,nComp+1),nScr2,
     &                    RSph(ipSph(iAng)),
     &                    iAng,Shells(iShll)%Transf,
     &                         Shells(iShll)%Prjct,
     &                    RSph(ipSph(jAng)),
     &                    jAng,Shells(jShll)%Transf,
     &                         Shells(jShll)%Prjct,
     &                    Scr1,iCmp*jCmp)
             Else
*             Transpose back to IJ,x,ab
              Call DGeTmO(Fnl(1,nComp+1),kk,kk,iBas*jBas,
     &                   Scr1,iBas*jBas)
             End If
             If (iPrint.ge.99)
     &          Call RecPrt(' Contracted Integrals in Sphericals',
     &                   ' ',Scr1,iBas*jBas,iCmp*jCmp)
*
*            At this point accumulate the batch of integrals onto the
*            final symmetry adapted integrals.
*
             If (iPrint.ge.99) Then
                Call RecPrt (' Accumulated SO integrals, so far...',
     &                               ' ',SO_Int,iBas*jBas,nSO)
             End If
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            nIC=1
            iIC=1
            If (mSO.ne.0)
     &         Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                     iShell,jShell,iShll,jShll,
     &                     iAO,jAO,Scr1,
     &                     iBas,jBas,nIC,iIC,SO_Int,mSO,nOp)
*
 140        Continue
 139        Continue
*
*           Multiply with factors due to projection operators
*
           If (Fact.ne.One) Call DScal_(nSO*iBas*jBas,Fact,SO_Int,1)
            If (iPrint.ge.99) Then
               Write (6,*) ' Scaling SO''s', Fact
               Call RecPrt(' Final SO integrals',' ',
     &                     SO_Int,iBas*jBas,mSO)
            End If
*
*-----------Accumulate contribution to the Hamiltonian.
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            Call SOAdd(SO_Int,iBas,jBas,mSO,h0,
     &                 n2Tri(iSmLbl),iSmLbl,
     &                 iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
*
            Call mma_deallocate(SO_Int)
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
      Call qExit('Drv2_RF')
      Return
      End
