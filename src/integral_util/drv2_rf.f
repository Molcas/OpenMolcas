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
*              GetMem                                                  *
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
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "lundio.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), Ccoor(3),
     &       Fldxyz((lMax+1)*(lMax+2)*(lMax+3)/6), h0(nh0)
      Character ChOper(0:7)*3
      Integer   nOp(2), iStabO(0:7),
     &          iDCRR(0:7), iDCRT(0:7), iStabM(0:7)
      Logical AeqB
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
      Call GetMem('Zeta','ALLO','REAL',iZeta,m2Max)
      Call GetMem('Zeta','ALLO','REAL',ipZI ,m2Max)
      Call GetMem('Kappa','ALLO','REAL',iKappa,m2Max)
      Call GetMem('PCoor','ALLO','REAL',iPCoor,m2Max*3)
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
         iCff   = iSD( 4,iS)
         iPrim  = iSD( 5,iS)
         iExp   = iSD( 6,iS)
         iAO    = iSD( 7,iS)
         ixyz   = iSD( 8,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         x1 = Work(ixyz)
         y1 = Work(ixyz+1)
         z1 = Work(ixyz+2)
         Do jS = 1, iS
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
            x2 = Work(jxyz)
            y2 = Work(jxyz+1)
            z2 = Work(jxyz+2)
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            nSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
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
            MemKrn=MemKer*m2Max
            Call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)
*
*           Allocate memory for the final integrals all in the
*           primitive basis.
            nComp = (lMax+1)*(lMax+2)*(lMax+3)/6
            lFinal = (nComp+1) * MaxPrm(iAng) * MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
            Call GetMem('Final','ALLO','REAL',ipFnl,lFinal)
            ipFnl1 = ipFnl +nComp * MaxPrm(iAng) * MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
*
*           Scratch area for contraction step
*
            nScr1 =  Max(MaxPrm(iAng),MaxPrm(jAng)) *
     &               Max(MaxBas(iAng),MaxBas(jAng)) *
     &               nComp*nElem(iAng)*nElem(jAng)
            Call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=nComp*MaxBas(iAng)*MaxBas(jAng)
     &           *nElem(iAng)*nElem(jAng)
            Call GetMem('ScrSph','ALLO','REAL',iScrt2,nScr2)
*
*           At this point we can compute Zeta.
*           This is now computed in the ij or ji order.
*
            Call ZXia(Work(iZeta),Work(ipZI),
     &                iPrim,jPrim,Work(iExp),Work(jExp))
*
            AeqB = iS.eq.jS
*
*           Allocate memory for SO integrals that will be generated by
*           this batch of AO integrals.
*
            Call GetMem(' SO ','ALLO','REAL',ipSO,nSO*iBas*jBas)
            call dcopy_(nSO*iBas*jBas,[Zero],0,Work(ipSO),1)
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,jStab(0,mdci),nStab(mdci),
     &                                  jStab(0,mdcj),nStab(mdcj),
     &                                                iDCRR,nDCRR)
*
*           Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                 iStabM,nStabM)
*
*           Find the DCR for M and S
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               iStabO,nStabO,iDCRT,nDCRT)
*
            If (iPrint.ge.19) Then
               Write (6,*)
               Write (6,*) ' g      =',nIrrep
               Write (6,*) ' u      =',nStab(mdci)
               Write (6,'(9A)') '(U)=',(ChOper(jStab(ii,mdci)),
     &               ii = 0, nStab(mdci)-1)
               Write (6,*) ' v      =',nStab(mdcj)
               Write (6,'(9A)') '(V)=',(ChOper(jStab(ii,mdcj)),
     &               ii = 0, nStab(mdcj)-1)
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
            iuv = nStab(mdci)*nStab(mdcj)
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
            A(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*x1
            A(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*y1
            A(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*z1
            nOp(1) = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
            if(jbas.lt.-99999) write(6,*) 'nDCRR=',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
             B(1) = DBLE(iPhase(1,iDCRR(lDCRR))*
     &              iPhase(1,iDCRT(lDCRT)))*x2
             B(2) = DBLE(iPhase(2,iDCRR(lDCRR))*
     &              iPhase(2,iDCRT(lDCRT)))*y2
             B(3) = DBLE(iPhase(3,iDCRR(lDCRR))*
     &              iPhase(3,iDCRT(lDCRT)))*z2
             nOp(2) = NrOpr(iEor(iDCRT(lDCRT),iDCRR(lDCRR)),iOper,
     &                nIrrep)
             If (iPrint.ge.49) Write (6,'(A,3(3F6.2,2X))')
     &             '***** Centers A, B, & C. *****',
     &             (A(i),i=1,3),(B(i),i=1,3),(Ccoor(i),i=1,3)
*
*            Compute kappa and P.
*
             Call Setup1(Work(iExp),iPrim,Work(jExp),jPrim,
     &                   A,B,Work(iKappa),Work(iPCoor),Work(ipZI))
*
*            Compute primitive integrals. Result is ordered ij,ab.
*
             Call RFInt(Work(iExp),iPrim,Work(jExp),jPrim,
     &                   Work(iZeta),Work(ipZI),
     &                   Work(iKappa),Work(iPcoor),
     &                   Work(ipFnl),iPrim*jPrim,nComp,
     &                   iAng,jAng,A,B,nOrder,Work(iKern),
     &                   MemKer,Ccoor,lMax)
             If (iPrint.ge.49)
     &          Call RecPrt(' Primitive Integrals',' ',
     &                      Work(ipFnl),iPrim*jPrim*
     &                      nElem(iAng)*nElem(jAng),nComp)
*
*-----------Accumulate contributions due to interaction between the
*           electric field and the multipole moments.
*
            nFnc=iPrim*jPrim*nElem(iAng)*nElem(jAng)
            call dcopy_(nFnc,[Zero],0,Work(ipFnl1),1)
            Call DNaXpY(nComp,nFnc,Fldxyz,1,
     &                  Work(ipFnl),1,nFnc,
     &                  Work(ipFnl1),1,0)
            If (iPrint.ge.99) Call RecPrt(' Solvation integrals',' ',
     &                        Work(ipFnl1),iPrim*jPrim,
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
     &                      Work(iCff),iPrim,iBas)
                Call RecPrt(' Right side contraction',' ',
     &                      Work(jCff),jPrim,jBas)
             End If
*
*            Transform ij,x,ab to j,xabI
             kk=nElem(iAng)*nElem(jAng)
             Call DGEMM_('T','N',
     &                   jPrim*kk,iBas,iPrim,
     &                   1.0d0,Work(ipFnl1),iPrim,
     &                   Work(iCff),iPrim,
     &                   0.0d0,Work(iScrt1),jPrim*kk)
*            Transform j,xabI to xab,IJ
             Call DGEMM_('T','N',
     &                   kk*iBas,jBas,jPrim,
     &                   1.0d0,Work(iScrt1),jPrim,
     &                   Work(jCff),jPrim,
     &                   0.0d0,Work(ipFnl1),kk*iBas)
*
             If (iPrint.ge.99) Call
     &          RecPrt(' Contracted integrals in cartesians',' ',
     &                     Work(ipFnl1),kk,iBas*jBas)
*
*            Transform to spherical gaussians if needed.
*
             If (Transf(iShll).or.Transf(jShll)) Then
*
*             Result comes back as IJxAB or IJxAb
              call dcopy_(kk*iBas*jBas,Work(ipFnl1),1,
     &                                Work(iScrt2),1)
              Call CarSph(Work(iScrt2),kk,iBas*jBas,
     &                    Work(ipFnl1),nScr2,RSph(ipSph(iAng)),
     &                    iAng,Transf(iShll),Prjct(iShll),
     &                    RSph(ipSph(jAng)),jAng,Transf(jShll),
     &                    Prjct(jShll),Work(iScrt1),iCmp*jCmp)
             Else
*             Transpose back to IJ,x,ab
              Call DGeTmO(Work(ipFnl1),kk,kk,iBas*jBas,
     &                   Work(iScrt1),iBas*jBas)
             End If
             If (iPrint.ge.99)
     &          Call RecPrt(' Contracted Integrals in Sphericals',
     &                   ' ',Work(iScrt1),iBas*jBas,iCmp*jCmp)
*
*            At this point accumulate the batch of integrals onto the
*            final symmetry adapted integrals.
*
             If (iPrint.ge.99) Then
                Call RecPrt (' Accumulated SO integrals, so far...',
     &                               ' ',Work(ipSO),iBas*jBas,nSO)
             End If
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            mSO=MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            nIC=1
            iIC=1
            If (mSO.ne.0)
     &         Call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                     iShell,jShell,iShll,jShll,Work(iScrt1),
     &                     iBas,jBas,nIC,iIC,Work(ipSO),mSO,nOp)
*
 140        Continue
 139        Continue
*
*           Multiply with factors due to projection operators
*
           If (Fact.ne.One) Call DScal_(nSO*iBas*jBas,Fact,Work(ipSO),1)
            If (iPrint.ge.99) Then
               Write (6,*) ' Scaling SO''s', Fact
               Call RecPrt(' Final SO integrals',' ',
     &                     Work(ipSO),iBas*jBas,mSO)
            End If
*
*-----------Accumulate contribution to the Hamiltonian.
*
            iSmLbl=llOper
            If (Prprt) iSmLbl=iAnd(1,iSmLbl)
            Call SOAdd(Work(ipSO),iBas,jBas,mSO,h0,
     &                 n2Tri(iSmLbl),iSmLbl,
     &                 iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
*
            Call GetMem('  SO ','FREE','REAL',ipSO,nSO*iBas*jBas)
            Call GetMem('ScrSph','FREE','REAL',iScrt2,nScr2)
            Call GetMem('Scrtch','FREE','REAL',iScrt1,nScr1)
            Call GetMem('Final','FREE','REAL',ipFnl,lFinal)
            Call GetMem('Kernel','FREE','REAL',iKern,MemKrn)
 131        Continue
         End Do
      End Do
 100  Continue
*
      Call GetMem('PCoor','FREE','REAL',iPCoor,n2Max*3)
      Call GetMem('Kappa','FREE','REAL',iKappa,n2Max)
      Call GetMem('Zeta','FREE','REAL',ipZI ,n2Max)
      Call GetMem('Zeta','FREE','REAL',iZeta,n2Max)
*
c     Call GetMem('Drv2_RF','CHEC','REAL',iDum,iDum)
*
      Call qExit('Drv2_RF')
      Return
      End
