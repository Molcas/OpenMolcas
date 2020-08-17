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
*               1994, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Dot1El(Kernel,KrnlMm,Hess,nHess,DiffOp,CCoor,
     &                 FD,nFD,lOper,nComp,Label)
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
* Called from: Seward                                                  *
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
*              SOGthr                                                  *
*              DesymD                                                  *
*              DScal    (ESSL)                                         *
*              TriPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January '90                                              *
*             Anders Bernhardsson  Dec '94                             *
*                                                                      *
*             Modified for general kernel routines January '91         *
*             Modified for nonsymmetrical operators February '91       *
*             Modified for gradients October '91                       *
*             Modified for Hermite-Gauss quadrature November '90       *
*             Modified for Rys quadrature November '90                 *
*             Modified for multipole moments November '90              *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
c#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), Ccoor(3,nComp), FD(nFD),
     &       RB(3), Hess(nHess)
      Character Label*80
c     Character ChOper(0:7)*3
      Integer iDCRR(0:7), iDCRT(0:7), iStabM(0:7),iCoM(0:7,0:7),
     &          IndHss(0:1,0:2,0:1,0:2,0:7), nOp(2),
     &           iStabO(0:7),lOper(nComp),IndGrd(0:2,0:1,0:7)
      Logical AeqB, EQ, TstFnc, DiffOp, IfHss(0:1,0:2,0:1,0:2),
     &        TF,Chck,ifgrd(0:2,0:1)
c     Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      itri(i1,i2)=MAX(i1,i2)*(MAX(i1,i2)-1)/2+MIN(i1,i2)
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
*
      call dcopy_(nHess,[Zero],0,Hess,1)
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
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*     Double loop over shells.
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
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
C        Do jS = 1, iS
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
C        write(6,*)
C    &  'iShll,iAng,iCmp,iBas,iPrim,iAO,ixyz,mdci,iShell'
C        write(6,*) (iSD(i,iS),i=0,11)
C        write(6,*)
C    &  'jShll,jAng,jCmp,jBas,jPrim,jAO,jxyz,mdcj,jShell'
C        write(6,*) (iSD(i,jS),i=0,11)
*
*       Call kernel routine to get memory requirement.
*
        nOrdOp = 0  ! not used in this implementation
        Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
        MemKrn=MemKer*m2Max
        Call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)
*
*       Allocate memory for the final integrals, all in the
*       primitive basis.
*
        lFinal = 21 * MaxPrm(iAng) * MaxPrm(jAng) *
     &           nElem(iAng)*nElem(jAng)
        Call GetMem('Final','ALLO','REAL',ipFnl,lFinal)
*
*       Scratch area for contraction step
*
        nScr1 =  MaxPrm(iAng)*MaxPrm(jAng) *
     &           nElem(iAng)*nElem(jAng)
        Call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)
*
*       Scratch area for the transformation to spherical gaussians
*
        nScr2=MaxPrm(iAng)*MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
        Call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)
*
          Call GetMem(' DAO ','Allo','Real',ipDAO,
     &                iPrim*jPrim*nElem(iAng)*nElem(jAng))
*
*         At this point we can compute Zeta.
*
          Call ZXia(Work(iZeta),Work(ipZI),
     &              iPrim,jPrim,Shells(iShll)%Exp,
     &                          Shells(jShll)%Exp)
*
            AeqB = iS.eq.jS

*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,iOper,nIrrep,
     &               jStab(0,mdci),nStab(mdci),
     &               jStab(0,mdcj),nStab(mdcj),iDCRR,nDCRR)
            If (.Not.DiffOp .and. nDCRR.eq.1 .and. EQ(A,B)) Go To 131
c           If (iPrint.ge.49) Write (6,'(10A)')
c    &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
*-----------Find the stabilizer for A and B
*
            Call Inter(jStab(0,mdci),nStab(mdci),
     &                 jStab(0,mdcj),nStab(mdcj),
     &                 iStabM,nStabM)
*
*          Generate all possible (left) CoSet
*          To the stabil. of A and B
*
            Do 433 i = 0, nIrrep-1
               Do 434 j = 0, nStabM-1
                  iCoM(i,j) = iEor(iOper(i),iStabM(j))
 434            Continue
 433         Continue
*           Order the Coset so we will have the unique ones first
            nMax = 1
            Do 435 j = 1, nIrrep-1
*              Check uniqueness
               Do 436 i = 0, nMax - 1
                  Do 437 ielem = 0, nStabM-1
                     If (iCoM(i,1).eq.iCoM(j,ielem))
     &                  Go To 435
 437               Continue
 436            Continue
*              Move unique CoSet
               nMax = nMax + 1
               Do 438 ielem = 0, nStabM-1
                  iTmp = iCoM(nMax-1,ielem)
                  iCoM(nMax-1,ielem) = iCoM(j,ielem)
                  iCoM(j,ielem) = iTmp
 438            Continue
               If (nMax.eq.nIrrep/nStabM) Go To 439
 435         Continue
 439         Continue
             nCoM=nIrrep/nStabM
             Call LCopy(36,[.false.],0,IfHss,1)
             Do 400 iAtom=0,1
               Do 410 iCar=0,2
                 Do 420 jAtom=0,iAtom
                   if (iAtom.eq.jAtom) Then
                      iStop=iCar
                   Else
                      iStop=2
                   End If
                   Do 430 jCar=0,iStop
                      iComp1=2**iCar
                      iComp2=2**jCar
                      iComp=iEOr(iComp1,iComp2)
                      Chck=TstFnc(iOper,nIrrep,iCoM,
     &                       nCoM,iChTbl,0,iComp,
     &                       nStabM)
                      If (Chck)
     &                  IfHss(iAtom,iCar,jAtom,jCar)=.true.
 430            Continue
 420          Continue
 410        Continue
 400       Continue
           Call ICopy(nirrep*36,[0],0,Indhss(0,0,0,0,0),1)
           Call ICopy(nirrep*6,[0],0,indgrd,1)
*
*          Determine which displacement in all IR's, each center is associated with
*
           nnIrrep=nIrrep
           If (sIrrep) nnIrrep=1

           Do 200   iIrrep=0,nnIrrep-1
             nDisp1 = IndDsp(mdci,iIrrep)
             nDisp2 = IndDsp(mdcj,iIrrep)
             Do 201 iCar = 0,2
              iComp = 2**iCar
              If ( TF(mdci,iIrrep,iComp)) Then
                 nDisp1 = nDisp1 + 1
                 IndGrd(iCar,0,iIrrep) = nDisp1
                 if (iIrrep.eq.0) IfGrd(iCar,0) = .True.
              Else
                 if (iIrrep.eq.0) IfGrd(iCar,0) = .False.
                 IndGrd(iCar,0,iIrrep)=0
              End If
              iComp = 2**iCar
              If ( TF(mdcj,iIrrep,iComp)) Then
                 nDisp2 = nDisp2 + 1
                 IndGrd(iCar,1,iIrrep) = nDisp2
                 if (iIrrep.eq.0) IfGrd(iCar,1) = .True.
              Else
                 IndGrd(iCar,1,iIrrep)=0
                 if (iIrrep.eq.0) IfGrd(iCar,1) = .True.
              End If
 201         Continue
 200       Continue
*
*          Determine index for each 2'nd derivative
*
*
           Do  iIrrep=0,nIrrep-1
             Do iAtom=0,1
              Do iCar=0,2
               Do jAtom=0,iAtom
               if (iAtom.eq.jAtom) Then
                 istop=iCar
               Else
                 iStop=2
                End If
                Do jCar=0,istop
                 If ((IndGrd(iCar,iAtom,iIrrep).gt.0).and.
     &               (IndGrd(jCar,jAtom,iIrrep).gt.0))
     &           Then
                   IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &                itri(IndGrd(iCar,iAtom,iIrrep),
     &                     IndGrd(jCar,jAtom,iIrrep))
                  Else
                      IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=0
                  End If
                 End Do
               End Do
              End Do
             End Do
           End Do
           If (.not.DiffOp) Then
               iAtom=1
               Do 440 iCar=0,2
               Do 441 jAtom=0,1
                   If (iAtom.eq.jAtom) Then
                     iStop=iCar
                   Else
                     iStop=2
                   End If
                   Do 460 jCar=0,iStop
                    If (IfHss(0,iCar,0,jCar).or.
     &                    IfHss(0,jCar,0,iCar)) Then
                    IfHss(iAtom,iCar,jAtom,jCar)=.false.
                       Do 445 iIrrep=0,nIrrep-1
                           IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &                          -IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
 445                  Continue
                    End If
 460               Continue
 441           Continue
 440           Continue
*
            End If
*           Allocate memory for the elements of the Fock or 1st order
*           denisty matrix which are associated with the current shell
*           pair.
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            If (nSO.eq.0) Go To 131
            Call GetMem('DSOpr ','ALLO','REAL',ipDSOp,nSO*iPrim*jPrim)
            Call GetMem('DSO ','ALLO','REAL',ipDSO,nSO*iPrim*jPrim)
            call dcopy_(nSO*iPrim*jPrim,[Zero],0,Work(ipDSO),1)
            call dcopy_(nSO*iPrim*jPrim,[Zero],0,Work(ipDSOp),1)
*
*           Gather the elements from 1st order density / Fock matrix.
*
            Call SOGthr(Work(ipDSO),iBas,jBas,nSO,FD,
     &                  n2Tri(iSmLbl),iSmLbl,
     &                  iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)
*
*           Project the Fock/1st order density matrix in AO
*           basis on to the primitive basis.
*
c           If (iPrint.ge.99) Then
c              Call RecPrt(' Left side contraction',' ',
c    &                     Shells(iShll)%pCff,iPrim,iBas)
c              Call RecPrt(' Right side contraction',' ',
c    &                     Shells(jShll)%pCff,jPrim,jBas)
c           End If
*
*           Transform IJ,AB to J,ABi
            Call DGEMM_('T','T',
     &                  jBas*nSO,iPrim,iBas,
     &                  1.0d0,Work(ipDSO),iBas,
     &                        Shells(iShll)%pCff,iPrim,
     &                  0.0d0,Work(ipDSOp),jBas*nSO)
*           Transform J,ABi to AB,ij
            Call DGEMM_('T','T',
     &                  nSO*iPrim,jPrim,jBas,
     &                  1.0d0,Work(ipDSOp),jBas,
     &                        Shells(jShll)%pCff,jPrim,
     &                  0.0d0,Work(ipDSO),nSO*iPrim)
*           Transpose to ij,AB
            Call DGeTmO(Work(ipDSO),nSO,nSO,iPrim*jPrim,Work(ipDSOp),
     &                  iPrim*jPrim)
            Call GetMem('DSO ','Free','Real',ipDSO,nSO*iBas*jBas)
*
c           If (iPrint.ge.99) Call
c    &         RecPrt(' Decontracted 1st order density/Fock matrix',
c    &                ' ',Work(ipDSOp),iPrim*jPrim,nSO)
*
*           Loops over symmetry operations.
*
            nOp(1) = NrOpr(0,iOper,nIrrep)
            if(jBas.lt.-999999) write(6,*) 'gcc overoptimization',nDCRR
            Do 140 lDCRR = 0, nDCRR-1
               RB(1)  = DBLE(iPhase(1,iDCRR(lDCRR)))*B(1)
               RB(2)  = DBLE(iPhase(2,iDCRR(lDCRR)))*B(2)
               RB(3)  = DBLE(iPhase(3,iDCRR(lDCRR)))*B(3)
               nOp(2) = NrOpr(iDCRR(lDCRR),iOper,nIrrep)
               If (EQ(A,RB).and. (.Not.DiffOp)) Go To 140
*
c              If (iPrint.ge.49) Then
c                 Write (6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),
c    &                  i=0,nStabM-1),')'
c              End If
*
               llOper = lOper(1)
               Do 5000 iComp = 2, nComp
                  llOper = iOr(llOper,lOper(iComp))
 5000          Continue
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
                  FactNd = sqrt(DBLE(iuv))*DBLE(nStabO)
     &                   / DBLE(nIrrep*LmbdT)
               End If
*
c              If (iPrint.ge.49) Then
c                 Write (6,'(A,/,2(3F6.2,2X))')
c    &            ' *** Centers A, RB ***',
c    &            ( A(i),i=1,3), (RB(i),i=1,3)
c              End If
*
*--------------Desymmetrize the matrix with which we will
*              contracte the trace.
*
               Call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,
     &                     iShell,jShell,iShll,jShll,
     &                     Work(ipDAO),iPrim,jPrim,
     &                     Work(ipDSOp),nSO,nOp,FactNd)
*
*--------------Project the spherical harmonic space onto the
*              cartesian space.
*
               kk = nElem(iAng)*nElem(jAng)
               If (Shells(iShll)%Transf.or.Shells(jShll)%Transf) Then
*
*-----------------ij,AB --> AB,ij
                  Call DGeTmO(Work(ipDAO),iPrim*jPrim,iPrim*jPrim,
     &                        iCmp*jCmp,Work(iScrt1),iCmp*jCmp)
*-----------------AB,ij --> ij,ab
                  Call SphCar(Work(iScrt1),iCmp*jCmp,iPrim*jPrim,
     &                        Work(iScrt2),nScr2,
     &                        RSph(ipSph(iAng)),iAng,
     &                        Shells(iShll)%Transf,
     &                        Shells(iShll)%Prjct,
     &                        RSph(ipSph(jAng)),jAng,
     &                        Shells(jShll)%Transf,
     &                        Shells(jShll)%Prjct,
     &                        Work(ipDAO),kk)
               End If
c              If (iPrint.ge.99) Call RecPrt(
c    &                  ' Decontracted FD in the cartesian space',
c    &                  ' ',Work(ipDAO),iPrim*jPrim,kk)
*
*--------------Compute kappa and P.
*
               Call Setup1(Shells(iShll)%Exp,iPrim,
     &                     Shells(jShll)%Exp,jPrim,
     &                     A,RB,Work(iKappa),Work(iPCoor),Work(ipZI))
*
*--------------Compute gradients of the primitive integrals and
*              trace the result.
*
*
CBS            write(6,*) 'Call the  Kernel'
*
               Call Kernel(Shells(iShll)%Exp,iPrim,
     &                     Shells(jShll)%Exp,jPrim,
     &                     Work(iZeta),Work(ipZI),
     &                     Work(iKappa),Work(iPcoor),
     &                     Work(ipFnl),iPrim*jPrim,
     &                     iAng,jAng,A,RB,nOrder,Work(iKern),
     &                     MemKer*iPrim*jPrim,Ccoor,
     &                     nOrdOp,Hess,nHess,
     &                     IfHss,IndHss,ifgrd,indgrd,Work(ipDAO),
     &                     mdci,mdcj,nOp,lOper,nComp,
     &                     iStabM,nStabM)

#ifdef _DEBUG_
               write(6,*)  'Hess after Kernel call in dot1el '
               Call HssPrt(Hess,nHess)
#endif
*
 140        Continue
*
            Call GetMem('DSOpr ','Free','REAL',ipDSOp,nSO*iPrim*jPrim)
 131        Continue
          Call GetMem(' DAO ','Free','Real',ipDAO,iPrim*jPrim*
     &                nElem(iAng)*nElem(jAng))
        Call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
        Call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
*
        Call GetMem('Final','Free','Real',ipFnl,lFinal)
        Call GetMem('Kernel','Free','Real',iKern,MemKrn)
*
C        End Do
C     End Do
      End Do
*
      Call Free_iSD()
*
      Call GetMem('Kappa','FREE','REAL',iKappa,n2Max)
      Call GetMem('PCoor','FREE','REAL',iPCoor,n2Max*3)
      Call GetMem('Zeta','FREE','REAL',ipZI ,n2Max)
      Call GetMem('Zeta','FREE','REAL',iZeta,n2Max)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_character(Label)
      End
