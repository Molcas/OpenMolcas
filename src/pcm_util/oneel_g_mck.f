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
      SubRoutine OneEl_g_mck(Kernel,KrnlMm,Grad,nGrad,DiffOp,CCoor,FD,
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
      use Basis_Info
      use Center_Info
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "angtp.fh"
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
CNIKO      Real*8 A(3), B(3), Ccoor(3,nComp), FD(nFD),
      Real*8 A(3), B(3), Ccoor(*), FD(nFD),
     &       RB(3), Grad(nGrad)
      Character ChOper(0:7)*3, Label*80
      Integer iDCRR(0:7), iDCRT(0:7), iStabM(0:7),
     &          IndGrd(3,2), nOp(2), iStabO(0:7), lOper(nComp)
      Logical AeqB, EQ, DiffOp, IfGrad(3,3)
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
      call dcopy_(nGrad,[Zero],0,Grad,1)
*
*     Auxiliary memory allocation.
*
      Call GetMem('Zeta','ALLO','REAL',iZeta,S%m2Max)
      Call GetMem('Zeta','ALLO','REAL',ipZI ,S%m2Max)
      Call GetMem('Kappa','ALLO','REAL',iKappa,S%m2Max)
      Call GetMem('PCoor','ALLO','REAL',iPCoor,S%m2Max*3)
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal)
      Call Setup_iSD()
*                                                                      *
************************************************************************
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
         iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)

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
*
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (nSO.eq.0) Go To 131
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
            If (.Not.DiffOp .and. nDCRR.eq.1 .and. EQ(A,B)) Go To 131
            If (iPrint.ge.49) Write (6,'(10A)')
     &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &         ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*
*           Call kernel routine to get memory requirement.
*
            Call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
            MemKrn=MemKer*S%m2Max
            Call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)
*
*           Allocate memory for the final integrals, all in the
*           primitive basis.
*
            lFinal = 6 * S%MaxPrm(iAng) * S%MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng) * nComp
            Call GetMem('Final','ALLO','REAL',ipFnl,lFinal)
*
*           Scratch area for contraction step
*
            nScr1 =S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
            Call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
            Call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)
*
            Call GetMem(' DAO ','Allo','Real',ipDAO,iPrim*jPrim*
     &                  nElem(iAng)*nElem(jAng))
*
*           At this point we can compute Zeta.
*
            Call ZXia(Work(iZeta),Work(ipZI),
     &                iPrim,jPrim,Shells(iShll)%Exp,
     &                            Shells(jShll)%Exp)
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
            Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,
     &                 dc(mdcj)%iStab,dc(mdcj)%nStab,
     &                             iStabM,nStabM)
*
*           Allocate memory for the elements of the Fock or 1st order
*           denisty matrix which are associated with the current shell
*           pair.
*
            Call GetMem('DSOpr ','ALLO','REAL',ipDSOp,nSO*iPrim*jPrim)
            Call GetMem('DSO ','ALLO','REAL',ipDSO,nSO*iPrim*jPrim)
*
*           Gather the elements from 1st order density / Fock matrix.
*
            Call SOGthr(Work(ipDSO),iBas,jBas,nSO,FD,
     &                  n2Tri(iSmLbl),iSmLbl,
     &                  iCmp,jCmp,iShell,jShell,
     &                  AeqB,iAO,jAO)
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
            If (iPrint.ge.99) Call
     &         RecPrt(' Decontracted 1st order density/Fock matrix',
     &                ' ',Work(ipDSOp),iPrim*jPrim,nSO)
*
*           Loops over symmetry operations.
*
            nOp(1) = NrOpr(0)
c VV: gcc bug: one has to use this if!
          if(nDCRR.ge.1) then
            Do 140 lDCRR = 0, nDCRR-1
               Call OA(iDCRR(lDCRR),B,RB)
               nOp(2) = NrOpr(iDCRR(lDCRR))
               If (EQ(A,RB).and. .Not.DiffOp) Go To 140
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
               Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
*
*--------------Compute normalization factor due the DCR symmetrization
*              of the two basis functions and the operator.
*
               iuv = dc(mdci)%nStab*dc(mdcj)%nStab
               FactNd = DBLE(iuv*nStabO) / DBLE(nIrrep**2 * LmbdT)
               If (MolWgh.eq.1) Then
                  FactNd = FactNd * DBLE(nIrrep)**2 / DBLE(iuv)
               Else If (MolWgh.eq.2) Then
                  FactNd = sqrt(DBLE(iuv))
     &                    * DBLE(nStabO)/DBLE(nIrrep*LmbdT)
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
     &                     iAO,jAO,Work(ipDAO),iPrim,jPrim,
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
     &                        RSph(ipSph(iAng)),
     &                        iAng,Shells(iShll)%Transf,
     &                             Shells(iShll)%Prjct,
     &                        RSph(ipSph(jAng)),
     &                        jAng,Shells(jShll)%Transf,
     &                             Shells(jShll)%Prjct,
     &                        Work(ipDAO),kk)
               End If
               If (iPrint.ge.99) Call RecPrt(
     &                  ' Decontracted FD in the cartesian space',
     &                  ' ',Work(ipDAO),iPrim*jPrim,kk)
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
               Call Kernel(Shells(iShll)%Exp,iPrim,
     &                     Shells(jShll)%Exp,jPrim,
     &                     Work(iZeta),Work(ipZI),
     &                     Work(iKappa),Work(iPcoor),
     &                     Work(ipFnl),iPrim*jPrim,
     &                     iAng,jAng,A,RB,nOrder,Work(iKern),
     &                     MemKer,Ccoor,nOrdOp,Grad,nGrad,
     &                     IfGrad,IndGrd,Work(ipDAO),
     &                     mdci,mdcj,nOp,lOper,nComp,
     &                     iStabM,nStabM)
               If (iPrint.ge.49) Call PrGrad_mck(' In Oneel',
     &             Grad,nGrad,ChDisp,5)
*
 140        Continue
          endif
            Call GetMem('DSOpr ','Free','REAL',ipDSOp,nSO*iPrim*jPrim)
            Call GetMem(' DAO ','Free','Real',ipDAO,iPrim*jPrim*
     &                nElem(iAng)*nElem(jAng))
            Call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
            Call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
            Call GetMem('Final','Free','Real',ipFnl,lFinal)
            Call GetMem('Kernel','Free','Real',iKern,MemKrn)
 131        Continue
C        End Do
C     End Do
      End Do
*
      Call Free_iSD()
      Call GetMem('Kappa','FREE','REAL',iKappa,S%m2Max)
      Call GetMem('PCoor','FREE','REAL',iPCoor,S%m2Max*3)
      Call GetMem('Zeta','FREE','REAL',ipZI ,S%m2Max)
      Call GetMem('Zeta','FREE','REAL',iZeta,S%m2Max)
*
      If (iPrint.ge.15)Call PrGrad_mck(Label,Grad,nGrad,ChDisp,5)
*
      Return
      End
