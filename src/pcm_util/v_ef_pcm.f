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
************************************************************************
      Subroutine V_EF_PCM(nAt,nTs,DoPot,DoFld,AtmC,Tessera,V,EF_n,EF_e)
      Implicit Real*8 (a-h,o-z)
      Dimension AtmC(3,nAt)
      Dimension V(*),EF_n(3,*),EF_e(3,*),Tessera(4,*)
      Logical DoPot, DoFld
*
*---- Compute potential on tesserae
*
      If(DoPot) then
        Call FZero(V,nTs)
        nOrdOp = 0
        Call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
      EndIf
*
*---- Compute electric field on tesserae
*
      If(DoFld) then
        Call FZero(EF_n,3*nTs)
        Call FZero(EF_e,3*nTs)
        nOrdOp = 1
        Call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
      EndIf
      Return
      End
c----------------------------------------------------------------------
      Subroutine Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
      Implicit Real*8 (A-H,O-Z)
      Dimension V(*),EF_n(3,*),EF_e(3,*)
      Dimension Tessera(4,*)
      Dimension AtmC(3,*)
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 Temp(3)
      Real*8, Allocatable :: Chrg(:)
      Logical Found
      Real*8, Allocatable :: D1ao(:)
*                                                                      *
      Call mma_allocate(Chrg,nat)
      Call Get_dArray('Nuclear charge',Chrg,nAt)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the desired multipole on the tiles
*
*     1) The nuclear contribution
*
      Do iTs = 1, nTs
         Call EFNuc(Tessera(1,iTs),Chrg,AtmC,nAt,Temp,nOrdOp)
         If(nOrdOp.eq.0) then
           V(iTs) = Temp(1)
         ElseIf(nOrdOp.eq.1) then
           EF_n(1,iTs) = Temp(1)
           EF_n(2,iTs) = Temp(2)
           EF_n(3,iTs) = Temp(3)
         EndIf
      End Do
*
      Call mma_deallocate(Chrg)
*
*     2) The electronic contribution
*
*
*     Get the total 1st order AO density matrix
*
      Call Qpg_dArray('D1ao',Found,nDens)
      If (Found .and. nDens/=0) Then
         Call mma_allocate(D1ao,nDens,Label='D1ao')
      Else
         Write (6,*) 'Mlt_pcm: D1ao not found.'
         Call Abend()
      End If
      Call Get_D1ao(D1ao,nDens)
*
      Call Allocate_Work(ipFactOp,nTs)
      Call Allocate_iWork(iplOper,nTs)
      call dcopy_(nTs,[One],0,Work(ipFactOp),1)
      Call ICopy(nTs,[255],0,iWork(iplOper),1)
*
      Call drv_ef_PCM(Work(ipFactOp),nTs,D1ao,nDens,
     &              Tessera,iWork(iplOper),
     &              EF_e,nOrdOp)
      If(nOrdOp.eq.0) then
        Do iTs = 1, nTs
          V(iTs) = Ef_e(1,iTs)
        EndDo
      EndIf
*
      Call Free_iWork(iplOper)
      Call Free_Work(ipFactOp)
      Call Free_Work(D1ao)
*
      Return
      End
*------------------------------------------------------------------------
      SubRoutine drv_ef_PCM(FactOp,nTs,FD,nFD,CCoor,lOper,VTessera,
     &                       nOrdOp)
************************************************************************
*                                                                      *
* Object: to compute the local multipole moment, desymmetrize the 1st  *
*         order density matrix and accumulate contributions to the     *
*         global multipole expansion.                                  *
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
      use Center_Info
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
      Real*8 A(3), B(3), C(3), FD(nFD), FactOp(nTs), CCoor(4,nTs),
     &       RB(3), TRB(3), TA(3), VTessera(3,nTs)
      Character ChOper(0:7)*3
      Integer   lOper(nTs), iStabO(0:7),
     &          iDCRR(0:7), iDCRT(0:7), iStabM(0:7), nOp(3)
      Logical AeqB
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
*
      iIrrep = 0
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
      Call Nr_Shells(nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      Do iS = 1, nSkal
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
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (nSO.eq.0) Go To 131
            If (iPrint.ge.19) Write (6,'(A,A,A,A,A)')
     &        ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
*
*           Call kernel routine to get memory requirement.
*
            Call EFMmP(nOrder,MemKer,iAng,jAng,nOrdOp)
*           Write (*,*)nOrder,MemKer,iAng,jAng,nOrdOp
            MemKrn=MemKer*S%m2Max
            Call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)
*
*           Allocate memory for the final integrals, all in the
*           primitive basis.
*
            nComp = (nOrdOp+1)*(nOrdOp+2)/2
            lFinal = S%MaxPrm(iAng) * S%MaxPrm(jAng)
     &             * nElem(iAng)*nElem(jAng)
     &             * nComp
            Call GetMem('Final','ALLO','REAL',ipFnl,lFinal)
*
*           Scratch area for contraction step
*
            nScr1 =  S%MaxPrm(iAng)*S%MaxPrm(jAng) *
     &               nElem(iAng)*nElem(jAng)
            Call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)
*
*           Scratch area for the transformation to spherical gaussians
*
            nScr2=S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
            Call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)
*
            nDAO =iPrim*jPrim*nElem(iAng)*nElem(jAng)
            Call GetMem(' DAO ','Allo','Real',ipDAO,nDAO)
*
*           At this point we can compute Zeta.
*
            Call ZXia(Work(iZeta),Work(ipZI),
     &                iPrim,jPrim,Shells(iShll)%Exp,
     &                            Shells(jShll)%Exp)
*
            AeqB = iS.eq.jS
*
*           Find the DCR for A and B
*
            Call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,
     &                     dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
            If (iPrint.ge.49) Write (6,'(10A)')
     &         ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
*
*-----------Find the stabilizer for A and B
*
            Call Inter(dc(mdci)%iStab,dc(mdci)%nStab,
     &                 dc(mdcj)%iStab,dc(mdcj)%nStab,
     &                 iStabM,nStabM)
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
            Do lDCRR = 0, nDCRR-1
               Call OA(iDCRR(lDCRR),B,RB)
*
*--------------Loop over operators
*
               Do 5000 iTile = 1, nTs
                  If (FactOp(iTile).eq.Zero) Go To 5000
                  call dcopy_(3,Ccoor(1,iTile),1,C,1)
*
*-----------------Generate stabilizer of the operator.
*
                  Call SOS(iStabO,nStabO,lOper(iTile))
*
*-----------------Find the DCR for M and S
*
                  Call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,
     &                     iDCRT,nDCRT)
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
                  iuv = dc(mdci)%nStab*dc(mdcj)%nStab
                  FactNd = DBLE(iuv*nStabO) / DBLE(nIrrep**2*LmbdT)
                  If (MolWgh.eq.1) Then
                     FactNd = FactNd * DBLE(nIrrep)**2 / DBLE(iuv)
                  Else If (MolWgh.eq.2) Then
                     FactNd = Sqrt(DBLE(iuv))*DBLE(nStabO) /
     &                        DBLE(nIrrep*LmbdT)
                  End If
                  FactNd = FactNd * FactOp(iTile)
*
                  Do lDCRT = 0, nDCRT-1
                     nOp(1) = NrOpr(iDCRT(lDCRT))
                     nOp(2) = NrOpr(iEor(iDCRT(lDCRT),iDCRR(lDCRR)))
                     nOp(3) = NrOpr(0)

                     Call OA(iDCRT(lDCRT),A,TA)
                     Call OA(iDCRT(lDCRT),RB,TRB)
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
     &                           iAO,jAO,Work(ipDAO),iPrim,jPrim,
     &                           Work(ipDSOp),nSO,nOp,FactNd)
*
*--------------------Project the spherical harmonic space onto the
*                    cartesian space.
*
                     kk = nElem(iAng)*nElem(jAng)
                     If (Shells(iShll)%Transf.or.
     &                   Shells(jShll)%Transf) Then
*
*-----------------------ij,AB --> AB,ij
                        Call DGeTmO(Work(ipDAO),iPrim*jPrim,iPrim*jPrim,
     &                              iCmp*jCmp,Work(iScrt1),iCmp*jCmp)
*-----------------------AB,ij --> ij,ab
                        Call SphCar(Work(iScrt1),iCmp*jCmp,iPrim*jPrim,
     &                              Work(iScrt2),nScr2,
     &                              RSph(ipSph(iAng)),
     &                              iAng,Shells(iShll)%Transf,
     &                                   Shells(iShll)%Prjct,
     &                              RSph(ipSph(jAng)),
     &                              jAng,Shells(jShll)%Transf,
     &                                   Shells(jShll)%Prjct,
     &                              Work(ipDAO),kk)
                     End If
                     If (iPrint.ge.99) Call RecPrt(
     &                        ' Decontracted FD in the cartesian space',
     &                        ' ',Work(ipDAO),iPrim*jPrim,kk)
*
*--------------------Compute kappa and P.
*
                     Call Setup1(Shells(iShll)%Exp,iPrim,
     &                           Shells(jShll)%Exp,jPrim,
     &                   TA,TRB,Work(iKappa),Work(iPCoor),Work(ipZI))
*
*
*--------------------Compute the potential at a tessera.
*
cpcm_solvent
c        write(6,*)'Work(iExp),iPrim,Work(jExp),jPrim'
c        write(6,*)Shells(iShll)%Exp(1),iPrim,Shells(jShll)%Exp(1),jPrim
c        write(6,*)'Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor)'
c        write(6,*)Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor)
c        write(6,*)'Work(ipFnl),iPrim*jPrim,nComp,iAng,jAng,norder'
c        write(6,*)Work(ipFnl),iPrim*jPrim,nComp,iAng,jAng,norder
cpcm_solvent end
                     Call EFPrm(Shells(iShll)%Exp,iPrim,
     &                          Shells(jShll)%Exp,jPrim,
     &                          Work(iZeta),Work(ipZI),
     &                          Work(iKappa),Work(iPcoor),
     &                          Work(ipFnl),iPrim*jPrim,nComp,
     &                          iAng,jAng,TA,TRB,nOrder,Work(iKern),
     &                          MemKer,C,nOrdOp)
                     If (iPrint.ge.49) Call RecPrt(' Final Integrals',
     &                                 ' ',Work(ipFnl),nDAO,nComp)
*
*--------------------Trace with 1st order density matrix and accumulate
*                    to the potenital at tessera iTile
*
                     If (iPrint.ge.49) Call RecPrt(
     &                        ' Decontracted FD in the cartesian space',
     &                        ' ',Work(ipDAO),nDAO,1)
                     ipFnlc=ipFnl
                     Do iComp = 1, nComp
                        If (iPrint.ge.49)
     &                     Call RecPrt('VTessera(iComp,iTile)',' ',
     &                                  VTessera(iComp,iTile),1,1)
                        VTessera(iComp,iTile)=
     &                      VTessera(iComp,iTile) +
     &                      DDot_(nDAO,Work(ipDAO),1,Work(ipFnlc),1)
                        If (iPrint.ge.49)
     &                     Call RecPrt('VTessera(iComp,iTile)',' ',
     &                                  VTessera(iComp,iTile),1,1)
                        ipFnlc=ipFnlc+nDAO
                     End Do
*
                  End Do
 5000          Continue
            End Do
*
            Call GetMem('DSOpr ','Free','REAL',ipDSOp,nSO*iPrim*jPrim)
            Call GetMem(' DAO ','Free','Real',ipDAO,iPrim*jPrim*
     &                  nElem(iAng)*nElem(jAng))
            Call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
            Call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
            Call GetMem('Final','Free','Real',ipFnl,lFinal)
            Call GetMem('Kernel','Free','Real',iKern,MemKrn)
 131        Continue
         End Do
      End Do
*
      Call GetMem('PCoor','FREE','REAL',iPCoor,S%m2Max*3)
      Call GetMem('Kappa','FREE','REAL',iKappa,S%m2Max)
      Call GetMem('Zeta','FREE','REAL',ipZI ,S%m2Max)
      Call GetMem('Zeta','FREE','REAL',iZeta,S%m2Max)
*
c     Call GetMem('drv_ef_PCM','CHEC','REAL',iDum,iDum)
*
      Return
      End
