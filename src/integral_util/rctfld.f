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
* Copyright (C) 1992,2002, Roland Lindh                                *
************************************************************************
      SubRoutine RctFld(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
************************************************************************
*                                                                      *
*     Driver for RctFld_                                               *
*                                                                      *
************************************************************************
      use PCM_arrays, only: MM
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
      Logical First, Dff, NonEq
*
      nComp=(lMax+1)*(lMax+2)*(lMax+3)/6
      Call GetMem('Vs','Allo','Real',ipVs,nComp*2)
      Call GetMem('QV','Allo','Real',ipQV,nComp*2)
*
      Call RctFld_(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,
     &             MM,nComp,Work(ipVs),Work(ipQV))
*
      Call GetMem('QV','Free','Real',ipQV,nComp*2)
      Call GetMem('Vs','Free','Real',ipVs,nComp*2)
*
      Return
      End
      SubRoutine RctFld_(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,
     &                   Q_solute,nComp,Vs,QV)
************************************************************************
*                                                                      *
* Object: to apply a modification to the one-electron hamiltonian due  *
*         the reaction field. The code here is direct!                 *
*         This subroutine works only if a call to GetInf has been      *
*         prior to calling this routine.                               *
*                                                                      *
*         h1: one-electron hamiltonian to be modified. Observe that    *
*             the contribution due to the reaction field is added to   *
*             this array, i.e. it should be set prior to calling this  *
*             routine.                                                 *
*                                                                      *
*         TwoHam: dito two-electron hamiltonian.                       *
*                                                                      *
*         D:  the first order density matrix                           *
*                                                                      *
*         h1, TwoHam and D are all in the SO basis.                    *
*                                                                      *
*         Observe the energy expression for the electric field -       *
*         charge distribution interaction!                             *
*                                                                      *
*         -1/2 Sum(nl) E(tot,nl)M(tot,nl)                              *
*                                                                      *
* Called from: DrvRF                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              DCopy  (ESSL)                                           *
*              RecPrt                                                  *
*              MltNuc                                                  *
*              Drv1                                                    *
*              AppFld                                                  *
*              DDot_  (ESSL)                                           *
*              Drv2                                                    *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             July '92                                                 *
*                                                                      *
*             Modified for nonequilibrum calculations January 2002 (RL)*
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1), Origin(3)
      Character*72 Label
      Character*8 Label2
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
      Logical First, Dff, NonEq
      Real*8 Q_solute(nComp,2), Vs(nComp,2), QV(nComp,2)
      Dimension FactOp(1), lOper(1)
*
*-----Statement Functions
*
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
*
      iRout = 1
      iPrint = nPrint(iRout)
      Call qEnter('RctFld')
*
      lOper(1)=1
      nOrdOp=lMax
*-----Set flag so only the diagonal blocks are computed
      Prprt=.True.
      Call FZero(Origin,3)
*
*-----Generate local multipoles in the primitive basis and accumulate to
*     global multipoles.
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*-----Add nuclear-nuclear contribution to the potential energy
*
      If (First) Then
*
         If (NonEq) Then
            Call Get_dArray('RCTFLD',QV,nComp*2)
            Call Get_dScalar('E_0_NN',E_0_NN)
         End If
*1)
*
*------- Compute M(nuc,nl), nuclear multipole moments
*
         Do iMax = 0, lMax
            ip = 1+iOff(iMax)
            Call RFNuc(Origin,Q_solute(ip,1),iMax)
         End Do

         if(lXF) Then
*
*------- Add contribution from XFIELD multipoles
*
*        Use Vs as temporary space, it will anyway be overwritten
            Call XFMoment(lMax,Q_solute,Vs,nComp,Origin)
         EndIf

         If (iPrint.ge.99) Call RecPrt('Nuclear Multipole Moments',
     &                                 ' ',Q_solute(1,1),1,nComp)
*
*--------Solve dielectical equation(nuclear contribution), i.e.
*        M(nuc,nl) -> E(nuc,nl)
*
         call dcopy_(nComp,Q_solute(1,1),1,Vs(1,1),1)
         Call AppFld(Vs(1,1),rds,Eps,lMax,EpsInf,NonEq)
*
         If (iPrint.ge.99) Call RecPrt('Nuclear Electric Field',
     &                                 ' ',Vs(1,1),1,nComp)
*
*--------Vnn = Vnn - 1/2 Sum(nl) E(nuc,nl)*M(nuc,nl)
*
         RepNuc = PotNuc -
     &            Half * DDot_(nComp,Q_solute(1,1),1,Vs(1,1),1)
*
*------- Add contibutions due to slow counter charges
*
         If (NonEq) RepNuc=RepNuc+E_0_NN
         If (iPrint.ge.99) Write (6,*) ' RepNuc=',RepNuc
*2)
*
*--------Compute contibution to the one-electron hamiltonian
*
*        hpq = hpq + Sum(nl) E(nuc,nl)*<p|M(nl)|q>
*
         If (iPrint.ge.19) Then
            Write (6,*) 'h1'
            lOff = 1
            Do iIrrep = 0, nIrrep-1
               n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
               If (n.gt.0) Then
                  Write (Label,'(A,I1)')
     &             'Diagonal Symmetry Block ',iIrrep+1
                  Call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
                  lOff = lOff + n
               End If
            End Do
         End If
*
*------- Add potential due to slow counter charges
*
         If (NonEq) Then
            call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
            Call AppFld_NonEQ_2(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
            Call DaXpY_(nComp,One,QV(1,2),1,Vs(1,1),1)
         End If
*
         Call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(1,1),lMax,h1,nh1)
*
         If (iPrint.ge.19) Then
            Write (6,*) 'h1(mod)'
            lOff = 1
            Do iIrrep = 0, nIrrep-1
               n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
               If (n.gt.0) Then
                  Write (Label,'(A,I1)')
     &             'Diagonal Symmetry Block ',iIrrep+1
                  Call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
                  lOff = lOff + n
               End If
            End Do
         End If
*
*------- Update h1 and RepNuc_save with respect to static contributions!
*
         Label2='PotNuc00'
         Call Put_Temp(Label2,[RepNuc],1)
         Label2='h1_raw  '
         Call Put_Temp(Label2,h1,nh1)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the electronic contribution to the charge distribution.
*
*3)
*     M(el,nl) =  - Sum(p,q) Dpq <p|M(nl)|q>
*
      nOpr=1
      FactOp(1)=One
*-----Reset array for storage of multipole moment expansion
      call dcopy_(nComp,[Zero],0,Q_solute(1,2),1)
      Do iMltpl = 1, lMax
         Do ix = iMltpl, 0, -1
            If (Mod(ix,2).eq.0) Then
               iSymX=1
            Else
               ixyz=1
               iSymX=2**IrrFnc(ixyz)
               If (Origin(1).ne.Zero) iSymX = iOr(iSymX,1)
            End If
            Do iy = iMltpl-ix, 0, -1
               If (Mod(iy,2).eq.0) Then
                  iSymY=1
               Else
                  ixyz=2
                  iSymY=2**IrrFnc(ixyz)
                  If (Origin(2).ne.Zero) iSymY = iOr(iSymY,1)
               End If
               iz = iMltpl-ix-iy
               If (Mod(iz,2).eq.0) Then
                  iSymZ=1
               Else
                  ixyz=4
                  iSymZ=2**IrrFnc(ixyz)
                  If (Origin(3).ne.Zero) iSymZ = iOr(iSymZ,1)
               End If
*
               iTemp = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
     &                            nIrrep),nIrrep)
               lOper(1)=iOr(lOper(1),iTemp)
            End Do
         End Do
      End Do
      If (iPrint.ge.19) Then
         Write (6,*) '1st order density'
         lOff = 1
         Do iIrrep = 0, nIrrep-1
            n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
            Write (Label,'(A,I1)')
     &       'Diagonal Symmetry Block ',iIrrep+1
            Call Triprt(Label,' ',D(lOff),nBas(iIrrep))
            lOff = lOff + n
         End Do
      End If
*
      Call Drv1_RF(FactOp,nOpr,D,nh1,Origin,lOper,Q_solute(1,2),lMax)
*
      If (iPrint.ge.99) Call RecPrt('Electronic Multipole Moments',
     &                              ' ',Q_solute(1,2),1,nComp)
*
*-----Solve dielectical equation(electronic contribution), i.e.
*     M(el,nl) -> E(el,nl)
*
      call dcopy_(nComp,Q_solute(1,2),1,Vs(1,2),1)
      Call AppFld(Vs(1,2),rds,Eps,lMax,EpsInf,NonEq)
      If (iPrint.ge.99) Call RecPrt('Electronic Electric Field',
     &                              ' ',Vs(1,2),1,nComp)
*4)
*
*-----Compute contribution to the two-electron hamiltonian.
*
*     T(D)pq = T(D)pq + Sum(nl) E(el,nl)*<p|M(nl)|q>
*
      Call Drv2_RF(lOper(1),Origin,nOrdOp,Vs(1,2),lMax,TwoHam,nh1)
*
      If (iPrint.ge.19) Then
         Write (6,*) 'h1(mod)'
         lOff = 1
         Do iIrrep = 0, nIrrep-1
            n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
            If (n.gt.0) Then
               Write (Label,'(A,I1)')
     &          'Diagonal Symmetry Block ',iIrrep+1
               Call Triprt(Label,' ',h1(lOff),nBas(iIrrep))
               lOff = lOff + n
            End If
         End Do
         Write (6,*) 'TwoHam(mod)'
         lOff = 1
         Do iIrrep = 0, nIrrep-1
            n = nBas(iIrrep)*(nBas(iIrrep)+1)/2
            Write (Label,'(A,I1)')
     &       'Diagonal Symmetry Block ',iIrrep+1
            Call Triprt(Label,' ',TwoHam(lOff),nBas(iIrrep))
            lOff = lOff + n
         End Do
         Write (6,*) ' RepNuc=',RepNuc
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Write information to be used for gradient calculations or for
*     non-equilibrium calculations
*
      If (.Not.NonEq) Then
*
*        Save total solute multipole moments and total potential
*        of the solution.
*
         call dcopy_(nComp,Q_solute(1,1),1,QV(1,1),1)
         call daxpy_(nComp,One,Q_solute(1,2),1,QV(1,1),1)
         call dcopy_(nComp,Vs(1,1),1,QV(1,2),1)
         call daxpy_(nComp,One,Vs(1,2),1,QV(1,2),1)
         Call Put_dArray('RCTFLD',QV,nComp*2)
*
*        Compute terms to be added to RepNuc for non-equilibrium
*        calculation.
*
         call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
         Call AppFld_NonEQ_1(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
         E_0_NN=-Half*DDot_(nComp,QV(1,1),1,QV(1,2),1)
*
         call dcopy_(nComp,QV(1,1),1,QV(1,2),1)
         Call AppFld_NonEQ_2(QV(1,2),rds,Eps,lMax,EpsInf,NonEq)
         E_0_NN=E_0_NN+DDot_(nComp,Q_solute(1,1),1,QV(1,2),1)
         Call Put_dScalar('E_0_NN',E_0_NN)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Call qExit('RctFld')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Dff)
      End
