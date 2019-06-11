************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine DrvPCM(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
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
      iRout = 1
      iPrint = nPrint(iRout)
      Call QEnter('DrvPCM')
*                                                                      *
************************************************************************
*                                                                      *
*
*     Get the total 1st order AO density matrix
*
*     (unused?)
      Call Get_D1ao(ipD1ao,nDens)
      If (nDens.ne.nh1) Then
         Call WarningMessage(2,'DrvPCM: nDens.ne.nh1')
         Write (6,*) nDens,nh1
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate list of all atoms
*
*     Cord: list of all atoms
*
      Call Get_nAtoms_All(MaxAto)
*
      Call GetMem('Cord','Allo','Real',ipCord,3*MaxAto)
      Call GetMem('Chrg','Allo','Real',ipChrg,MaxAto)
*
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         Z = Charge(jCnttp)
         mCnt = nCntr(jCnttp)
         If (AuxCnttp(jCnttp)) mCnt = 0
         jxyz = ipCntr(jCnttp)
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            x1 = Work(jxyz)
            y1 = Work(jxyz+1)
            z1 = Work(jxyz+2)
            Do i = 0, nIrrep/nStab(ndc) - 1
               iFacx=iPhase(1,iCoset(i,0,ndc))
               iFacy=iPhase(2,iCoset(i,0,ndc))
               iFacz=iPhase(3,iCoset(i,0,ndc))
               Work(ipCord+(nc-1)*3  ) = x1*DBLE(iFacx)
               Work(ipCord+(nc-1)*3+1) = y1*DBLE(iFacy)
               Work(ipCord+(nc-1)*3+2) = z1*DBLE(iFacz)
               Work(ipChrg+(nc-1)) = Z
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('PCM-Chg','Allo','Real',ip_Q,2*nTs)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('V_Tile','Allo','Real',ip_V_Ts,2*nTs)
      Call GetMem('V_Save','Allo','Real',ip_VSave,2*nTs)
      Call GetMem('Q_Slow','Allo','Real',ip_QSlow,nTs)
      Call GetMem('V_Slow','Allo','Real',ip_VSlow,nTs)
      Call DrvPCM_(h1,TwoHam,D,RepNuc,nh1,First,NonEq,
     &             Work(ipChrg),Work(ipCord),MaxAto,
     &             Work(ip_Tess),Work(ip_DM),Work(ip_V_Ts),
     &             Work(ip_VSave),Work(ip_Q),Work(ip_QSlow),
     &             Work(ip_VSlow),nTs,Eps,EpsInf)
      Call GetMem('V_Slow','Free','Real',ip_VSlow,nTs)
      Call GetMem('Q_Slow','Free','Real',ip_QSlow,nTs)
      Call GetMem('V_Save','Free','Real',ip_VSave,2*nTs)
      Call GetMem('V_Tile','Free','Real',ip_V_Ts,2*nTs)
*                                                                      *
************************************************************************
*                                                                      *
*---- Put the current set of PCM charges on the run file.
*
      Call Put_dArray('PCM Charges',Work(ip_Q),2*nTs)
      Call GetMem('PCM-Chg','Free','Real',ip_Q,2*nTs)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('Chrg','Free','Real',ipChrg,MaxAto)
      Call GetMem('Cord','Free','Real',ipCord,3*MaxAto)
      Call GetMem('D1ao','Free','Real',ipD1ao,nDens)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('DrvPCM')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Dff)
      End
      SubRoutine DrvPCM_(h1,TwoHam,D,RepNuc,nh1,First,NonEq,
     &                   Z_Nuc,Cord,MaxAto,Tessera,DMat,VTessera,
     &                   VSave,QTessera,QTessera_Slow,VSlow,nTs,Eps,
     &                   EpsInf)
      Implicit Real*8 (A-H,O-Z)
      External PCMInt, NaMem
      Real*8 h1(nh1), TwoHam(nh1), D(nh1), Z_Nuc(MaxAto),
     &       Cord(3,MaxAto), Tessera(4,nTs), DMat(nTs,nTs),
     &       VTessera(2,nTs), VSave(2,Nts), QTessera(2,nTs),
     &       QTessera_Slow(nTs),VSlow(nTs), Origin(3)
#include "SysDef.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "pcm_pointers.fh"
      Logical First
      Character*8 Label
      Logical Save_tmp, NonEq
      Real*8 RepNuc_Save
      Common /RF/RepNuc_Save
      Dimension ip(1)
*                                                                      *
************************************************************************
*                                                                      *
*-----Reaction field a la PCM
*
*     Set up some parameters.
*
      rHrmt=One
      nOrdOp=0
      nComp=1
      lOper=255
      kOper=iChBas(1)
      Origin(1)=Zero
      Origin(2)=Zero
      Origin(3)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*---- Pick up slow charges and W's. These are present if we are doing
*     the final state.
*
      If (NonEq) Then
*
*------- Read slow components originating from the initial state
*
*        Write (*,*) 'Rd:',QTessera_Slow(1)
         Call Get_dArray('RCTFLD',QTessera_Slow,nTs)
         Call Get_dScalar('W_or_el',W_0_or_el)
         Call Get_dScalar('W_or_Inf',W_0_or_Inf)
*
*       Compute the electrostatic potential due to slow charges
*
        call dcopy_(nTs,[Zero],0,VSlow,1)
        Do iTile = 1, nTs
          XI=Tessera(1,iTile)
          YI=Tessera(2,iTile)
          ZI=Tessera(3,iTile)
          Do jTile = 1, nTs
            If (jTile.eq.iTile) Then
              Dij = 1.0694D0 * Two*Sqrt(Pi/Tessera(4,iTile))
            Else
              XJ=Tessera(1,jTile)
              YJ=Tessera(2,jTile)
              ZJ=Tessera(3,jTile)
              RIJ = Sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
              Dij = One / RIJ
            End If
            VSlow(iTile) = VSlow(iTile) + Dij * QTessera_Slow(jTile)
          End Do
        End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Evaluate the potential field on the tiles.
*
      Save_tmp=PrPrt
      PrPrt=.True.
*
*---- Do the nuclear contribution
*
      Do iTile = 1, nTs
         Call EFNuc(Tessera(1,iTile),Z_Nuc,Cord,MaxAto,
     &             VTessera(1,iTile),nOrdOp)
         VTessera(2,iTile)=Zero
      End Do
*
*---- Do the electronic contribution
*
      Call Allocate_Work(ipFactOp,nTs)
      Call Allocate_iWork(iplOper,nTs)
      call dcopy_(nTs,[One],0,Work(ipFactOp),1)
      Call ICopy(nTs,[255],0,iWork(iplOper),1)
*
      Call Drv1_PCM(Work(ipFactOp),nTs,D,nh1,Tessera,iWork(iplOper),
     &              VTessera,nOrdOp)
*
      Call Free_iWork(iplOper)
      Call Free_Work(ipFactOp)
*
*     Save the electrostatic potential
*
      call dcopy_(2*nTs,VTessera,1,VSave,1)
*
*     Add the slow charge contribution to the nuclear potential
*
      If (NonEq) Call DaXpY_(nTs,One,VSlow,1,VTessera(1,1),2)
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the charges on the cavity boundary, nuclear and
*     electronic.
*
      Call PCM_Driver(iPrint,DMat,VTessera,QTessera,nTs)
*
*---- Make the slow charges (also called orientational charges or
*     frozen charges). This is always done regardless if they ever
*     will be used.  The charges can be used in non-equilibrium
*     calculations.
*
      If (EpsInf.gt.Zero.and..Not.NonEq) Then
        Fact=(Eps-EpsInf)/(Eps-One)
        call dcopy_(nTs,    QTessera(1,1),2,QTessera_Slow,1)
        Call DaXpY_(nTs,One,QTessera(2,1),2,QTessera_Slow,1)
        Call DScal_(nTs,Fact,QTessera_Slow,1)
*
*      Compute the electrostatic potential due to slow charges
*
       call dcopy_(nTs,[Zero],0,VSlow,1)
       Do iTile = 1, nTs
         XI=Tessera(1,iTile)
         YI=Tessera(2,iTile)
         ZI=Tessera(3,iTile)
         Do jTile = 1, nTs
           If (jTile.eq.iTile) Then
             Dij = 1.0694D0 * Two*Sqrt(Pi/Tessera(4,iTile))
           Else
             XJ=Tessera(1,jTile)
             YJ=Tessera(2,jTile)
             ZJ=Tessera(3,jTile)
             RIJ = Sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
             Dij = One / RIJ
           End If
           VSlow(iTile) = VSlow(iTile) + Dij * QTessera_Slow(jTile)
         End Do
       End Do

      End If
*
*     Recover the nuclear potential (discarding the possible variations
*     due to slow charges)
*
      call dcopy_(nTs,VSave(1,1),2,VTessera(1,1),2)
*
      W_or_el     = Zero
      W_or_nuc    = Zero
      W_or_Inf    = Zero
      W_or_InfNuc = zero
      Do iTile = 1, nTs
        W_or_el  = W_or_el  + QTessera_Slow(iTile)*VTessera(2,iTile)
        W_or_nuc = W_or_nuc + QTessera_Slow(iTile)*VTessera(1,iTile)
        If(NonEq) then
          QInf = QTessera(1,iTile) + QTessera(2,iTile)
        Else
          Fact = (Eps - Epsinf) / (Eps - One)
          QInf = (QTessera(1,iTile) + QTessera(2,iTile)) * (One - Fact)
        EndIf
        W_or_Inf = W_or_Inf + QInf * VSlow(iTile)
        W_or_InfNuc = W_or_InfNuc + QTessera(1,iTile) * VSlow(iTile)
      End Do
*
*---- Write out the slow charges and constants if this is an
*     equilibrium calculations.
*
      If (EpsInf.gt.Zero.and..Not.NonEq) Then
*        Write (*,*) 'Wr:',QTessera_Slow(1)
         Call Put_dArray('RCTFLD',QTessera_Slow,nTs)
         Call Put_dScalar('W_or_el',W_or_el)
         Call Put_dScalar('W_or_Inf',W_or_Inf)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Now add terms to RepNuc
*
*     Interaction terms:
*     nuclei-nuclear solvation charge       (ENN)
*     nuclei-electronic solvation charge    (ENE)
*     electrons-nuclear solvation charge    (EEN)
*     electrons-electronic solvation charge (EEE)
*
      ENN = Zero
      ENE = Zero
      EEN = Zero
      EEE = Zero
      Do iTile = 1, nTs
        ENN = ENN + QTessera(1,iTile) * VTessera(1,iTile)
        ENE = ENE + QTessera(1,iTile) * VTessera(2,iTile)
        EEN = EEN + QTessera(2,iTile) * VTessera(1,iTile)
        EEE = EEE + QTessera(2,iTile) * VTessera(2,iTile)
      End Do
      If (First) then
         RepNuc = RepNuc + Half * ENN
         If(NonEq)
     &     RepNuc = RepNuc
     &            + Half * W_or_nuc
     &            + Half * W_or_InfNuc
     &            - Half * W_0_or_el
     &            - Half * W_0_or_Inf
         RepNuc_Save=RepNuc
         Label='PotNuc00'
         Call Put_Temp(Label,[RepNuc],1)
      End If

*                                                                      *
************************************************************************
*                                                                      *
*     Now add terms to h1 and TwoHam!
*
      Call Allocate_Work(ipC,3*nTs)
      nTiles=nTs
      call dcopy_(nTiles,Tessera(1,1),4,Work(ipC  ),3)
      call dcopy_(nTiles,Tessera(2,1),4,Work(ipC+1),3)
      call dcopy_(nTiles,Tessera(3,1),4,Work(ipC+2),3)
      Call Allocate_Work(ipQ,nTs)
      Label='<Q|V>'
*
*
      If (First) then
*
*------- PCM-integrals weighted by Q(1)
*        h1 + correction
*
         call dcopy_(nTs,QTessera,2,Work(ipQ),1)
         If (NonEq) Call DaXpY_(nTs,One,QTessera_Slow,1,Work(ipQ),1)
         Call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,
     &                        Origin,nOrdOp,rHrmt,[kOper])
         nInt=n2Tri(lOper)
         Call CmpInt(Work(ip(1)),nInt,nBas,nIrrep,lOper)
         Alpha=One
         Call DaXpY_(nInt,Alpha,Work(ip(1)),1,h1,1)
         Call Free_Work(ip(1))
*
*------  Save the modified h1 matrix
*
         Label='h1_raw  '
         Call Put_Temp(Label,h1,nh1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- PCM-integrals weighted by Q(2)
*     TwoHam + correction
*
      call dcopy_(nTs,QTessera(2,1),2,Work(ipQ),1)
      Call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,Origin,
     &                     nOrdOp,rHrmt,[kOper])
      nInt=n2Tri(lOper)
      Call CmpInt(Work(ip(1)),nInt,nBas,nIrrep,lOper)
      Alpha=One
      Call DaXpY_(nInt,Alpha,Work(ip(1)),1,TwoHam,1)
      Call Free_Work(ip(1))
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipQ)
      Call Free_Work(ipC)
      PrPrt=Save_tmp
c     Call PrMtrx('h1',1,1,ip_of_Work(h1),Work)
c     Call PrMtrx('TwoHam',1,1,ip_of_Work(TwoHam),Work)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
