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
      SubRoutine DrvPCM_(h1,TwoHam,D,RepNuc,nh1,First,NonEq,
     &                   Z_Nuc,Cord,MaxAto,Tessera,DMat,VTessera,
     &                   VSave,QTessera,QTessera_Slow,VSlow,nTs,Eps,
     &                   EpsInf)
      use PCM_arrays
      use Symmetry_Info, only: iChBas
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
#include "stdalloc.fh"
      Real*8, Allocatable:: FactOp(:), Integrals(:)
      Integer, Allocatable:: lOper2(:)
      Logical First
      Character*8 Label
      Logical Save_tmp, NonEq
      Real*8 RepNuc_Save
      Common /RF/RepNuc_Save
      Dimension ip(1)
#include "oneel_interface.fh"
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
        VSlow(1:nTs)=Zero
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
      call mma_allocate(FactOp,nTs,Label='FactOp')
      call mma_allocate(lOper2,nTs,Label='lOper2')
      FactOp(:)=One
      lOper2(:)=255
*
      Call Drv1_PCM(FactOp,nTs,D,nh1,Tessera,lOper2,VTessera,nOrdOp)
*
      Call mma_deallocate(lOper2)
      Call mma_deallocate(FactOp)
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
      Call mma_allocate(C_Tessera,3,nTs,Label='C_Tessera')
      nTiles=nTs
      C_Tessera(:,:)=Tessera(1:3,:)
      Call mma_allocate(Q_Tessera,nTs,Label='Q_Tessera')
      Label='<Q|V>'
*
*
      If (First) then
*
*------- PCM-integrals weighted by Q(1)
*        h1 + correction
*
         Q_Tessera(:)=QTessera(1,:)
         If (NonEq) Call DaXpY_(nTs,One,QTessera_Slow,1,Q_Tessera,1)
         Call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,
     &                        Origin,nOrdOp,rHrmt,[kOper],Integrals)
         nInt=n2Tri(lOper)
         Call CmpInt(Integrals(ip(1)),nInt,nBas,nIrrep,lOper)
         Alpha=One
         Call DaXpY_(nInt,Alpha,Integrals(ip(1)),1,h1,1)
         Call mma_deallocate(Integrals)
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
      Q_Tessera(:)=QTessera(2,:)
      Call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,Origin,
     &                     nOrdOp,rHrmt,[kOper],Integrals)
      nInt=n2Tri(lOper)
      Call CmpInt(Integrals(ip(1)),nInt,nBas,nIrrep,lOper)
      Alpha=One
      Call DaXpY_(nInt,Alpha,Integrals(ip(1)),1,TwoHam,1)
      Call mma_deallocate(Integrals)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Q_Tessera)
      Call mma_deallocate(C_Tessera)
      PrPrt=Save_tmp
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
