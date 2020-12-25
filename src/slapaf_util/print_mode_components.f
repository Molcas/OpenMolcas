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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
* Print_Mode_Components
*
*> @brief
*> Print the contributions from the primitive internal coordinates to the
*> vibrational modes.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Compute and print the contributions from the primitive internal coordinates
*> (stretches, angles, dihedrals, out-of-planes) to the vibrational modes
*>
*> @param[in] Modes Vibrational modes, as computed by e.g. ::freqanal
*> @param[in] Freq Vibrational frequencies
*> @param[in] nModes Number of modes
*> @param[in] lModes Size of \p Modes
*> @param[in] lDisp Number of displacements per irrep
************************************************************************
      Subroutine Print_Mode_Components(Modes,Freq,nModes,lModes,lDisp)
      use Symmetry_Info, only: nIrrep
      use Slapaf_Info, only: Cx, Gx, Gx0, NAC, Q_nuclear, dMass, Coor,
     &                       Grd, Weights, ANr, Shift, GNrm, Lambda,
     &                       Energy, Energy0, DipM, MF, qInt, dqInt,
     &                       RefGeo, BM, iBM, dBM, idBM, nqBM, BMx,
     &                       Degen, nStab, jStab, iCoSet, AtomLbl,
     &                       Smmtrc, Lbl, mRowH
      use Slapaf_Parameters, only: iRow, iRow_c, iInt, nFix,
     &                             ddV_Schlegel, HWRS, iOptH, HUpMet,
     &                             HrmFrq_Show, IRC, Curvilinear,
     &                             Redundant, FindTS, nBVec, nDimBC,
     &                             User_Def, Analytic_Hessian, MaxItr,
     &                             iOptC, UpMeth, HSet, BSet, rHidden,
     &                             CnstWght, PrQ, lOld, Numerical, Beta,
     &                             Beta_Disp, Line_Search,
     &                             TSConstraints, GNrm_Threshold, Mode,
     &                             ThrEne, ThrGrd, nLambda, iRef,
     &                             ThrCons, ThrMEP, Baker, eMEPTest,
     &                             rMEP, MEP, nMEP, MEPNum, MEPCons,
     &                             dMEPStep, MEP_Type, MEP_Algo,
     &                             Header, Max_Center, mTROld, RtRnc,
     &                             Delta, rFuzz, lNmHss, Cubic,
     &                             Request_Alaska, Request_RASSI,
     &                             lOld_Implicit, CallLast, lSoft,
     &                             lCtoF, Track, TwoRunFiles, isFalcon
      use thermochem
      Implicit None
#include "backup_info.fh"
#include "print.fh"
#include "stdalloc.fh"
#include "real.fh"
      Real*8 rDum(1)
      Real*8 :: Modes(*), Freq(*), Mx, MinComp
      Integer :: LuInput, iRout,
     &           iPrint, nX, i, j, nB, iq, nAll_Atoms, nUnique_Atoms,
     &           iB, lDisp(nIrrep), nModes, lModes, LuIC, ii, im, nK,
     &           iErr, PLback, nQQ, nsAtom
      Real*8, Dimension(:,:), Allocatable :: KMtrx, KTrsp, KKtB, IntMod,
     &                                       NMod
      Integer, Dimension(:), Allocatable :: Sort
      Integer, Parameter :: nLbl=10*MxAtom
      Integer, External :: IsFreeUnit, iPrintLevel, AixRm
      Logical :: Cartesian, Found
      Character(Len=8) :: Filename
      Character(Len=16) :: StdIn
      Character(Len=24), Allocatable :: Label(:)
      Character(Len=180), External :: Get_Ln_EOF
      Character(LEN=180):: Line
      Real*8, External :: DDot_

      Real*8, Allocatable:: Bk_Energy(:)
      Real*8, Allocatable:: Bk_Energy0(:)
      Real*8, Allocatable:: Bk_DipM(:,:)
      Real*8, Allocatable:: Bk_GNrm(:)
      Real*8, Allocatable:: Bk_Cx(:,:,:)
      Real*8, Allocatable:: Bk_Gx(:,:,:)
      Real*8, Allocatable:: Bk_Gx0(:,:,:)
      Real*8, Allocatable:: Bk_MF(:,:)
      Real*8, Allocatable:: Bk_Lambda(:,:)

      Real*8, Allocatable:: Bk_NAC(:,:)
      Real*8, Allocatable:: Bk_Q_nuclear(:)
      Real*8, Allocatable:: Bk_dMass(:)
      Real*8, Allocatable:: Bk_Coor(:,:)
      Real*8, Allocatable:: Bk_Grd(:,:)
      Real*8, Allocatable:: Bk_Weights(:)
      Real*8, Allocatable:: Bk_Shift(:,:)
      Real*8, Allocatable:: Bk_qInt(:,:)
      Real*8, Allocatable:: Bk_dqInt(:,:)
      Real*8, Allocatable:: Bk_RefGeo(:,:)
      Real*8, Allocatable:: Bk_BMx(:,:)
      Real*8, Allocatable:: Bk_Degen(:,:)

      Real*8, Allocatable:: Bk_BM(:)
      Real*8, Allocatable:: Bk_dBM(:)
      Integer, Allocatable:: Bk_iBM(:)
      Integer, Allocatable:: Bk_idBM(:)
      Integer, Allocatable:: Bk_nqBM(:)

      Integer, Allocatable:: Bk_ANr(:)
      Integer, Allocatable:: Bk_jStab(:,:)
      Integer, Allocatable:: Bk_nStab(:)
      Integer, Allocatable:: Bk_iCoSet(:,:)
      Character(LEN=LENIN), Allocatable:: Bk_AtomLbl(:)
      Logical, Allocatable:: Bk_Smmtrc(:,:)
      Character(LEN=8), Allocatable:: Bk_Lbl(:)
      Integer, Allocatable:: Bk_mRowH(:)
*
*
*---- Ugly hack: backup all "global" slapaf variables in case this is
*                called from inside slapaf
*
*     Note, this routine might be called from outside the Slapaf
*     environment. In which case there is no backup to be made.
*
      nsAtom=SIZE(Cx,2)
      If (Allocated(Cx)) Then
         Call mma_allocate(Bk_Cx,3,nsAtom,MaxItr+1,Label='Bk_Cx')
         Bk_Cx(:,:,:) = Cx(:,:,:)
         Call mma_deallocate(Cx)
      End If
      If (Allocated(Gx)) Then
         Call mma_allocate(Bk_Gx,3,nsAtom,MaxItr+1,Label='Bk_Gx')
         Bk_Gx(:,:,:) = Gx(:,:,:)
         Call mma_deallocate(Gx)
      End If
      If (Allocated(Gx0)) Then
         Call mma_allocate(Bk_Gx0,3,nsAtom,MaxItr+1,Label='Bk_Gx0')
         Bk_Gx0(:,:,:) = Gx0(:,:,:)
         Call mma_deallocate(Gx0)
      End If
      If (Allocated(NAC)) Then
         Call mma_allocate(Bk_NAC,3,nsAtom,Label='Bk_NAC')
         Bk_NAC(:,:) = NAC(:,:)
         Call mma_deallocate(NAC)
      End If
      If (Allocated(Q_nuclear)) Then
         Call mma_allocate(Bk_Q_nuclear,nsAtom,Label='Bk_Q_nuclear')
         Bk_Q_nuclear(:) = Q_nuclear(:)
         Call mma_deallocate(Q_nuclear)
      End If
      If (Allocated(dMass)) Then
         Call mma_allocate(Bk_dMass,nsAtom,Label='Bk_dMass')
         Bk_dMass(:) = dMass(:)
         Call mma_deallocate(dMass)
      End If
      If (Allocated(Coor)) Then
         Call mma_allocate(Bk_Coor,3,nsAtom,Label='Bk_Coor')
         Bk_Coor(:,:) = Coor(:,:)
         Call mma_deallocate(Coor)
      End If
      If (Allocated(Grd)) Then
         Call mma_allocate(Bk_Grd,3,nsAtom,Label='Bk_Grd')
         Bk_Grd(:,:) = Grd(:,:)
         Call mma_deallocate(Grd)
      End If
      If (Allocated(ANr)) Then
         Call mma_allocate(Bk_ANr,nsAtom,Label='Bk_ANr')
         Bk_ANr(:) = ANr(:)
         Call mma_deallocate(ANr)
      End If
      If (Allocated(Weights)) Then
         Call mma_allocate(Bk_Weights,SIZE(Weights),Label='Bk_Weights')
         Bk_Weights(:) = Weights(:)
         Call mma_deallocate(Weights)
      End If
      If (Allocated(Shift)) Then
         Call mma_allocate(Bk_Shift,SIZE(Shift,1),MaxItr,
     &                     Label='Bk_Shift')
         Bk_Shift(:,:) = Shift(:,:)
         Call mma_deallocate(Shift)
      End If
      If (Allocated(GNrm)) Then
         Call mma_allocate(Bk_GNrm,MaxItr+1,Label='Bk_GNrm')
         Bk_GNrm(:) = GNrm(:)
         Call mma_deallocate(GNrm)
      End If
      If (Allocated(Lambda)) Then
         Call mma_allocate(Bk_Lambda,SIZE(Lambda,1),MaxItr+1,
     &                     Label='Bk_Lambda')
         Bk_Lambda(:,:) = Lambda(:,:)
         Call mma_deallocate(Lambda)
      End If
      If (Allocated(Energy)) Then
         Call mma_allocate(Bk_Energy,MaxItr+1,Label='Bk_Energy')
         Bk_Energy(:) = Energy(:)
         Call mma_deallocate(Energy)
      End If
      If (Allocated(Energy0)) Then
         Call mma_allocate(Bk_Energy0,MaxItr+1,Label='Bk_Energy0')
         Bk_Energy0(:) = Energy0(:)
         Call mma_deallocate(Energy0)
      End If
      If (Allocated(MF)) Then
         Call mma_allocate(Bk_MF,3,nsAtom,Label='Bk_MF')
         Bk_MF(:,:) = MF(:,:)
         Call mma_deallocate(MF)
      End If
      If (Allocated(DipM)) Then
         Call mma_allocate(Bk_DipM,3,MaxItr+1,Label='Bk_DipM')
         Bk_DipM(:,:) = DipM(:,:)
         Call mma_deallocate(DipM)
      End If
      If (Allocated(qInt)) Then
         Call mma_allocate(Bk_qInt,SIZE(qInt,1),MaxItr,
     &                     Label='Bk_qInt')
         Bk_qInt(:,:) = qInt(:,:)
         Call mma_deallocate(qInt)
      End If
      If (Allocated(dqInt)) Then
         Call mma_allocate(Bk_dqInt,SIZE(dqInt,1),MaxItr,
     &                     Label='Bk_dqInt')
         Bk_dqInt(:,:) = dqInt(:,:)
         Call mma_deallocate(dqInt)
      End If
      If (Allocated(RefGeo)) Then
         Call mma_allocate(Bk_RefGeo,3,nsAtom,Label='Bk_RefGeo')
         Bk_RefGeo(:,:) = RefGeo(:,:)
         Call mma_deallocate(RefGeo)
      End If
      If (Allocated(BM)) Then
         Call mma_allocate(Bk_BM,SIZE(BM),Label='Bk_BM')
         Bk_BM(:) = BM(:)
         Call mma_deallocate(BM)
      End If
      If (Allocated(dBM)) Then
         Call mma_allocate(Bk_dBM,SIZE(dBM),Label='Bk_dBM')
         Bk_dBM(:) = dBM(:)
         Call mma_deallocate(dBM)
      End If
      If (Allocated(iBM)) Then
         Call mma_allocate(Bk_iBM,SIZE(iBM),Label='Bk_iBM')
         Bk_iBM(:) = iBM(:)
         Call mma_deallocate(iBM)
      End If
      If (Allocated(idBM)) Then
         Call mma_allocate(Bk_idBM,SIZE(idBM),Label='Bk_idBM')
         Bk_idBM(:) = idBM(:)
         Call mma_deallocate(idBM)
      End If
      If (Allocated(nqBM)) Then
         Call mma_allocate(Bk_nqBM,SIZE(nqBM),Label='Bk_nqBM')
         Bk_nqBM(:) = nqBM(:)
         Call mma_deallocate(nqBM)
      End If
      If (Allocated(BMx)) Then
         Call mma_allocate(Bk_BMx,SIZE(BMx,1),SIZE(BMx,2),
     &                     Label='Bk_BMx')
         Bk_BMx(:,:) = BMx(:,:)
         Call mma_deallocate(BMx)
      End If
      If (Allocated(Degen)) Then
         Call mma_allocate(Bk_Degen,SIZE(Degen,1),SIZE(Degen,2),
     &                     Label='Bk_Degen')
         Bk_Degen(:,:) = Degen(:,:)
         Call mma_deallocate(Degen)
      End If
      If (Allocated(jStab)) Then
         Call mma_allocate(Bk_jStab,[0,7],[1,SIZE(jStab,2)],
     &                     Label='Bk_jStab')
         Bk_jStab(:,:) = jStab(:,:)
         Call mma_deallocate(jStab)
      End If
      If (Allocated(iCoSet)) Then
         Call mma_allocate(Bk_iCoSet,[0,7],[1,SIZE(iCoSet,2)],
     &                     Label='Bk_iCoSet')
         Bk_iCoSet(:,:) = iCoSet(:,:)
         Call mma_deallocate(iCoSet)
      End If
      If (Allocated(nStab)) Then
         Call mma_allocate(Bk_nStab,SIZE(nStab),Label='Bk_nStab')
         Bk_nStab(:) = nStab(:)
         Call mma_deallocate(nStab)
      End If
      If (Allocated(AtomLbl)) Then
         Call mma_allocate(Bk_AtomLbl,SIZE(AtomLbl),Label='Bk_AtomLbl')
         Bk_AtomLbl(:) = AtomLbl(:)
         Call mma_deallocate(AtomLbl)
      End If
      If (Allocated(Smmtrc)) Then
         Call mma_allocate(Bk_Smmtrc,3,SIZE(Smmtrc,2),Label='Bk_Smmtrc')
         Bk_Smmtrc(:,:) = Smmtrc(:,:)
         Call mma_deallocate(Smmtrc)
      End If
      If (Allocated(Lbl)) Then
         Call mma_allocate(Bk_Lbl,SIZE(Lbl),Label='Bk_Lbl')
         Bk_Lbl(:) = Lbl(:)
         Call mma_deallocate(Lbl)
      End If
      If (Allocated(mRowH)) Then
         Call mma_allocate(Bk_mRowH,SIZE(mRowH),Label='Bk_mRowH')
         Bk_mRowH(:) = mRowH(:)
         Call mma_deallocate(mRowH)
      End If

      Bk_Header(:)=Header(:)
      Bk_iRef=iRef
      Bk_NmIter=NmIter
      Bk_iter=iter
      Bk_nDimBC=nDimBC
      Bk_MxItr=MxItr
      Bk_Max_Center=Max_Center
      Bk_iOptC=iOptC
      Bk_mode=mode
      Bk_mTROld=mTROld
      Bk_nWndw=nWndw
      Bk_iOptH=iOptH
      Bk_nLambda=nLambda
      Bk_nMEP=nMEP
      Bk_nBVec=nBVec
      Bk_IRC=IRC
      Bk_mTtAtm=mTtAtm
      Bk_MEPnum=MEPnum
      Bk_RootMap(:)=RootMap(:)
      Bk_Stop=Stop
      Bk_lOld=lOld
      Bk_CurviLinear=CurviLinear
      Bk_HSet=HSet
      Bk_BSet=BSet
      Bk_lNmHss=lNmHss
      Bk_Cubic=Cubic
      Bk_Baker=Baker
      Bk_DDV_Schlegel=DDV_Schlegel
      Bk_Line_Search=Line_Search
      Bk_HWRS=HWRS
      Bk_Analytic_Hessian=Analytic_Hessian
      Bk_FindTS=FindTS
      Bk_MEP=MEP
      Bk_User_Def=User_Def
      Bk_rMEP=rMEP
      Bk_lOld_Implicit=lOld_Implicit
      Bk_HrmFrq_Show=HrmFrq_Show
      Bk_eMEPtest=eMEPtest
      Bk_Redundant=Redundant
      Bk_lCtoF=lCtoF
      Bk_lSoft=lSoft
      Bk_CallLast=CallLast
      Bk_TwoRunFiles=TwoRunFiles
      Bk_TSConstraints=TSConstraints
      Bk_MEPCons=MEPCons
      Bk_Track=Track
      Bk_Request_Alaska=Request_Alaska
      Bk_Request_RASSI=Request_RASSI
      Bk_cMass(:)=cMass(:)
      Bk_ThrEne=ThrEne
      Bk_ThrGrd=ThrGrd
      Bk_Beta=Beta
      Bk_Beta_Disp=Beta_Disp
      Bk_Delta=Delta
      Bk_Rtrnc=Rtrnc
      Bk_rHidden=rHidden
      Bk_ThrCons=ThrCons
      Bk_GNrm_Threshold=GNrm_Threshold
      Bk_CnstWght=CnstWght
      Bk_dMEPStep=dMEPStep
      Bk_rFuzz=rFuzz
      Bk_ThrMEP=ThrMEP
      Bk_lTherm=lTherm
      Bk_lDoubleIso=lDoubleIso
      Bk_nUserPT=nUserPT
      Bk_nsRot=nsRot
      Bk_UserT(:)=UserT(:)
      Bk_UserP=UserP
      Bk_HUpMet=HUpMet
      Bk_UpMeth=UpMeth
      Bk_MEP_Type=MEP_Type
      Bk_MEP_Algo=MEP_Algo
      Bk_isFalcon=isFalcon
      Bk_mB_Tot=mB_Tot
      Bk_mdB_Tot=mdB_Tot
      Bk_mq=mq
      Bk_iSBS=iSBS
      Bk_WeightedConstraints=WeightedConstraints
      Bk_NADC=NADC
      Bk_EDiffZero=EDiffZero
      Bk_ApproxNADC=ApproxNADC
      Bk_iState(:)=iState(:)
      Bk_iRow = iRow
      Bk_iRow_c=iRow_c
      Bk_nFix = nFix
      Bk_iInt = iInt
*
      iRout = 55
      iPrint=nPrint(iRout)
*
*---- Make a backup of the runfile, since we are going to change the
*     internal coordinates definition.
*
      Call fCopy('RUNFILE','RUNBCK2',iErr)
      If (iErr.ne.0) Call Abend()
*
*---- Remove the Hessian and disable translational and rotational
*     invariance
*
      Call Put_dArray('Hess',rDum,0)
      Call Get_iScalar('System BitSwitch',iSBS)
      iSBS=iOr(iSBS,2**7)
      iSBS=iOr(iSBS,2**8)
      Call Put_iScalar('System BitSwitch',iSBS)
*                                                                      *
************************************************************************
* Call Slapaf to build the B matrix and get the displacement vectors   *
* corresponding to the primitive internal coordinates (bonds and       *
* angles)                                                              *
************************************************************************
*                                                                      *
*---- Process the input
*
      PLback=iPrintLevel(-1)
      i=iPrintLevel(0)
      LuInput=21
      LuInput=IsFreeUnit(LuInput)
      Call StdIn_Name(StdIn)
      Call Molcas_open(LuInput,StdIn)
*
      iRow=0
      iRow_c=0
      nFix=0
      iInt=0
      Call RdCtl_Slapaf(LuInput,.True.)
      Curvilinear = .True.
      Cartesian = .Not. Curvilinear
      Numerical = .False.
*
      Close(LuInput)
      i=iPrintLevel(PLback)
*                                                                      *
************************************************************************
*                                                                      *
      BSet=.True.
      HSet=.False.
      PrQ=.False.
      nWndw=iter
      iRef=0
      Call BMtrx(SIZE(Coor,2),Coor,iter,mTtAtm,nWndw)
      nQQ = SIZE(Shift,1)
*                                                                      *
************************************************************************
*                                                                      *
*---- First get the (transposed) K matrix, (dQ/dq)^T
*
      Call mma_allocate(KMtrx,mq,nQQ,label="KMtrx")
      Call mma_allocate(KTrsp,nQQ,mq,label="KTrsp")
      Call Qpg_dArray('K',Found,nK)
      If (.not.Found .or. (nK.ne.mq*nQQ)) Call Abend()
      Call Get_dArray('K',KMtrx,nK)
      KTrsp(:,:)=Transpose(KMtrx)
      Call mma_deallocate(KMtrx)
*
*---- Form the full B matrix in the redundant internal coordinates
*     (primitives) and multiply by KK(t) on the fly
*
      nX = 3*mTtAtm
      Call mma_allocate(KKtB,mq,nX,label="KKtB")
      Call FZero(KKtB,nX*mq)
      i=1
      Do iq=1,mq
         nB=nqBM(iq)
         Do iB=i,i+nB-1
            j=iBM(iB)
            Do ii=1,mq
               KKtB(ii,j)=KKtB(ii,j)+BM(iB)*
     &            DDot_(nQQ,KTrsp(1,ii),1,KTrsp(1,iq),1)
            End Do
         End Do
         i=i+nB
      End Do
      Call mma_deallocate(KTrsp)
*
*---- Get the Cartesian normal modes
*
      Call Get_iScalar('Unique atoms',nUnique_Atoms)
      Call Get_nAtoms_All(nAll_Atoms)
      If (3*nAll_Atoms.ne.nX) Call Abend()
      Call mma_allocate(NMod,nX,nModes,label="NMod")
      Call FZero(NMod,nX*nModes)
      Call Get_NMode_All(Modes,lModes,nModes,nUnique_Atoms,
     &                   NMod,nAll_Atoms,lDisp)
*
*---- Compute the overlaps with the primitive displacements
*
      Call mma_allocate(IntMod,mq,nModes,label="IntMod")
      Call DGeMM_('N','N',mq,nModes,nX,One,KKtB,mq,
     &            NMod,nX,Zero,IntMod,mq)
      Call mma_deallocate(KKtB)
      Call mma_deallocate(NMod)
*
*---- "Normalize" the maximum value to 1
*
      Do i=1,nModes
         Mx=Zero
         Do j=1,mq
           Mx=Max(Mx,Abs(IntMod(j,i)))
         End Do
         If (Mx.gt.1.0D-10) Call DScal_(mq,One/Mx,IntMod(1,i),1)
      End Do
*
*---- Print the overlaps
*
      Call mma_allocate(Label,nLbl,label='Label')
      Filename='INTCOR'
      LuIC=21
      LuIC=IsFreeUnit(LuIC)
      Call Molcas_open(LuIC,Filename)
      i=1
      Line=Get_Ln_EOF(LuIC)
      Do While (Line .ne. 'EOF')
         j=Index(Line,'=')+1
         Label(i)=AdjustL(Line(j:))
         i=i+1
         Line=Get_Ln_EOF(LuIC)
      End Do
      Close(LuIC)
*
      MinComp=Half
      Call CollapseOutput(1,'Principal components of the normal modes')
      Write(6,'(3X,A)') '----------------------------------------'
      Write(6,*)
      Write(6,'(3X,A,F4.2,A)') '(Only contributions larger than ',
     &                         MinComp,' times the maximum are printed)'
      Write(6,*)
      Call mma_allocate(Sort,mq,label="Sort")
      Do i=1,nModes
         Write(6,*)
         Write(6,'(6X,A,1X,I6)') 'Mode',i
         Write(Line,'(5X,F10.2)') Freq(i)
         j=Index(Line,'-')
         If (j.gt.0) Line(j:j)='i'
         Write(6,'(8X,A)') 'Frequency: '//Trim(Line)//' cm-1'
         Write(6,'(6X,A)') '---------------------------------'
         Do j=1,mq
            Sort(j)=j
         End Do
         Do j=1,mq
            Do ii=j+1,mq
               If (Abs(IntMod(Sort(ii),i)).gt.
     &             Abs(IntMod(Sort(j),i))) Then
                  im=Sort(j)
                  Sort(j)=Sort(ii)
                  Sort(ii)=im
               End If
            End Do
            If (Abs(IntMod(Sort(j),i)).lt.MinComp) Exit
            Write(6,'(8X,A,F7.4)') Label(Sort(j)),IntMod(Sort(j),i)
         End Do
         Write(6,'(6X,A)') '---------------------------------'
      End Do
      Call CollapseOutput(0,'Principal components of the normal modes')
*
*---- Clean up
*
      Call mma_deallocate(Label)
      Call mma_deallocate(IntMod)
      Call mma_deallocate(Sort)
*                                                                      *
************************************************************************
*                                                                      *
*---- Restore the runfile and the "global" variables
*
      Call fCopy('RUNBCK2','RUNFILE',iErr)
      If (iErr.ne.0) Call Abend()
      If (AixRm('RUNBCK2').ne.0) Call Abend
*
      Header(:)=Bk_Header(:)
      iRef=Bk_iRef
      NmIter=Bk_NmIter
      iter=Bk_iter
      nDimBC=Bk_nDimBC
      MxItr=Bk_MxItr
      Max_Center=Bk_Max_Center
      iOptC=Bk_iOptC
      mode=Bk_mode
      mTROld=Bk_mTROld
      nWndw=Bk_nWndw
      iOptH=Bk_iOptH
      nLambda=Bk_nLambda
      nMEP=Bk_nMEP
      nBVec=Bk_nBVec
      IRC=Bk_IRC
      mTtAtm=Bk_mTtAtm
      MEPnum=Bk_MEPnum
      RootMap(:)=Bk_RootMap(:)
      Stop=Bk_Stop
      lOld=Bk_lOld
      CurviLinear=Bk_CurviLinear
      HSet=Bk_HSet
      BSet=Bk_BSet
      lNmHss=Bk_lNmHss
      Cubic=Bk_Cubic
      Baker=Bk_Baker
      DDV_Schlegel=Bk_DDV_Schlegel
      Line_Search=Bk_Line_Search
      HWRS=Bk_HWRS
      Analytic_Hessian=Bk_Analytic_Hessian
      FindTS=Bk_FindTS
      MEP=Bk_MEP
      User_Def=Bk_User_Def
      rMEP=Bk_rMEP
      lOld_Implicit=Bk_lOld_Implicit
      HrmFrq_Show=Bk_HrmFrq_Show
      eMEPtest=Bk_eMEPtest
      Redundant=Bk_Redundant
      lCtoF=Bk_lCtoF
      lSoft=Bk_lSoft
      CallLast=Bk_CallLast
      TwoRunFiles=Bk_TwoRunFiles
      TSConstraints=Bk_TSConstraints
      MEPCons=Bk_MEPCons
      Track=Bk_Track
      Request_Alaska=Bk_Request_Alaska
      Request_RASSI=Bk_Request_RASSI
      cMass(:)=Bk_cMass(:)
      ThrEne=Bk_ThrEne
      ThrGrd=Bk_ThrGrd
      Beta_Disp=Bk_Beta_Disp
      Delta=Bk_Delta
      Rtrnc=Bk_Rtrnc
      rHidden=Bk_rHidden
      ThrCons=Bk_ThrCons
      GNrm_Threshold=Bk_GNrm_Threshold
      CnstWght=Bk_CnstWght
      dMEPStep=Bk_dMEPStep
      rFuzz=Bk_rFuzz
      ThrMEP=Bk_ThrMEP
      lTherm=Bk_lTherm
      lDoubleIso=Bk_lDoubleIso
      nUserPT=Bk_nUserPT
      nsRot=Bk_nsRot
      UserT(:)=Bk_UserT(:)
      UserP=Bk_UserP
      HUpMet=Bk_HUpMet
      UpMeth=Bk_UpMeth
      MEP_Type=Bk_MEP_Type
      MEP_Algo=Bk_MEP_Algo
      isFalcon=Bk_isFalcon
      mB_Tot=Bk_mB_Tot
      mdB_Tot=Bk_mdB_Tot
      mq=Bk_mq
      iSBS=Bk_iSBS
      WeightedConstraints=Bk_WeightedConstraints
      NADC=Bk_NADC
      EDiffZero=Bk_EDiffZero
      ApproxNADC=Bk_ApproxNADC
      iState(:)=Bk_iState(:)
      iRow = Bk_iRow
      iRow_c=Bk_iRow_c
      nFix = Bk_nFix
      iInt = Bk_iInt
*
*     Process arrays that is always allocated.
*
      If (Allocated(Bk_Cx)) Then
         Cx(:,:,:) = Bk_Cx(:,:,:)
         Call mma_deallocate(Bk_Cx)
      Else
         Call mma_deallocate(Cx)
      End If
      If (Allocated(Bk_Gx)) Then
         Gx(:,:,:) = Bk_Gx(:,:,:)
         Call mma_deallocate(Bk_Gx)
      Else
         Call mma_deallocate(Gx)
      End If
      If (Allocated(Bk_Gx0)) Then
         Gx0(:,:,:) = Bk_Gx0(:,:,:)
         Call mma_deallocate(Bk_Gx0)
      Else
         Call mma_deallocate(Gx0)
      End If
      If (Allocated(Bk_Q_nuclear)) Then
         Q_nuclear(:) = Bk_Q_nuclear(:)
         Call mma_deallocate(Bk_Q_nuclear)
      Else
         Call mma_deallocate(Q_nuclear)
      End If
      If (Allocated(Bk_dMass)) Then
         dMass(:) = Bk_dMass(:)
         Call mma_deallocate(Bk_dMass)
      Else
         Call mma_deallocate(dMass)
      End If
      If (Allocated(Bk_Coor)) Then
         Coor(:,:) = Bk_Coor(:,:)
         Call mma_deallocate(Bk_Coor)
      Else
        If (Allocated(Coor)) Call mma_deallocate(Coor)
      End If
      If (Allocated(Bk_Grd)) Then
         Grd(:,:) = Bk_Grd(:,:)
         Call mma_deallocate(Bk_Grd)
      Else
        If (Allocated(Grd)) Call mma_deallocate(Grd)
      End If
      If (Allocated(Bk_ANr)) Then
         ANr(:) = Bk_ANr(:)
         Call mma_deallocate(Bk_ANr)
      Else
         Call mma_deallocate(ANr)
      End If
      If (Allocated(Bk_Weights)) Then
         Weights(:) = Bk_Weights(:)
         Call mma_deallocate(Bk_Weights)
      Else
         Call mma_deallocate(Weights)
      End If
      If (Allocated(Bk_Shift)) Then

         If (SIZE(Shift,1)/=SIZE(Bk_Shift,1)) Then
            Call mma_deallocate(Shift)
            Call mma_allocate(Shift,SIZE(Bk_Shift,1),MaxItr,
     &                        Label='Shift')
         End If
         Shift(:,:) = Bk_Shift(:,:)
         Call mma_deallocate(Bk_Shift)
      Else
         Call mma_deallocate(Shift)
      End If
      If (Allocated(Bk_GNrm)) Then
         GNrm(:) = Bk_GNrm(:)
         Call mma_deallocate(Bk_GNrm)
      Else
         Call mma_deallocate(GNrm)
      End If
      If (Allocated(Bk_Energy)) Then
         Energy(:) = Bk_Energy(:)
         Call mma_deallocate(Bk_Energy)
      Else
         Call mma_deallocate(Energy)
      End If
      If (Allocated(Bk_Energy0)) Then
         Energy0(:) = Bk_Energy0(:)
         Call mma_deallocate(Bk_Energy0)
      Else
         Call mma_deallocate(Energy0)
      End If
      If (Allocated(Bk_MF)) Then
         MF(:,:) = Bk_MF(:,:)
         Call mma_deallocate(Bk_MF)
      Else
         Call mma_deallocate(MF)
      End If
      If (Allocated(Bk_DipM)) Then
         DipM(:,:) = Bk_DipM(:,:)
         Call mma_deallocate(Bk_DipM)
      Else
        If (Allocated(DipM)) Call mma_deallocate(DipM)
      End If
      If (Allocated(Bk_RefGeo)) Then
         RefGeo(:,:) = Bk_RefGeo(:,:)
         Call mma_deallocate(Bk_RefGeo)
      Else
        If (Allocated(RefGeo)) Call mma_deallocate(RefGeo)
      End If
*
*     Process arrays that is allocated optionally.
*
      If (Allocated(Bk_NAC)) Then
         NAC(:,:) = Bk_NAC(:,:)
         Call mma_deallocate(Bk_NAC)
      Else
        If (Allocated(NAC)) Call mma_deallocate(NAC)
      End If
      If (Allocated(Bk_Lambda)) Then
         If (.NOT.Allocated(Lambda)) Then
            Call mma_allocate(Lambda,SIZE(Bk_Lambda,1),MaxItr+1,
     &                        Label='Lambda')
         End If
         Lambda(:,:) = Bk_Lambda(:,:)
         Call mma_deallocate(Bk_Lambda)
      Else
        If (Allocated(Lambda)) Call mma_deallocate(Lambda)
      End If
      If (Allocated(Bk_qInt)) Then
         If (.NOT.Allocated(qInt)) Then
            Call mma_allocate(qInt,SIZE(Bk_qInt,1),MaxItr,
     &                        Label='qInt')
         End If
         qInt(:,:) = Bk_qInt(:,:)
         Call mma_deallocate(Bk_qInt)
      Else
        If (Allocated(qInt)) Call mma_deallocate(qInt)
      End If
      If (Allocated(Bk_dqInt)) Then
         If (.NOT.Allocated(dqInt)) Then
            Call mma_allocate(dqInt,SIZE(Bk_dqInt,1),MaxItr,
     &                        Label='dqInt')
         End If
         dqInt(:,:) = Bk_dqInt(:,:)
         Call mma_deallocate(Bk_dqInt)
      Else
        If (Allocated(dqInt)) Call mma_deallocate(dqInt)
      End If

      If (Allocated(Bk_BM)) Then
         If (Allocated(BM)) Call mma_deallocate(BM)
         Call mma_allocate(BM,SIZE(Bk_BM),Label='BM')
         BM(:) = Bk_BM(:)
         Call mma_deallocate(Bk_BM)
      Else
         If (Allocated(BM)) Call mma_deallocate(BM)
      End If
      If (Allocated(Bk_dBM)) Then
         If (Allocated(dBM)) Call mma_deallocate(dBM)
         Call mma_allocate(dBM,SIZE(Bk_dBM),Label='dBM')
         dBM(:) = Bk_dBM(:)
         Call mma_deallocate(Bk_dBM)
      Else
         If (Allocated(dBM)) Call mma_deallocate(dBM)
      End If
      If (Allocated(Bk_iBM)) Then
         If (Allocated(iBM)) Call mma_deallocate(iBM)
         Call mma_allocate(iBM,SIZE(Bk_iBM),Label='iBM')
         iBM(:) = Bk_iBM(:)
         Call mma_deallocate(Bk_iBM)
      Else
         If (Allocated(iBM)) Call mma_deallocate(iBM)
      End If
      If (Allocated(Bk_idBM)) Then
         If (Allocated(idBM)) Call mma_deallocate(idBM)
         Call mma_allocate(idBM,SIZE(Bk_idBM),Label='idBM')
         idBM(:) = Bk_idBM(:)
         Call mma_deallocate(Bk_idBM)
      Else
         If (Allocated(idBM)) Call mma_deallocate(idBM)
      End If
      If (Allocated(Bk_nqBM)) Then
         If (Allocated(nqBM)) Call mma_deallocate(nqBM)
         Call mma_allocate(nqBM,SIZE(Bk_nqBM),Label='nqBM')
         nqBM(:) = Bk_nqBM(:)
         Call mma_deallocate(Bk_nqBM)
      Else
         If (Allocated(nqBM)) Call mma_deallocate(nqBM)
      End If
      If (Allocated(Bk_BMx)) Then
         If (Allocated(BMx)) Call mma_deallocate(BMx)
         Call mma_allocate(BMx,SIZE(Bk_BMx,1),SIZE(Bk_BMx,2),
     &                     Label='BMx')
         BMx(:,:) = Bk_BMx(:,:)
         Call mma_deallocate(Bk_BMx)
      Else
         If (Allocated(BMx)) Call mma_deallocate(BMx)
      End If
      If (Allocated(Bk_Degen)) Then
         If (Allocated(Degen)) Call mma_deallocate(Degen)
         Call mma_allocate(Degen,SIZE(Bk_Degen,1),SIZE(Bk_Degen,2),
     &                     Label='Degen')
         Degen(:,:) = Bk_Degen(:,:)
         Call mma_deallocate(Bk_Degen)
      Else
         If (Allocated(Degen)) Call mma_deallocate(Degen)
      End If
      If (Allocated(Bk_jStab)) Then
         If (Allocated(jStab)) Call mma_deallocate(jStab)
         Call mma_allocate(jStab,[0,7],[1,SIZE(Bk_jStab,2)],
     &                     Label='jStab')
         jStab(:,:) = Bk_jStab(:,:)
         Call mma_deallocate(Bk_jStab)
      Else
         If (Allocated(jStab)) Call mma_deallocate(jStab)
      End If
      If (Allocated(Bk_iCoSet)) Then
         If (Allocated(iCoSet)) Call mma_deallocate(iCoSet)
         Call mma_allocate(iCoSet,[0,7],[1,SIZE(Bk_iCoSet,2)],
     &                     Label='iCoSet')
         iCoSet(:,:) = Bk_iCoSet(:,:)
         Call mma_deallocate(Bk_iCoSet)
      Else
         If (Allocated(iCoSet)) Call mma_deallocate(iCoSet)
      End If
      If (Allocated(Bk_nStab)) Then
         If (Allocated(nStab)) Call mma_deallocate(nStab)
         Call mma_allocate(nStab,SIZE(Bk_nStab,1),Label='nStab')
         nStab(:) = Bk_nStab(:)
         Call mma_deallocate(Bk_nStab)
      Else
         If (Allocated(nStab)) Call mma_deallocate(nStab)
      End If
      If (Allocated(Bk_AtomLbl)) Then
         If (Allocated(AtomLbl)) Call mma_deallocate(AtomLbl)
         Call mma_allocate(AtomLbl,SIZE(Bk_AtomLbl,1),Label='AtomLbl')
         AtomLbl(:) = Bk_AtomLbl(:)
         Call mma_deallocate(Bk_AtomLbl)
      Else
         If (Allocated(AtomLbl)) Call mma_deallocate(AtomLbl)
      End If
      If (Allocated(Bk_Smmtrc)) Then
         If (Allocated(Smmtrc)) Call mma_deallocate(Smmtrc)
         Call mma_allocate(Smmtrc,3,SIZE(Bk_Smmtrc,2),Label='Smmtrc')
         Smmtrc(:,:) = Bk_Smmtrc(:,:)
         Call mma_deallocate(Bk_Smmtrc)
      Else
         If (Allocated(Smmtrc)) Call mma_deallocate(Smmtrc)
      End If
      If (Allocated(Bk_Lbl)) Then
         If (Allocated(Lbl)) Call mma_deallocate(Lbl)
         Call mma_allocate(Lbl,SIZE(Bk_Lbl),Label='Lbl')
         Lbl(:) = Bk_Lbl(:)
         Call mma_deallocate(Bk_Lbl)
      Else
         If (Allocated(Lbl)) Call mma_deallocate(Lbl)
      End If
      If (Allocated(Bk_mRowH)) Then
         If (Allocated(mRowH)) Call mma_deallocate(mRowH)
         Call mma_allocate(mRowH,SIZE(Bk_mRowH),Label='mRowH')
         mRowH(:) = Bk_mRowH(:)
         Call mma_deallocate(Bk_mRowH)
      Else
         If (Allocated(mRowH)) Call mma_deallocate(mRowH)
      End If
*
      End Subroutine Print_Mode_Components
