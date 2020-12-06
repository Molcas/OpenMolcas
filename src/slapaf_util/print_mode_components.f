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
     &                       RefGeo
      Implicit None
#include "backup_info.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
      Real*8 rDum(1)
      Real*8 :: Modes(*), Freq(*), Mx, MinComp
      Integer :: LuInput, iRow, iInt, nFix, iRout,
     &           iPrint, nX, i, j, nB, iq, nAll_Atoms, nUnique_Atoms,
     &           iB, lDisp(nIrrep), nModes, lModes, LuIC, ii, im, nK,
     &           iErr, PLback
      Real*8, Dimension(:,:), Allocatable :: KMtrx, KTrsp, KKtB, IntMod,
     &                                       NMod
      Integer, Dimension(:), Allocatable :: Sort
      Integer, Parameter :: nLbl=10*MxAtom
      Integer, External :: IsFreeUnit, iPrintLevel, AixRm
      Logical :: Cartesian, Numerical, PrQ, Found
      Character(Len=8) :: Lbl(nLbl),Filename
      Character(Len=16) :: StdIn
      Character(Len=24), Allocatable :: Label(:)
      Character(Len=180), External :: Get_Ln_EOF
      Real*8, External :: DDot_
      Integer, External :: ip_of_Work

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

      Integer, Allocatable:: Bk_ANr(:)
*
*
*---- Ugly hack: backup all "global" slapaf variables in case this is
*                called from inside slapaf
*
*     Note, this routine might be called from outside the Slapaf
*     environment. In which case there is no backup to be made.
*
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

      Bk_iSym(:)=iSym(:)
      Bk_iCoSet(0:7,:)=iCoSet(0:7,:)
      Bk_nStab(:)=nStab(:)
      Bk_ipAtom=ipAtom
      Bk_ipNSup=ipNSup
      Bk_ipR12=ipR12
      Bk_iRef=iRef
      Bk_nQQ=nQQ
      Bk_mRowH(:)=mRowH(:)
      Bk_nRowH=nRowH
      Bk_ipEner=ipEner
      Bk_ipGnrm=ipGnrm
      Bk_ipH=ipH
      Bk_NmIter=NmIter
      Bk_ipShf=ipShf
      Bk_ipGrd=ipGrd
      Bk_ipRlx=ipRlx
      Bk_ipStat=ipStat
      Bk_iter=iter
      Bk_Lngth=Lngth
      Bk_nDimBC=nDimBC
      Bk_MxItr=MxItr
      Bk_mInt=mInt
      Bk_ipB=ipB
      Bk_ipBt=ipBt
      Bk_ipBVec=ipBVec
      Bk_Max_Center=Max_Center
      Bk_ipQ0=ipQ0
      Bk_ipVal=ipVal
      Bk_nSupSy=nSupSy
      Bk_ipBOld=ipBOld
      Bk_lif=lif
      Bk_ipCx= ipCx
      Bk_ipGx=ipGx
      Bk_ipANr=ipANr
      Bk_iOptC=iOptC
      Bk_mode=mode
      Bk_mTROld=mTROld
      Bk_nWndw=nWndw
      Bk_iOptH=iOptH
      Bk_jStab(0:7,:)=jStab(0:7,:)
      Bk_ipMF=ipMF
      Bk_nLambda=nLambda
      Bk_iRow_c=iRow_c
      Bk_ipL=ipL
      Bk_ipEner0=ipEner0
      Bk_ipGx0=ipGx0
      Bk_nMEP=nMEP
      Bk_ipRef=ipRef
      Bk_ipGradRef=ipGradRef
      Bk_ipDipM=ipDipM
      Bk_ipK_Ref=ipK_Ref
      Bk_nBVec=nBVec
      Bk_IRC=IRC
      Bk_ipCM=ipCM
      Bk_ipCoor=ipCoor
      Bk_mTtAtm=mTtAtm
      Bk_nsAtom=nsAtom
      Bk_MEPnum=MEPnum
      Bk_RootMap(:)=RootMap(:)
      Bk_Smmtrc(:)=Smmtrc(:)
      Bk_Stop=Stop
      Bk_lWrite=lWrite
      Bk_Exist=Exist
      Bk_Change=Change
      Bk_lSup=lSup
      Bk_lOld=lOld
      Bk_CurviLinear=CurviLinear
      Bk_lRowH=lRowH
      Bk_HSet=HSet
      Bk_BSet=BSet
      Bk_lNmHss=lNmHss
      Bk_Cubic=Cubic
      Bk_PDH=PDH
      Bk_Baker=Baker
      Bk_Schlegel=Schlegel
      Bk_DDV_Schlegel=DDV_Schlegel
      Bk_Line_Search=Line_Search
      Bk_HWRS=HWRS
      Bk_Analytic_Hessian=Analytic_Hessian
      Bk_FirstCall=FirstCall
      Bk_FindTS=FindTS
      Bk_MEP=MEP
      Bk_Ref_Geom=Ref_Geom
      Bk_lRP=lRP
      Bk_User_Def=User_Def
      Bk_Ref_Grad=Ref_Grad
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
      Bk_Degen(:)=Degen(:)
      Bk_cMass(:)=cMass(:)
      Bk_Trial(:)=Trial(:)
      Bk_ThrEne=ThrEne
      Bk_ThrGrd=ThrGrd
      Bk_Beta=Beta
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
      Bk_Char=Char
      Bk_Header(:)=Header(:)
      Bk_Line=Line
      Bk_BLine=BLine
      Bk_HUpMet=HUpMet
      Bk_UpMeth=UpMeth
      Bk_AtomLbl(:)=AtomLbl(:)
      Bk_MEP_Type=MEP_Type
      Bk_MEP_Algo=MEP_Algo
      Bk_isFalcon=isFalcon
      Bk_ip_B=ip_B
      Bk_ip_dB=ip_dB
      Bk_ip_iB=ip_iB
      Bk_ip_idB=ip_idB
      Bk_ip_nqB=ip_nqB
      Bk_mB_Tot=mB_Tot
      Bk_mdB_Tot=mdB_Tot
      Bk_mq=mq
      Bk_iSBS=iSBS
      Bk_ipWeights=ipWeights
      Bk_WeightedConstraints=WeightedConstraints
      Bk_NADC=NADC
      Bk_EDiffZero=EDiffZero
      Bk_ApproxNADC=ApproxNADC
      Bk_iState(:)=iState(:)
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
      Call RdCtl_Slapaf(iRow,iInt,nFix,LuInput,.True.)
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
      nFix=0
      nWndw=iter
      iRef=0
      Call BMtrx(iRow,nBVec,ipB,nsAtom,mInt,Lbl,
     &           Coor,nDimBC,AtomLbl,Smmtrc,Degen,BSet,HSet,iter,
     &           mTtAtm,iOptH,
     &           User_Def,nStab,jStab,Curvilinear,Numerical,
     &           DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &           iCoSet,lOld,rHidden,nFix,nQQ,iRef,Redundant,
     &           MaxItr,nWndw)
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
      i=0
      Do iq=1,mq
         nB=iWork(ip_nqB+iq-1)
         Do iB=i,i+nB-1
            j=iWork(ip_iB+iB)
            Do ii=1,mq
               KKtB(ii,j)=KKtB(ii,j)+Work(ip_B+iB)*
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
      Call GetMem('BM','Free','Real',ip_B,mB_Tot)
      Call GetMem('iBM','Free','Inte',ip_iB,mB_Tot)
      Call GetMem('nqB','Free','Inte',ip_nqB,mq)
      Call GetMem(' B ',    'Free','Real',ipB,   (nsAtom*3)**2)
      Call GetMem('Relax',  'Free','Real',ipRlx,    Lngth)
*                                                                      *
************************************************************************
*                                                                      *
*---- Restore the runfile and the "global" variables
*
      Call fCopy('RUNBCK2','RUNFILE',iErr)
      If (iErr.ne.0) Call Abend()
      If (AixRm('RUNBCK2').ne.0) Call Abend
*
      iSym(:)=Bk_iSym(:)
      iCoSet(0:7,:)=Bk_iCoSet(0:7,:)
      nStab(:)=Bk_nStab(:)
      ipAtom=Bk_ipAtom
      ipNSup=Bk_ipNSup
      ipR12=Bk_ipR12
      iRef=Bk_iRef
      nQQ=Bk_nQQ
      mRowH(:)=Bk_mRowH(:)
      nRowH=Bk_nRowH
      ipEner=Bk_ipEner
      ipGnrm=Bk_ipGnrm
      ipH=Bk_ipH
      NmIter=Bk_NmIter
      ipShf=Bk_ipShf
      ipGrd=Bk_ipGrd
      ipRlx=Bk_ipRlx
      ipStat=Bk_ipStat
      iter=Bk_iter
      Lngth=Bk_Lngth
      nDimBC=Bk_nDimBC
      MxItr=Bk_MxItr
      mInt=Bk_mInt
      ipB=Bk_ipB
      ipBt=Bk_ipBt
      ipBVec=Bk_ipBVec
      Max_Center=Bk_Max_Center
      ipQ0=Bk_ipQ0
      ipVal=Bk_ipVal
      nSupSy=Bk_nSupSy
      ipBOld=Bk_ipBOld
      lif=Bk_lif
      ipCx=Bk_ipCx
      ipGx=Bk_ipGx
      ipANr=Bk_ipANr
      iOptC=Bk_iOptC
      mode=Bk_mode
      mTROld=Bk_mTROld
      nWndw=Bk_nWndw
      iOptH=Bk_iOptH
      jStab(0:7,:)=Bk_jStab(0:7,:)
      ipMF=Bk_ipMF
      nLambda=Bk_nLambda
      iRow_c=Bk_iRow_c
      ipL=Bk_ipL
      ipEner0=Bk_ipEner0
      ipGx0=Bk_ipGx0
      nMEP=Bk_nMEP
      ipRef=Bk_ipRef
      ipGradRef=Bk_ipGradRef
      ipDipM=Bk_ipDipM
      ipK_Ref=Bk_ipK_Ref
      nBVec=Bk_nBVec
      IRC=Bk_IRC
      ipCM=Bk_ipCM
      ipCoor=Bk_ipCoor
      mTtAtm=Bk_mTtAtm
      nsAtom=Bk_nsAtom
      MEPnum=Bk_MEPnum
      RootMap(:)=Bk_RootMap(:)
      Smmtrc(:)=Bk_Smmtrc(:)
      Stop=Bk_Stop
      lWrite=Bk_lWrite
      Exist=Bk_Exist
      Change=Bk_Change
      lSup=Bk_lSup
      lOld=Bk_lOld
      CurviLinear=Bk_CurviLinear
      lRowH=Bk_lRowH
      HSet=Bk_HSet
      BSet=Bk_BSet
      lNmHss=Bk_lNmHss
      Cubic=Bk_Cubic
      PDH=Bk_PDH
      Baker=Bk_Baker
      Schlegel=Bk_Schlegel
      DDV_Schlegel=Bk_DDV_Schlegel
      Line_Search=Bk_Line_Search
      HWRS=Bk_HWRS
      Analytic_Hessian=Bk_Analytic_Hessian
      FirstCall=Bk_FirstCall
      FindTS=Bk_FindTS
      MEP=Bk_MEP
      Ref_Geom=Bk_Ref_Geom
      lRP=Bk_lRP
      User_Def=Bk_User_Def
      Ref_Grad=Bk_Ref_Grad
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
      Degen(:)=Bk_Degen(:)
      cMass(:)=Bk_cMass(:)
      Trial(:)=Bk_Trial(:)
      ThrEne=Bk_ThrEne
      ThrGrd=Bk_ThrGrd
      Beta=Bk_Beta
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
      Char=Bk_Char
      Header(:)=Bk_Header(:)
      Line=Bk_Line
      BLine=Bk_BLine
      HUpMet=Bk_HUpMet
      UpMeth=Bk_UpMeth
      AtomLbl(:)=Bk_AtomLbl(:)
      MEP_Type=Bk_MEP_Type
      MEP_Algo=Bk_MEP_Algo
      isFalcon=Bk_isFalcon
      ip_B=Bk_ip_B
      ip_dB=Bk_ip_dB
      ip_iB=Bk_ip_iB
      ip_idB=Bk_ip_idB
      ip_nqB=Bk_ip_nqB
      mB_Tot=Bk_mB_Tot
      mdB_Tot=Bk_mdB_Tot
      mq=Bk_mq
      iSBS=Bk_iSBS
      ipWeights=Bk_ipWeights
      WeightedConstraints=Bk_WeightedConstraints
      NADC=Bk_NADC
      EDiffZero=Bk_EDiffZero
      ApproxNADC=Bk_ApproxNADC
      iState(:)=Bk_iState(:)
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
*
      End Subroutine
