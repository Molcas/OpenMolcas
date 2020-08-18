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
      Implicit None
#include "backup_info.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "real.fh"
      Real*8 :: Modes(*), Freq(*), Mx, MinComp
      Integer :: LuInput, iRow, iInt, nFix, iRout,
     &           iPrint, nX, i, j, nB, iq, nAll_Atoms, nUnique_Atoms,
     &           iB, lDisp(nSym), nModes, lModes, LuIC, ii, im, nK,
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
*
      Call QEnter('Print_Mode_Components')
*
*---- Ugly hack: backup all "global" slapaf variables in case this is
*                called from inside slapaf
*
      Bk_iOper(0:7)=iOper(0:7)
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
      Bk_ipdqInt=ipdqInt
      Bk_ipH=ipH
      Bk_ipQInt=ipQInt
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
      Bk_ipCx=ipCx
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
      Bk_nqInt=nqInt
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
      Bk_nSym=nSym
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
      Bk_Stat(0:MaxItr)=Stat(0:MaxItr)
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
      Bk_ipNADC=ipNADC
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
      Call Put_dArray('Hess',Work(ip_Dummy),0)
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
      Call BMtrx(iRow,nBVec,ipB,nsAtom,mInt,ipqInt,Lbl,
     &           Work(ipCoor),nDimBC,Work(ipCM),AtomLbl,nSym,iOper,
     &           Smmtrc,Degen,BSet,HSet,iter,ipdqInt,ipShf,
     &           Work(ipGx),Work(ipCx),mTtAtm,iWork(ipANr),iOptH,
     &           User_Def,nStab,jStab,Curvilinear,Numerical,
     &           DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &           iCoSet,lOld,rHidden,nFix,nQQ,iRef,Redundant,nqInt,
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
      If (ipNADC.ne.ip_Dummy) Call Free_Work(ipNADC)
      If (ipqInt.ne.ip_Dummy) Then
         Call GetMem('dqInt', 'Free','Real',ipdqInt, nqInt)
         Call GetMem('qInt',  'Free','Real',ipqInt,  nqInt)
      End If
      If (BSet) Then
         Call GetMem('Shift','Free','Real',ipShf,nQQ*iter)
      End If
      If (Ref_Geom) Then
         Call GetMem('ipRef',  'Free','Real',ipRef,    3*nsAtom)
      End If
      Call GetMem('Relax',  'Free','Real',ipRlx,    Lngth)
      Call GetMem('Grad',   'Free','Real',ipGrd,    3*nsAtom)
      Call GetMem('Coord',  'Free','Real',ipCoor,   3*nsAtom)
      Call GetMem('Anr',    'Free','Inte',ipANr,    nsAtom)
      Call GetMem('Charge', 'Free','Real',ipCM,     nsAtom)
      Call GetMem('Weights','Free','Real',ipWeights,nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*---- Restore the runfile and the "global" variables
*
      Call fCopy('RUNBCK2','RUNFILE',iErr)
      If (iErr.ne.0) Call Abend()
      If (AixRm('RUNBCK2').ne.0) Call Abend
*
      iOper(0:7)=Bk_iOper(0:7)
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
      ipdqInt=Bk_ipdqInt
      ipH=Bk_ipH
      ipQInt=Bk_ipQInt
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
      nqInt=Bk_nqInt
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
      nSym=Bk_nSym
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
      Stat(0:MaxItr)=Bk_Stat(0:MaxItr)
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
      ipNADC=Bk_ipNADC
      iState(:)=Bk_iState(:)
*
      Call QExit('Print_Mode_Components')
*
      End Subroutine
