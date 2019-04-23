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
      Subroutine Init_SlapAf(iRow)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "sbs.fh"
#include "nadc.fh"
#include "db.fh"
#include "print.fh"
      Integer   iAdd(0:7) , jPrmt(0:7)
      Logical Same, Do_ESPF, Exist_2, Found, Reduce_Prt
      External Reduce_Prt
      Character *8 CMAX
      Integer Columbus
#include "SysDef.fh"
      Character*100 Get_SuperName, SuperName
      External Get_SuperName
      Data jPrmt/1,-1,-1,1,-1,1,1,-1/
*
*     Statement function
*
      iPrmt(i,j) = jPrmt(iAnd(i,j))
*
      Call QEnter('Init')
************************************************************************
************************** StartUp section   ***************************
************************************************************************
*                                                                      *
      Call Get_iScalar('System BitSwitch',iSBS)
*                                                                      *
************************************************************************
*                                                                      *
*     Set the default value of iterations from MOLCAS_MAXITER if it
*     has been defined.
*
      Call GetEnvf('MOLCAS_MAXITER', CMAX)
*     Write (*,'(3A)') 'CMAX="',CMAX,'"'
      If (CMAX.ne.' ') Then
         Read (CMAX,'(I8)') iMAX
         MxItr = Min(MaxItr,iMax)
      Else
         MxItr = MaxItr
      End If
*                                                                      *
************************************************************************
*                                                                      *
      IRC=0
      iRow=0
      iRow_c=0
      nBVec=0
      lif = 0
      jPrint=10
      lOld = .False.
      lOld_Implicit = .False.
      Stop  = .False.
      lSup = .False.
      Baker = .False.
      Schlegel=.False.
      FindTS=.False.
      DDV_Schlegel=.False.
      Curvilinear=.True.
      Ref_Geom=.False.
      HWRS=.True.
      nLambda=0
      MEP = .False.
      rMEP= .False.
      nMEP=MaxItr
      Ref_Grad=.False.
      rHidden = Zero
      HrmFrq_Show=.False.
      eMEPTest=.True.
      MEP_Type='SPHERE'
      dMEPStep=0.1D0
      MEP_Algo='GS'
      MEPnum=0
      NmIter=0
      HSet=.False.
      BSet=.False.
      Redundant=.False.
      lSoft=.False.
      rFuzz=0.5D0
      isFalcon=.False.
      CallLast=.True.
      TwoRunFiles=.False.
      TSConstraints=.False.
      MEPCons=.False.
      Track=.False.
      Request_Alaska=.False.
      Request_RASSI=.False.
*
      Call Kriging_Init()
*                                                                      *
************************************************************************
*                                                                      *
*.... Optimization method. DO NOT EVER GO BEYOND BIT 30!!!
*
*      iOptC=000000000 (  0) No optimization
*   0  iOptC=000000001 (  1) Quasi Newton-Raphson
*   1  iOptC=000000010 (  2) c1-DIIS
*   2  iOptC=000000100 (  4) c2-DIIS
*   3  iOptC=000001000 (  8) RS-RFO
*   4  iOptC=00001.... ( 16) DIIS, <dx|dx>
*   5  iOptC=00010.... ( 32) DIIS, <dx|g>
*   6  iOptC=00100.... ( 64) DIIS, <g|g>
*   7  iOptC=01....... (128) Minimum, if not set TS search
*   8  iOptC=10....... (256) Optimization with constraint
*   9  iOptC           (512) set: RS-I-RFO, unset: RS-P-RFO
*  10  iOptC          (1024) HMF augmented with weak interactions
*  11  iOptC          (2048) augmented HMF used for selection of
*                           internal coordinates
*  12  iOptC          (4096) set if FindTS
*  13  iOptC          (8192) set if FindTS and in TS regime
*
      UpMeth='  RF  '
      iOptC=8
      iOptC=iOr(iOptC,64 )
      iOptC=iOr(iOptC,128)
      iOptC=iOr(iOptC,512)
      iOptC=iOr(iOptC,1024)
      iOptC=iOr(iOptC,2048)
*                                                                      *
************************************************************************
*                                                                      *
*     Hessian update
* 1   iOptH=00000001 (  1) Meyer (disabled)
* 2   iOptH=00000010 (  2) BP (disabled)
* 3   iOptH=00000100 (  4) BFGS
* 4   iOptH=00001000 (  8) None
* 5   iOptH=00010000 ( 16) MPS, for TS search
* 6   iOptH=-.1..... ( 32) Not used
* 7   iOptH=01000000 ( 64) EU, for TS search
* 8   iOptH=10000000 (128) TS-BFGS, for TS search
*
      iOptH=4
      HUpMet=' None '
*                                                                      *
************************************************************************
*                                                                      *
      RtRnc=Three
      Max_Center=15
      mode=-99999
      ipAtom = -1
      ipNSup = -1
      lRowH  = .False.
      lNmHss = .False.
* --- ThermoChemistry for Numerical Hessian
      lTherm = .False.
      lDoubleIso = .False.
      nUserPT= 0
      UserP = 1.0d0
      Do i=1, 64
        UserT(i) = 0.0d0
      EndDo
      nsRot = 0
*
      Cubic  = .False.
      PDH    = .True.
      Beta = 0.30D0
      GNrm_Threshold=0.2D0
      CnstWght=1.0D0
      Call DecideOnESPF(Do_ESPF)
      If (Do_ESPF) Then
         ThrGrd = 0.003D0
         ThrEne = 1.0D-5
         Line_Search=.False.
      Else
         ThrGrd = 0.0003D0
         ThrEne = 1.0D-6
         Line_Search=.True.
      End If
      ThrCons = 1.0D10
      Delta  = 1.0D-2
      nWndw = 5
      Call ICopy(mxdc,[0],0,nStab,1)
      Call ICopy(8*mxdc,[0],0,iCoSet,1)
      Call ICopy(8*mxdc,[0],0,jStab,1)
      call dcopy_(3*mxdc,[Zero],0,Degen,1)
      do 180 i = 1, 3*mxdc
         Smmtrc(i) = .False.
 180  Continue
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (iPL.eq.2) Then
         iPL=5
      Else If (iPL.eq.3) Then
         iPL=6
      Else If (iPL.eq.4) Then
         iPL=99
      Else If (iPL.eq.5) Then
         iPL=99
      End If
      Do iRout = 1, nRout
         nPrint(iRout) = iPL
      End Do
*
*     Reduced print level of Slapaf parameters after the first iteration
*
      If (Reduce_Prt().and.iPL.le.5) Then
         Do iRout = 1, nRout
            nPrint(iRout) = iPL-1
         End Do
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
      BLine = ' '
*                                                                      *
************************************************************************
*                                                                      *
*     Get Molecular data
*
*...  Read the title
*
      Call Get_cArray('Seward Title',Header,144)
*
*...  Read the number of irreps
*
      Call Get_iScalar('nSym',nSym)
      nIrrep=nSym
*
*...  Read number of atoms, charges, coordinates, gradients and
*     atom labels
*
      Call Get_Molecule(ipCM,ipCoor,ipGrd,AtomLbl,nsAtom,mxdc)
*                                                                      *
************************************************************************
*                                                                      *
*...  Read the symmetry operators
*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*                                                                      *
************************************************************************
*                                                                      *
      ipNADC=ip_Dummy
      NADC=.False.
      ApproxNADC=.False.
      Call Get_iScalar('Columbus',Columbus)
      If (Columbus.eq.1) Then
*
*        C&M mode
*
         Call Get_iScalar('ColGradMode',iMode)
         If (iMode.eq.3) NADC=.True.
      Else
*
*        M mode
*
*        ISPIN should only be found for RASSCF-based
*        methods, so no CI mode for SCF, MP2, etc. (or that's the idea)
*
C        Write (6,*) 'See if CI'
         Call Qpg_iScalar('ISPIN',Found)
         If (Found) Then
            Call Get_iScalar('ISPIN',ISPIN1)
            Call Get_iScalar('LSYM',LSYM1)
         Else
            ISPIN1=0
            LSYM1=0
         End If
C        Write (6,*) 'iSpin=',ISPIN1
C        Write (6,*) 'lSym=',LSYM1
*
         Call f_Inquire('RUNFILE2',Exist_2)
C        Write (6,*) 'Exist_2=',Exist_2
         If (Exist_2) Then
            Call NameRun('RUNFILE2')
            Call Qpg_iScalar('ISPIN',Found)
            If (Found) Then
               Call Get_iScalar('ISPIN',ISPIN2)
               Call Get_iScalar('LSYM',LSYM2)
            Else
               ISPIN2=0
               LSYM2=0
            End If
            Call NameRun('RUNFILE')
         Else
            ISPIN2 = ISPIN1
            LSYM2 = LSYM1
         End If
C        Write (6,*) 'iSpin=',ISPIN1,ISPIN2
C        Write (6,*) 'lSym=',LSYM1,LSYM2
*
*
*        Do not add the constraint at the NumGrad stage
*
         SuperName=Get_Supername()
         If (SuperName.ne.'numerical_gradient') Then
            If ((ISPIN1.ne.0).and.(LSYM1.ne.0))
     &         NADC= (ISPIN1.eq.ISPIN2) .and. (LSYM1.eq.LSYM2)
C           NADC= .False. ! for debugging
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*...  Read or initialize the root map
*
      Call Qpg_iArray('Root Mapping',Found,nRM)
      If (nRM.gt.0) Then
         Call Get_iArray('Root Mapping',RootMap,nRM)
      Else
         Call Qpg_iScalar('Number of roots',Found)
         nRoots = 1
         If (Found) Call Get_iScalar('Number of roots',nRoots)
         Call iCopy(MxRoot,[0],0,RootMap,1)
         Do i=1,nRoots
            RootMap(i)=i
         End Do
      End If
*
*...  Check if there is an analytic Hessian
      Call Get_AnalHess(ipHess,nHess)
      Analytic_Hessian=nHess.ne.0
      If (Analytic_Hessian) Call Free_Work(ipHess)
      If (.Not.Analytic_Hessian) Then
         Call NameRun('RUNOLD')
         Call Get_AnalHess(ipHess,nHess)
         Analytic_Hessian=nHess.ne.0
         If (Analytic_Hessian) Call Free_Work(ipHess)
         Call NameRun('#Pop')
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---  Set up of symmetry and degeneracy
*
      Call ICopy(3,[0],0,iSym,1)
      Do 600 iIrrep = 1, Min(4,nSym-1)
         jIrrep = iIrrep
         If (iIrrep.eq.3) jIrrep = 4
         Do 601 k = 1, 3
            If (iAnd(iOper(jIrrep),2**(k-1)).ne.0) iSym(k) = 2**(k-1)
 601     Continue
 600  Continue
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute the number of total symmetric displacements
*
      nDimbc = 0
*...  Loop over the unique atoms
      Do 610 isAtom = 1, nsAtom
*...     Find character of center
         iAdr = ipCoor -1 + (isAtom-1)*3
         iChxyz=0
         Do i = 1, 3
            If (Work(iAdr+i).ne.Zero) Then
               Do iIrrep= 0, nSym-1
                  If (iAnd(2**(i-1),iOper(iIrrep)).ne.0)
     &               iChxyz=iOr(iChxyz,2**(i-1))
               End Do
            End If
         End Do
         nStb = 0
         Do iIrrep = 0, nSym-1
            If (iAnd(iChxyz,iOper(iIrrep)).eq.0) Then
               jStab(nStb,isAtom)=iOper(iIrrep)
               nStb = nStb + 1
            End If
         End Do
         nStab(isAtom)=nStb
*...     Find the coset representatives
         iCoSet(0,nsAtom) = 0      ! Put in the unit operator
         nCoSet = 1
         Do iIrrep = 1, nSym-1
            itest=iAnd(iChxyz,iOper(iIrrep))
            Same=.False.
            Do jCoSet = 0, nCoSet-1
               jTest = iAnd(iChxyz,iCoSet(jCoSet,isAtom))
               Same = jTest.eq.iTest
               If (Same) Go To 7777
            End Do
 7777       Continue
            If (.Not.Same) Then
               nCoSet = nCoSet + 1
               iCoSet(nCoSet-1,isAtom) = iOper(iIrrep)
            End If
         End Do
         If (nSym/nStb.ne.nCoSet) Then
            Call WarningMessage(2,' Error while doing cosets.')
            Call Abend()
         End If
         Do 611 i = 1, 3
            iComp = 2**(i-1)
            Call ICopy(nCoSet,[0],0,iAdd,1)
            Do 640 iIrrep = 0, nSym-1
*...           find the stabilizer index
               iTest=iAnd(iChxyz,iOper(iIrrep))
               n=-1
               Do 641 jCoset = 0, nCoset-1
                  jTest=iAnd(iChxyz,iCoset(jCoSet,isAtom))
                  If (iTest.eq.jTest) n = jCoset
 641           Continue
               If (n.lt.0 .or. n.gt.nCoset-1) Then
                  Call WarningMessage(2,' Error finding coset element')
                  Call Abend()
               End If
               iAdd(n) = iAdd(n) + iPrmt(iOper(iIrrep),iComp)
 640        Continue
            Do 645 jCoSet = 0, nCoSet-1
               If (iAdd(jCoSet).eq.0) Go To 611
 645        Continue
            nDimbc = nDimbc + 1
            Smmtrc(3*(isAtom-1)+i)=.True.
 611     Continue
 610  Continue
*                                                                      *
************************************************************************
*                                                                      *
*     Transform charges to masses (C=12)
*
      ii = ipCM
      Call GetMem('Mass','Allo','Real',ip_xMass,nsAtom)
      Call Get_Mass(Work(ip_xMass),nsAtom)
*     Call RecPrt(' Charges',' ',Work(ipCM),nsAtom,1)
      Call GetMem('ANr','Allo','Inte',ipANr,nsAtom)
      jj = ipANr
      Do 110 isAtom = 1, nsAtom
         ind = Int(Work(ii))
         If (ind.le.0) Then
*        If (ind.eq.0) Then
*           Work(ii) = Zero
*        Else If (ind.eq.-1) Then
*           Work(ii) = 1.0D99
            Work(ii) = 1.0D-10
         Else
            Work(ii) = Work(ip_xMass+isAtom-1)
         End If
         iWork(jj)=ind
         ii = ii + 1
         jj = jj + 1
 110  Continue
      Call GetMem('Mass','Free','Real',ip_xMass,nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the multiplicities of the cartesian coordinates.
*
      mTtAtm=0
      Do 4100 isAtom = 1, nsAtom
         iOff = 3*(isAtom-1) + ipCoor
         mTtAtm=mTtAtm+iDeg(Work(iOff),iOper,nSym)
         tmp = DBLE(iDeg(Work(iOff),iOper,nSym))
         i=(isAtom-1)*3+1
         Degen(i)=tmp
         i=(isAtom-1)*3+2
         Degen(i)=tmp
         i=(isAtom-1)*3+3
         Degen(i)=tmp
 4100 Continue
*     Call RecPrt('Degen',' ',Degen,1,3*nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*     Call qpg_dArray('Transverse',lRP,nRP)
*     If (lRP) Then
*        Call Allocate_Work(ipR12,nRP)
*        Call Get_dArray('Transverse',Work(ipR12),nRP)
*     End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute center of mass and molecular mass. The molecule is
*     translated so origin and center of mass is identical.
*
      If (jPrint.ge.99) Call
     &     Prlist('Symmetry Distinct Nuclear Coordinates / Bohr',
     &                   AtomLbl,nsAtom,Work(ipCoor),3,nsAtom)
      LWrite = .False.
      If (jPrint.ge.99) lWrite=.True.
      Call CofMss(Work(ipCoor),Work(ipCM),iOper,nSym,
     &            nsAtom,LWrite,cMass,iSym)
      LWrite = .False.
      If (jPrint.ge.99) Call
     &     PrList('Symmetry Distinct Nuclear Forces / au',
     &                   AtomLbl,nsAtom,Work(ipGrd),3,nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
      ip_B=ip_Dummy
      ip_dB=ip_Dummy
      ip_iB=ip_iDummy
      ip_idB=ip_iDummy
      ip_nqB=ip_iDummy
      mB_Tot=0
      mdB_Tot=0
      mq=0
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Init')
      Return
      End
