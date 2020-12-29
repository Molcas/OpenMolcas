***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Init_SlapAf()
      use Symmetry_Info, only: nIrrep, iOper
      use Slapaf_Info, only: q_nuclear, dMass, Coor, Grd, ANr, Degen,
     &                       jStab, nStab, iCoSet, AtomLbl, Smmtrc,
     &                       RootMap
*     use Slapaf_Info, only: R12
      use Slapaf_Parameters, only: nDimBC, Analytic_Hessian, MaxItr,
     &                             Line_Search, ThrEne, ThrGrd, ThrCons,
     &                             ThrMEP, Header, MxItr, mTtAtm
*     use Slapaf_Parameters, only: lRP
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "sbs.fh"
#include "nadc.fh"
#include "db.fh"
#include "print.fh"
#include "stdalloc.fh"
      Integer   iAdd(0:7)
      Integer :: jPrmt(0:7)=[1,-1,-1,1,-1,1,1,-1]
      Logical Same, Do_ESPF, Exist_2, Found, Reduce_Prt
      External Reduce_Prt
      Character(LEN=8) CMAX
      Integer Columbus
#include "SysDef.fh"
      Character(LEN=100) SuperName
      Character(LEN=100), External:: Get_SuperName
      Real*8, Allocatable:: xMass(:)
*
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
      jPrint=10
*                                                                      *
************************************************************************
*                                                                      *
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
      ThrMEP = ThrGrd
      ThrCons = 1.0D10
      nWndw = 5
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
*     Get Molecular data
*
*...  Read the title
*
      Call Get_cArray('Seward Title',Header,144)
*
*...  Read number of atoms, charges, coordinates, gradients and
*     atom labels
*
      Call Get_Molecule()
*                                                                      *
************************************************************************
*                                                                      *
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
         Call mma_allocate(RootMap,nRM,Label='RootMap')
         Call Get_iArray('Root Mapping',RootMap,nRM)
      Else
         Call Qpg_iScalar('Number of roots',Found)
         nRoots = 1
         If (Found) Call Get_iScalar('Number of roots',nRoots)
         Call mma_allocate(RootMap,nRoots,Label='RootMap')
         RootMap(:)=Zero
         Do i=1,nRoots
            RootMap(i)=i
         End Do
      End If
*
*...  Check if there is an analytic Hessian
      Call qpg_dArray('Analytic Hessian',Analytic_Hessian,nHess)

      If (.Not.Analytic_Hessian) Then
         Call NameRun('RUNOLD')
         Call qpg_dArray('Analytic Hessian',Analytic_Hessian,nHess)
         Call NameRun('#Pop')
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute the number of total symmetric displacements
*
      Call mma_allocate(jStab ,[0,7],[1,SIZE(Coor,2)],Label='jStab ')
      Call mma_allocate(nStab ,      [1,SIZE(Coor,2)],Label='nStab ')
      Call mma_allocate(iCoSet,[0,7],[1,SIZE(Coor,2)],Label='iCoSet')
      Call mma_allocate(Smmtrc,3,SIZE(Coor,2),Label='Smmtrc')
      jStab(:,:)=0
      nStab(:)=0
      iCoSet(:,:)=0
      Smmtrc(:,:)=.False.

      nDimbc = 0
*...  Loop over the unique atoms
      Do 610 isAtom = 1, SIZE(Coor,2)
*...     Find character of center
         iChxyz=0
         Do i = 1, 3
            If (Coor(i,isAtom).ne.Zero) Then
               Do iIrrep= 0, nIrrep-1
                  If (iAnd(2**(i-1),iOper(iIrrep)).ne.0)
     &               iChxyz=iOr(iChxyz,2**(i-1))
               End Do
            End If
         End Do
         nStb = 0
         Do iIrrep = 0, nIrrep-1
            If (iAnd(iChxyz,iOper(iIrrep)).eq.0) Then
               jStab(nStb,isAtom)=iOper(iIrrep)
               nStb = nStb + 1
            End If
         End Do
         nStab(isAtom)=nStb
*...     Find the coset representatives
         iCoSet(0,SIZE(Coor,2)) = 0      ! Put in the unit operator
         nCoSet = 1
         Do iIrrep = 1, nIrrep-1
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
         If (nIrrep/nStb.ne.nCoSet) Then
            Call WarningMessage(2,' Error while doing cosets.')
            Call Abend()
         End If
         Do 611 i = 1, 3
            iComp = 2**(i-1)
            Call ICopy(nCoSet,[0],0,iAdd,1)
            Do 640 iIrrep = 0, nIrrep-1
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
               iAdd(n) = iAdd(n) + jPrmt(iAnd(iOper(iIrrep),iComp))
 640        Continue
            Do 645 jCoSet = 0, nCoSet-1
               If (iAdd(jCoSet).eq.0) Go To 611
 645        Continue
            nDimbc = nDimbc + 1
            Smmtrc(i,isAtom)=.True.
 611     Continue
 610  Continue
*                                                                      *
************************************************************************
*                                                                      *
*     Transform charges to masses (C=12)
*
      Call mma_allocate(dMass,SIZE(Coor,2),Label='dMass')
      Call mma_allocate(xMass,SIZE(Coor,2),Label='xMass')
      Call Get_Mass(xMass,SIZE(Coor,2))
*     Call RecPrt(' Charges',' ',Q_nuclear,SIZE(Coor,2),1)
      Call mma_allocate(ANr,SIZE(Coor,2),Label='ANr')
      Do isAtom = 1, SIZE(Coor,2)
         ind = Int(Q_nuclear(isAtom))
         If (ind.le.0) Then
            dMass(isAtom) = 1.0D-10
         Else
            dMass(isAtom) = xMass(isAtom)
         End If
         ANr(isAtom)=ind
      End Do
      Call mma_deallocate(xMass)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the multiplicities of the cartesian coordinates and the
*     total number of atoms.
*
      mTtAtm=0
      Call mma_Allocate(Degen,3,SIZE(Coor,2),Label='Degen')
      Do isAtom = 1, SIZE(Coor,2)
         mTtAtm=mTtAtm+iDeg(Coor(:,isAtom))
         tmp = DBLE(iDeg(Coor(:,isAtom)))
         Do i = 1, 3
            Degen(i,isAtom)=tmp
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Degen',' ',Degen,3,SIZE(Coor,2))
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Call qpg_dArray('Transverse',lRP,nRP)
*     If (lRP) Then
*        Call mma_allocate(R12,3,nRP/3,Label='R12')
*        Call Get_dArray('Transverse',R12,nRP)
*     End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute center of mass and molecular mass. The molecule is
*     translated so origin and center of mass is identical.
*
      If (jPrint.ge.99) Call
     &     Prlist('Symmetry Distinct Nuclear Coordinates / Bohr',
     &                   AtomLbl,SIZE(Coor,2),Coor,3,SIZE(Coor,2))
      If (jPrint.ge.99) Call
     &     PrList('Symmetry Distinct Nuclear Forces / au',
     &                   AtomLbl,SIZE(Coor,2),Grd,3,SIZE(Coor,2))
*                                                                      *
************************************************************************
*                                                                      *
      mB_Tot=0
      mdB_Tot=0
      mq=0
      Force_dB=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
