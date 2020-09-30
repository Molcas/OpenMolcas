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
* Copyright (C) 1991,1992, Roland Lindh                                *
************************************************************************
      SubRoutine Inputg(LuSpool)
************************************************************************
*                                                                      *
* Object: input module for the gradient code                           *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
*                                                                      *
*             Modified to complement GetInf, January '92.              *
************************************************************************
      use Alaska_Info, only: Am
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iChTbl, iOper, lIrrep, lBsFnc
      use Temporary_Parameters
      use Real_Info, only: CutInt
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "disp.fh"
#include "iavec.fh"
#include "stdalloc.fh"
#include "columbus_gamma.fh"
#include "exterm.fh"
#include "nac.fh"
#include "alaska_root.fh"
      Logical TstFnc, Type, Slct, T_Only, No_Input_OK
      Real*8, Allocatable:: Tmp(:), C(:,:), Scr(:,:), Temp(:,:)
      Integer, Allocatable:: IndCar(:)

      Character(LEN=1) :: xyz(0:2)=['x','y','z']
      Character(LEN=80) KWord, Key
      Integer iSym(3), iTemp(3*MxAtom)
      Logical timings,Reduce_Prt
      External Reduce_Prt
      COMMON /CHOTIME / timings

      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first

      COMMON  / OFembed_R1/ Xsigma
      COMMON  / OFembed_R2/ dFMD
      Character*16  OFE_KSDFT
      COMMON  / OFembed_C / OFE_KSDFT
*
      iRout = 99
      iPrint = nPrint(iRout)
      Do i = 1, nRout
         nPrint(i) = 5
      End Do
      If (ForceNAC) isNAC = .True.
      DoCSF  = .True.
      isCSF  = .False.
      Auto   = .False.
      Test   = .False.
      T_Only = .False.
      TRSymm = .False.
      lEq    = .False.
      Slct   = .False.
      l2DI   = .True.
      HF_Force=.False.
      NO_NUC = .False.
      State  =  0d0
      Timings_default = Timings
      Xsigma=1.0d4
      dFMD=0.0d0
      Do_OFemb=.false.
      KEonly=.false.
      OFE_first=.true.
      Show=.True.
      LuWr=6
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=iPL-1
      If (iPL.eq.0) Then
         jPrint=0
      Else If (iPL.eq.1) Then
         jPrint=0
      Else If (iPL.eq.2) Then
         jPrint=6
      Else If (iPL.eq.3) Then
         jPrint=6
      Else If (iPL.eq.4) Then
         jPrint=49
      Else
*     Else If (iPL.eq.5) Then
         jPrint=98
      End If
*
      Do i = 1, nRout
         nPrint(i) = jPrint
      End Do
*
*     First CutGrd can not be more accurate than CutInt!
      CutGrd = Max(1.0D-07,CutInt)
*     Second CutInt should now locally for Alaska be reset to the value
*     of CutInt/100!
      CutInt=CutGrd*1.0D-2
      Do 109 i = 1, 3*MxAtom
         IndxEq(i) = i
 109  Continue
      Do 1500 ldsp = 1, 3*MxAtom
         Direct(ldsp) = .True.
 1500 Continue
*                                                                      *
************************************************************************
*                                                                      *
*     KeyWord directed input
*
      Rewind(LuSpool)
      No_Input_OK=.True.
      Call RdNLst_(LuSpool,'ALASKA',No_Input_OK)
      KWord=' &ALASKA'
 998  Read (LuSpool,'(A72)',END=997,ERR=988) Key
      KWord = Key
      Call UpCase(KWord)
      If (KWord(1:1).eq.'*')    Go To 998
      If (KWord.eq.'')       Go To 998
      If (KWord(1:4).eq.'VERB') Go To 912
      If (KWord(1:4).eq.'PRIN') Go To 930
      If (KWord(1:4).eq.'EQUI') Go To 935
      If (KWord(1:4).eq.'CUTO') Go To 942
      If (KWord(1:4).eq.'HF-F') Go To 993
      If (KWord(1:4).eq.'NOIN') Go To 953
      If (KWord(1:4).eq.'SELE') Go To 960
      If (KWord(1:4).eq.'2DOP') Go To 965
      If (KWord(1:4).eq.'2DIP') Go To 966
      If (KWord(1:4).eq.'ONEO') Go To 990
      If (KWord(1:4).eq.'TEST') Go To 991
      If (KWord(1:4).eq.'SHOW') Go To 992
      If (KWord(1:4).eq.'PNEW') Go To 994
      If (KWord(1:4).eq.'POLD') Go To 995
      If (KWord(1:4).eq.'NONU') Go To 996
      If (KWord(1:4).eq.'EXTR') Go To 971
      If (KWord(1:4).eq.'CHOI') Go To 972
      If (KWord(1:4).eq.'OFEM') Go To 973
      If (KWord(1:4).eq.'KEON') Go To 974
      If (KWord(1:4).eq.'DFMD') Go To 975
* Keyword 'NUMErical' checked earlier - forces numerical gradients
* Keyword 'DELTa' selects the scaling factor for the displacements
*                 in the numerical_gradient module
* Keyword 'KEEP' does not remove the old gradient
* Keyword 'INVErt' inverts the treatment of constraints
* Here it's only included for consistency
      If (KWord(1:4).eq.'NUME') Go To 998
      If (KWord(1:4).eq.'DELT') Go To 998
      If (KWord(1:4).eq.'KEEP') Go To 998
      If (KWord(1:4).eq.'INVE') Go To 998
      If (KWord(1:4).eq.'ROOT') Go To 976
      If (KWord(1:4).eq.'NAC ') Go To 977
      If (KWord(1:4).eq.'NOCS') Go To 978
      If (KWord(1:4).eq.'AUTO') Go To 979
      If (KWord(1:4).eq.'END ') Go To 997
      Call WarningMessage(2,'Error in InputG')
      Write (LuWr,*) 'Inputg: Illegal keyword'
      Write (LuWr,'(A,A)') 'KWord=',KWord
      Call Quit_OnUserError()
*
 988  Call WarningMessage(2,'Error in InputG')
      Write (LuWr,*) 'Inputg: Error reading the input'
      Write (LuWr,'(A,A)') 'Last read line=',KWord
      Call Quit_OnUserError()
*                                                                      *
************************************************************************
*                                                                      *
*     Print level
*
 930  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 930
      If (KWord.eq.'')    Go To 930
      Read(KWord,*,Err=988) n
      Do 931 i = 1, n
 9301    Read(LuSpool,'(A)',Err=988) KWord
         If (KWord(1:1).eq.'*') Go To 9301
         If (KWord.eq.'')    Go To 9301
         Read(KWord,*,Err=988) jRout, iPrint
         nPrint(jRout)=iPrint
 931  Continue
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     Equivalence option
*
 935  Continue
      If (T_Only) Then
         Call WarningMessage(2,'Error in InputG')
         Write (LuWr,*)'EQUI option does not ork with RF calculations!'
         Call Quit_OnUserError()
      End If
      lEq=.True.
 936  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 936
      If (KWord.eq.'')    Go To 936
      Read(KWord,*) nGroup
      Do 937 iGroup = 1, nGroup
 938     Read(LuSpool,'(A)',Err=988) KWord
         If (KWord(1:1).eq.'*') Go To 938
         If (KWord.eq.'')    Go To 938
         Read(KWord,*) nElem,(iTemp(iElem),iElem=1,nElem)
         Do 939 iElem=2,nElem
            IndxEq(iTemp(iElem)) = iTemp(1)
            Direct(iTemp(iElem)) = .False.
939      Continue
937   Continue
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     Cutoff for computing primitive gradients
*
 942  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 942
      If (KWord.eq.'')    Go To 942
      Read(KWord,*,Err=988) CutGrd
      CutGrd = Abs(CutGrd)
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     Disable the utilization of translational and
*     rotational invariance of the energy in the
*     computation of the molecular gradient.
*
 953  TRSymm=.False.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     selection option
*
 960  Continue
      If (T_Only) Then
         Call WarningMessage(2,'Error in InputG')
         Write (LuWr,*)'SELE option does not work with RF calculations!'
         Call Quit_OnUserError()
      End If
      Slct = .True.
      If (lEq) Then
         Call WarningMessage(2,'Error in InputG')
         Write (LuWr,*) ' The Selection option must preceed the',
     &                  ' Equivalence option to work together.'
         Call Quit_OnUserError()
      End If
      Do 961 i = 1, 3*MxAtom
         Direct(i) = .False.
 961  Continue
 962  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 962
      If (KWord.eq.'')    Go To 962
      Read(KWord,*) nSlct
*
 963  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 963
      If (KWord.eq.'')    Go To 963
      Read(KWord,*) (iTemp(iElem),iElem=1,nSlct)
      Do 964 iElem=1,nSlct
         Direct(iTemp(iElem)) = .True.
964   Continue
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     Change default for the prescreening.
*
 965  l2DI  = .False.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     Change default for the prescreening.
*
 966  l2DI  = .True.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     Do not compute two electron integrals.
*
 990  Onenly = .TRUE.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     Process only the input.
*
 991  Test = .TRUE.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*-----Raise the printlevel to show gradient contributions
*
 992  Continue
      If (iPL.ge.2) Then
         nPrint(112) = 15
         nPrint(1)   = 15
         nPrint(33)  = 15
      End If
      Go To 998
*                                                                      *
****** PNEW ************************************************************
*                                                                      *
*
*-----Print gradient in NEW human-readable format
*
 994  Continue
      nPrint(1)   =  4
      Go To 998
*                                                                      *
****** POLD ************************************************************
*                                                                      *
*
*-----Print gradient in OLD format
*
 995  Continue
      nPrint(1)   =  5
      Go To 998
*                                                                      *
****** VERB ************************************************************
*                                                                      *
*
*----- Verbose mode.
*
 912  Continue
      nPrint(80)  =  6
      nPrint( 1)  =  6
      nPrint( 9)  =  6
      nPrint(99)  =  6
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*
*     Compute Hellmann-Feynman forces
*
 993  HF_Force = .TRUE.
      Go To 998
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*
*     Do not compute the nuclear charge contribution
*
 996  NO_NUC = .TRUE.
      Go To 998

************************************************************************
*                                                                      *
*
*     Put the program name and the time stamp onto the extract file
*
971   Write (LuWr,*)'InputG: EXTRACT option is redundant and is',
     &              ' ignored!'
      Go To 998
*                                                                      *
************************************************************************
*
*     Cholesky input section
*
 972  Continue
      Call Cho_alaska_rdInp(LuSpool)
      Go To 998
*                                                                      *
************************************************************************
*
*     Orbital-Free Embedding (OFE) input section
*
 973  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 973
      If (KWord.eq.'')    Go To 973
      Call UpCase(KWord)
      Call LeftAd(KWord)
      Read(KWord,'(A)') OFE_KSDFT
      Do_OFemb=.true.
      Go To 998
*                                                                      *
************************************************************************
*
*     Mode "Kinetic Energy Only" for OFE input section
*
 974  Continue
      KEonly=.true.
      Go To 998
*                                                                      *
************************************************************************
*
*     Mode "Kinetic Energy Only" for OFE input section
*
 975  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 975
      If (KWord.eq.'')    Go To 975
      Read(KWord,*) dFMD, Xsigma
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     Root keyword, now also for analytical gradient
*     This is a dummy, the keyword is already read in chk_numerical
*
 976  Read(LuSpool,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 976
      If (KWord.eq.'')    Go To 976
      Read(KWord,*) iRoot
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     NAC keyword: compute non-adiabatic couplings between 2 states
*     The keyword is already read in chk_numerical
*
 977  Read(LuSpool,'(A)',Err=988) KWord
      isNAC=.True.
      If (KWord(1:1).eq.'*') Go To 977
      If (KWord.eq.'')    Go To 977
      Read(KWord,*) NACstates(1),NACstates(2)
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     NOCSF keyword, to neglect the CSF contribution to the NAC,
*     which is the cause for translational variance
*
 978  DoCSF=.False.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     AUTO keyword, used by SLAPAF, to signal this is an automated
*     call to ALASKA
*
 979  Auto=.True.
      Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*                          End of input section.                       *
*                                                                      *
************************************************************************
 997  Continue
*
* NAC could have been activated through explicit input or through
* a previous MCLR
*
      If (isNAC) Then
         No_Nuc=.True.
*Get the state energies
         Call Get_iScalar('Number of roots',nRoots)
         Call mma_Allocate(Tmp,nRoots,Label='Tmp')
         Call Get_dArray('Last energies',Tmp,nRoots)
         Ediff=Tmp(NACstates(1))-Tmp(NACstates(2))
         Call mma_deallocate(Tmp)
      End If
*
      nCnttp_Valence=0
      Do iCnttp = 1, nCnttp
         If (dbsc(iCnttp)%Aux) Exit
         nCnttp_Valence = nCnttp_Valence+1
      End Do
*
      If (lEq)  TRSymm=.False.
      If (Slct) TRSymm=.False.
      iPrint=nPrint(iRout)
*
      TRsymm=(TRsymm.or.T_Only) .and. .Not.Test
*
*---- Compute number of centers and displacements. Ignore pseudo centers.
*     If any pseudo centers disable use of translational and rotational
*     invariance.
*
      mDisp = 0
      mdc = 0
      Do 10 iCnttp = 1, nCnttp_Valence
         If (dbsc(iCnttp)%pChrg) Then
             TRSymm=.False.
             mdc = mdc + dbsc(iCnttp)%nCntr
             Go To 10
         Else If(dbsc(iCnttp)%nFragType.gt.0.or.dbsc(iCnttp)%Frag) Then
           TRSymm = .false.
         End If
         Do 20 iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            mDisp = mDisp + 3*(nIrrep/dc(mdc)%nStab)
 20      Continue
 10   Continue
*
      If (HF_Force.and.Show.and.iPrint.ge.6) Then
         Write (LuWr,*)
         Write (LuWr,'(A)') '            O B S E R V E ! '
         Write (LuWr,'(A)') '            Option for computation of '//
     &                   'interstate couling vector or'
         Write (LuWr,'(A)') '            Hellmann-Feynman gradient '//
     &                   'is active.'

         Write (LuWr,*)
      End If
      If (Show.and.iPrint.ge.6) Then
         Write (LuWr,*)
         Write (LuWr,'(20X,A,E10.3)')
     &     ' Threshold for contributions to the gradient:',CutGrd
         Write (LuWr,*)
      End If
*
*     Generate symmetry adapted cartesian displacements
*
      If (Show.and.iPrint.ge.6) Then
         Write (LuWr,*)
         Write (LuWr,'(20X,A)')
     &              '********************************************'
         Write (LuWr,'(20X,A)')
     &              '* Symmetry Adapted Cartesian Displacements *'
         Write (LuWr,'(20X,A)')
     &           '********************************************'
         Write (LuWr,*)
      End If
*
      Call ICopy(MxAtom*8,[0],0,IndDsp,1)
      Call ICopy(MxAtom*3,[0],0,InxDsp,1)
      call dcopy_(3*MxSym*MxAtom,[One],0,Disp_Fac,1)
      Call ICopy(3*MxAtom,[1],0,mult_Disp,1)
      nDisp = 0
      Do iIrrep = 0, nIrrep-1
         lDisp(iIrrep) = 0
         Type = .True.
*        Loop over basis function definitions
         mdc = 0
         mc = 1
         Do iCnttp = 1, nCnttp_Valence
*           Loop over unique centers associated with this basis set.
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               IndDsp(mdc,iIrrep) = nDisp
*              Loop over the cartesian components
               Do iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(dc(mdc)%iCoSet,
     &                iIrrep,iComp,dc(mdc)%nStab) .and.
     &                .Not.dbsc(iCnttp)%pChrg ) Then
                      nDisp = nDisp + 1
                      If (iIrrep.eq.0) InxDsp(mdc,iCar+1) = nDisp
                      lDisp(iIrrep) = lDisp(iIrrep) + 1
                      mult_Disp(nDisp)=nIrrep/dc(mdc)%nStab
                      If (Type) Then
      If (Show.and.iPrint.ge.6) then
                         Write (LuWr,*)
                         Write (LuWr,'(10X,A,A)')
     &                    ' Irreducible representation : ',
     &                      lIrrep(iIrrep)
                         Write (LuWr,'(10X,2A)')
     &                      ' Basis function(s) of irrep: ',
     &                       lBsFnc(iIrrep)
                         Write (LuWr,*)
                         Write (LuWr,'(A)')
     &                   ' Basis Label        Type   Center Phase'
      End If
                         Type = .False.
                      End If
                      If (iIrrep.eq.0) Then
                         Do jOper = 0, nIrrep-1
                            Disp_Fac(iCar+1,jOper,mdc)=
     &                        DBLE(iPrmt( jOper ,iComp) *
     &                             iChTbl(iIrrep,jOper))
                         End Do
                      End If
      If (Show.and.iPrint.ge.6) Then
                      Write (LuWr,'(I4,3X,A8,5X,A1,7X,8(I3,4X,I2,4X))')
     &                      nDisp,dc(mdc)%LblCnt,xyz(iCar),
     &                      (mc+iCo,iPrmt(
     &                      NrOpr(dc(mdc)%iCoSet(iCo,0)),iComp)*
     &                      iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0))),
     &                      iCo=0,nIrrep/dc(mdc)%nStab-1 )
      End If
      Write (ChDisp(nDisp),'(A,1X,A1)')
     &      dc(mdc)%LblCnt,xyz(iCar)
                  End If
*
               End Do
               mc = mc + nIrrep/dc(mdc)%nStab
            End Do
         End Do
*
      End Do
*
      If (nDisp.ne.mDisp) Then
         Call WarningMessage(2,'Error in InputG')
         Write (LuWr,*)
     &      ' Wrong number of symmetry adapted displacements',
     &       nDisp,'=/=',mDisp
         Call Abend()
      End If
*
*     Set up data for the utilization of the translational
*     and rotational invariance of the energy.
*
      If (TRSymm) Then
         iSym(1) = 0
         iSym(2) = 0
         iSym(3) = 0
         Do i = 1, Min(nIrrep-1,5)
            j = i
            If (i.eq.3) j = 4
            Do k = 1, 3
               If (iAnd(iOper(j),2**(k-1)).ne.0) iSym(k) = 2**(k-1)
            End Do
         End Do
         nTR = 0
*--------Translational equations
         Do i = 1, 3
            If (iSym(i).eq.0) nTR = nTR + 1
         End Do
         If (iPrint.ge.99) Write (LuWr,*) ' nTR=',nTR
*--------Rotational equations
         If (.Not.T_Only) Then
            Do i = 1,3
               j = i+1
               If (j.gt.3) j = j-3
               k = i+2
               If (k.gt.3) k = k-3
               ijSym = iEor(iSym(j),iSym(k))
               If (ijSym.eq.0) nTR = nTR + 1
            End Do
         End If
         If (nTR.eq.0) Then
            TRSymm = .False.
            Go To 9876
         End If
         If (iPrint.ge.99) Write (LuWr,*) ' nTR=',nTR
         Call mma_allocate(Am,nTR,lDisp(0),Label='Am')
         Call mma_allocate(Temp,nTR,nTR,Label='Temp')
         Call mma_allocate(C,4,lDisp(0),Label='C')
         Call mma_allocate(IndCar,lDisp(0),Label='IndCar')
*
         Am(:,:)=Zero
         C(:,:)=Zero
*
*        Generate temporary information of the symmetrical
*        displacements.
*
         ldsp = 0
         mdc = 0
         iIrrep = 0
         Do 2100 iCnttp = 1, nCnttp_Valence
            Do 2200 iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
*              Call RecPrt(' Coordinates',' ',
*    &                     dbsc(iCnttp)%Coor(1,iCnt),1,3)
               Fact = Zero
               iComp = 0
               If (dbsc(iCnttp)%Coor(1,iCnt).ne.Zero)
     &             iComp = iOr(iComp,1)
               If (dbsc(iCnttp)%Coor(2,iCnt).ne.Zero)
     &             iComp = iOr(iComp,2)
               If (dbsc(iCnttp)%Coor(3,iCnt).ne.Zero)
     &             iComp = iOr(iComp,4)
               Do jIrrep = 0, nIrrep-1
                  If ( TstFnc(dc(mdc)%iCoSet,
     &                  jIrrep,iComp,dc(mdc)%nStab) ) Then
                     Fact = Fact + One
                  End If
               End Do
               Do iCar = 1, 3
                  iComp = 2**(iCar-1)
                  If ( TstFnc(dc(mdc)%iCoSet,
     &                  iIrrep,iComp,dc(mdc)%nStab) ) Then
                     ldsp = ldsp + 1
                     Direct(lDsp)=.True.
*--------------------Transfer the coordinates
                     call dcopy_(3,dbsc(iCnttp)%Coor(1:3,iCnt),1,
     &                           C(1,ldsp),1)
*--------------------Transfer the multiplicity factor
                     C(4,ldsp)= Fact
                     IndCar(ldsp) = iCar
                  End If
                End Do
 2200       Continue
 2100    Continue
         If (iPrint.ge.99) Then
            Call RecPrt(' Information',' ',C,4,lDisp(0))
            Write (LuWr,*) (IndCar(i),i=1,lDisp(0))
         End If
*
*--------Set up coefficient for the translational equations
*
         iTR = 0
         Do i = 1,3
            If (iSym(i).eq.0) Then
               iTR = iTR + 1
               Do ldsp = 1, lDisp(0)
                  If (IndCar(ldsp).eq.i) Then
                     Am(iTR,ldsp) = C(4,ldsp)
                  End If
               End Do
            End If
         End Do
*
*--------Set up coefficient for the rotational invariance
*
         If (.Not.T_Only) Then
            Do i = 1, 3
               j = i + 1
               If (j.gt.3) j = j - 3
               k = i + 2
               If (k.gt.3) k = k - 3
               ijSym = iEor(iSym(j),iSym(k))
               If (ijSym.ne.0) Go To 1210
               iTR = iTR + 1
               Do ldsp = 1, lDisp(0)
                  If (IndCar(ldsp).eq.j) Then
                     Fact = C(4,ldsp) * C(k,ldsp)
                     Am(iTR,ldsp) = Fact
                  Else If (IndCar(ldsp).eq.k) Then
                     Fact =-C(4,ldsp) * C(j,ldsp)
                     Am(iTR,ldsp) = Fact
                  End If
               End Do
 1210          Continue
            End Do
         End If
         If (iPrint.ge.99)
     &      Call RecPrt(' The A matrix',' ',Am,nTR,lDisp(0))
*
*--------Now, transfer the coefficient of those gradients which will
*        not be computed directly.
*        The matrix to compute the inverse of is determined via
*        a Gram-Schmidt procedure.
*
*--------Pick up the other vectors
         Do iTR = 1, nTR
*           Write (LuWr,*) ' Looking for vector #',iTR
            ovlp = Zero
            kTR = 0
*-----------Check all the remaining vectors
            Do 1231 ldsp = 1, lDisp(0)
               Do jTR = 1, iTR-1
                  If (iTemp(jTR).eq.ldsp) Go To 1231
               End Do
*              Write (LuWr,*) ' Checking vector #', ldsp
               call dcopy_(nTR,Am(1,ldsp),1,Temp(1,iTR),1)
*              Call RecPrt(' Vector',' ',Temp(1,iTR),nTR,1)
*--------------Gram-Schmidt orthonormalize against accepted vectors
               Do lTR = 1, iTR-1
                  alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,lTR),1)
*                 Write (LuWr,*) ' <x|y> =', alpha
                  Call DaXpY_(nTR,-alpha,Temp(1,lTR),1,Temp(1,iTR),1)
               End Do
*              Call RecPrt(' Remainings',' ',Temp(1,iTR),nTR,1)
               alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,iTR),1)
*              Write (LuWr,*) ' Remaining overlap =', alpha
*--------------Check the remaining magnitude of vector after Gram-Schmidt
               If (alpha.gt.ovlp) Then
                  kTR = ldsp
                  ovlp = alpha
               End If
               If (.Not.Direct(ldsp).and.alpha.gt.1.0D-2) Then
                  kTR = ldsp
                  ovlp = 1.0D99
               End If
 1231       Continue
            If (kTR.eq.0) Then
               Call WarningMessage(2,'Error in InputG')
               Write (LuWr,*) ' No Vector found!'
               Call Abend()
            End If
*           Write (LuWr,*) ' Selecting vector #', kTR
*-----------Pick up the "best" vector
            call dcopy_(nTR,Am(1,kTR),1,Temp(1,iTR),1)
            Do lTR = 1, iTR-1
               alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,lTR),1)
               Call DaXpY_(nTR,-alpha,Temp(1,lTR),1,Temp(1,iTR),1)
            End Do
            alpha = DDot_(nTR,Temp(1,iTR),1,Temp(1,iTR),1)
            Call DScal_(nTR,One/Sqrt(alpha),Temp(1,iTR),1)
            iTemp(iTR) = kTR
         End Do
         Do iTR = 1, nTR
            call dcopy_(nTR,Am(1,iTemp(iTR)),1,Temp(1,iTR),1)
            Am(:,iTemp(iTR))=Zero
         End Do
         If (iPrint.ge.99) Then
            Call RecPrt(' The A matrix',' ',Am,nTR,lDisp(0))
            Call RecPrt(' The T matrix',' ',Temp,nTR,nTR)
            Write (LuWr,*) (iTemp(iTR),iTR=1,nTR)
         End If
*
*        Compute the inverse of the T matrix
*
         Call MatInvert(Temp,nTR)
         If (IPrint.ge.99)
     &      Call RecPrt(' The T-1 matrix',' ',Temp,nTR,nTR)
         Call DScal_(nTR**2,-One,Temp,1)
*
*        Generate the complete matrix
*
         Call mma_allocate(Scr,nTR,lDisp(0),Label='Scr')
         Call DGEMM_('N','N',
     &               nTR,lDisp(0),nTR,
     &               1.0d0,Temp,nTR,
     &               Am,nTR,
     &               0.0d0,Scr,nTR)
         If (IPrint.ge.99) Call RecPrt(' A-1*A',' ',Scr,nTR,lDisp(0))
         Call mma_deallocate(Am)
         Call mma_allocate(Am,lDisp(0),lDisp(0),Label='Am')
         Am(:,:)=Zero
         Do i = 1, lDisp(0)
            Am(i,i)=One
         End Do
         Do 1250 iTR = 1, nTR
            ldsp = iTemp(iTR)
            call dcopy_(lDisp(0),Scr(1,iTR),nTR,Am(1,lDisp),lDisp(0))
 1250    Continue
         If (iPrint.ge.99)
     &      Call RecPrt('Final A matrix',' ',Am,lDisp(0),lDisp(0))
*
*
         Call mma_deallocate(Scr)
         Call mma_deallocate(IndCar)
         Call mma_deallocate(C)
         Call mma_deallocate(Temp)
         Do 1501 iTR = 1, nTR
            ldsp = iTemp(iTR)
            Direct(ldsp)=.False.
 1501    Continue
*
         Write (LuWr,*)
         Write (LuWr,'(20X,A,A)')
     &      ' Automatic utilization of translational and',
     &      ' rotational invariance of the energy is employed.'
         Write (LuWr,*)
         Do 7000 i = 1, lDisp(0)
            If (Direct(i)) Then
               Write (LuWr,'(25X,A,A)') Chdisp(i), ' is independent'
            Else
               Write (LuWr,'(25X,A,A)') Chdisp(i), ' is dependent'
            End If
 7000    Continue
         Write (LuWr,*)
*
      Else
         nTR = 0
         If (Show.and.iPrint.ge.6) Then
            Write (LuWr,*)
            Write (LuWr,'(20X,A,A)')
     &         ' No automatic utilization of translational and',
     &         ' rotational invariance of the energy is employed.'
            Write (LuWr,*)
         End If
      End If
*
      If (Slct) Then
         Write (LuWr,*)
         Write (LuWr,'(20X,A)') ' The Selection option is used'
         Write (LuWr,*)
         Do 7100 i = 1, lDisp(0)
            If (Direct(i)) Then
               Write (LuWr,'(25X,A,A)') Chdisp(i), ' is computed'
            Else
               Write (LuWr,'(25X,A,A)') Chdisp(i), ' is set to zero'
            End If
 7100    Continue
         Write (LuWr,*)
      End If
*
 9876 Continue
*
*     Set up the angular index vector
*
      i = 0
      Do 1000 iR = 0, iTabMx
         Do 2000 ix = iR, 0, -1
            Do 3000 iy = iR-ix, 0, -1
               iz = iR-ix-iy
               i = i + 1
               ixyz(1,i) = ix
               ixyz(2,i) = iy
               ixyz(3,i) = iz
 3000       Continue
 2000    Continue
 1000 Continue
*
      Onenly = HF_Force
*
      Return
      End
