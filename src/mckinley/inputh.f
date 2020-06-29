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
*               1996, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Inputh(Run_MCLR)
************************************************************************
*                                                                      *
* Object: input module for the gradient code                           *
*                                                                      *
* Called from: McKinley                                                *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              RecPrt                                                  *
*              DaXpY   (ESSL)                                          *
*              DDot_   (ESSL)                                          *
*              DScal   (ESSL)                                          *
*              DGEMM_  (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September 1991                                           *
*                                                                      *
*             Modified to complement GetInf, January 1992              *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "disp.fh"
#include "disp2.fh"
#include "iavec.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "SysDef.fh"
      Logical TstFnc, Type, Slct
c      Logical DoCholesky
      Character*1 xyz(0:2)
      Character*8 Label,labelop
      Character*32 Label2
      character*4 lab
      Logical Run_MCLR
      Character*80  KWord, Key
      Integer iSym(3), iTemp(3*mxdc)
      Data xyz/'x','y','z'/
*
      Call QEnter('InputH')

c      Call DecideOnCholesky(DoCholesky)
c      If (DoCholesky) Then
c       write(6,*)'** Cholesky or RI/DF not yet implemented in McKinley '
c       call abend()
c      EndIf

      iRout=99
      Do i = 1, nRout
         nPrint(i) = 5
      End Do
      show=.false.
      Onenly = .False.
      nmem=0
      Test   = .False.
      TRSymm = .false.
      lEq    = .False.
      Slct   = .False.
      PreScr   = .true.
      lGrd=.True.
      lHss=.True.
      Nona=.false.
      Run_MCLR=.True.
      CutInt = 1.0D-07
      ipert=2
      iCntrl=1
      Call lCopy(mxpert,[.true.],0,lPert,1)
      sIrrep=.false.
      iprint=0
      Do 109 i = 1, 3*mxdc
         IndxEq(i) = i
 109  Continue
*
*     KeyWord directed input
*
      LuRd=5
      Call RdNLst(LuRd,'MCKINLEY')
 998  Read(5,'(A72)',END=977,ERR=988) Key
      KWord = Key
      Call UpCase(KWord)
      If (KWord(1:1).eq.'*')    Go To 998
      If (KWord.eq.BLine)       Go To 998
*     If (KWord(1:4).eq.'EQUI') Go To 935
*     If (KWord(1:4).eq.'MEMO') Go To 951
*     If (KWord(1:4).eq.'NOTR') Go To 952
*     If (KWord(1:4).eq.'NOIN') Go To 953
      If (KWord(1:4).eq.'SHOW') Go To 992
      If (KWord(1:4).eq.'MEM ') Go To 697
*
      If (KWord(1:4).eq.'CUTO') Go To 942
      If (KWord(1:4).eq.'VERB') Go To 912
      If (KWord(1:4).eq.'NOSC') Go To 965
      If (KWord(1:4).eq.'ONEO') Go To 990
      If (KWord(1:4).eq.'SELE') Go To 960
      If (KWord(1:4).eq.'REMO') Go To 260
      If (KWord(1:4).eq.'PERT') Go To 975
      If (KWord(1:4).eq.'TEST') Go To 991
      If (KWord(1:4).eq.'EXTR') Go To 971
      If (KWord(1:4).eq.'NONA') Go To 972
      If (KWord(1:4).eq.'NOMC') Go To 973
      If (KWord(1:4).eq.'END ') Go To 997
      Write (6,*) 'InputH: Illegal keyword'
      Write (6,'(A,A)') 'KWord=',KWord
      Call QTrace
      Call Abend()
 977  Write (6,*) 'InputH: end of input file.'
      Write (6,'(A,A)') 'Last command=',KWord
      Call QTrace
      Call Abend()
 988  Write (6,*) 'InputH: error reading input file.'
      Write (6,'(A,A)') 'Last command=',KWord
      Call QTrace
      Call Abend()
*                                                                      *
****** MEM  ************************************************************
*                                                                      *
 697  Read(5,*) nmem
      nmem=nmem*1048576/rtob
      goto 998
*                                                                      *
****** PERT ************************************************************
*                                                                      *
*     Select which part of the Hessian will be compiuted.
*
975   Read(5,'(A)',Err=988) KWord
      If (KWord(1:1).eq.'*') Go To 975
      If (KWord.eq.BLine)    Go To 975
      Call UpCase(KWord)
      Lab=KWORD(1:4)
      If (KWORD(1:4).eq.'HESS') Then
         ipert=2
      Else If (KWORD(1:4).eq.'GEOM') Then
         ipert=1
      Else
         Write (6,*) 'InputH: Illegal perturbation keyword'
         Write (6,'(A,A)') 'KWord=',KWord
         Call QTrace
         Call Abend()
      End If

      Goto 998
*                                                                      *
****** EQUI ************************************************************
*                                                                      *
*     Equivalence option
*
*935  Continue
*     lEq=.True.
*936  Read(5,'(A)',Err=988) KWord
*     If (KWord(1:1).eq.'*') Go To 936
*     If (KWord.eq.BLine)    Go To 936
*     Read(KWord,*) nGroup
*     Do 937 iGroup = 1, nGroup
*938     Read(5,'(A)',Err=988) KWord
*        If (KWord(1:1).eq.'*') Go To 938
*        If (KWord.eq.BLine)    Go To 938
*        Read(KWord,*) nElem,(iTemp(iElem),iElem=1,nElem)
*        Do 939 iElem=2,nElem
*           IndxEq(iTemp(iElem)) = iTemp(1)
*           Direct(iTemp(iElem)) = .False.
*939      Continue
*937   Continue
*     Go To 998
*                                                                      *
****** CUTO ************************************************************
*                                                                      *
*     Cuttoff for computing primitive gradients
*
 942  Read(5,*) Cutint
*     If (KWord(1:1).eq.'*') Go To 942
*     If (KWord.eq.BLine)    Go To 942
*     Read(KWord,*,Err=988) CutInt
      CutInt = Abs(CutInt)
      Go To 998
*                                                                      *
****** MEMO ************************************************************
*                                                                      *
*     Screen off memory
*
*951  Read(5,'(A)',Err=988) KWord
*     If (KWord(1:1).eq.'*') Go To 951
*     If (KWord.eq.BLine)    Go To 951
*     Read(KWord,*,Err=988) MemHid
*     If (MemHid.le.0) MemHid = 1
*     Go To 998
*                                                                      *
****** NOIN ************************************************************
*                                                                      *
*     Disable the utilization of translational and
*     rotational invariance of the energy in the
*     computation of the molecular gradient.
*
*953  TRSymm=.False.
*     Go To 998
*                                                                      *
****** SELE ************************************************************
*                                                                      *
*     selection option
*
 960  Continue
      slct=.true.
      Call lCopy(mxpert,[.false.],0,lPert,1)
*962  Continue
      Read(5,*) nslct
*     If (KWord(1:1).eq.'*') Go To 962
*     If (KWord.eq.BLine)    Go To 962
*     Read(KWord,*) nSlct
*
      Read(5,*) (iTemp(iElem),iElem=1,nSlct)
      Do 964 iElem=1,nSlct
         lpert(iTemp(iElem)) = .True.
964   Continue
      Go To 998
*                                                                      *
****** REMO ************************************************************
*                                                                      *
 260  Continue
      Slct=.true.
      Read(5,*) nslct
*
      Read(5,*) (iTemp(iElem),iElem=1,nSlct)
      Do 264 iElem=1,nSlct
         lpert(iTemp(iElem)) = .false.
264   Continue
      Go To 998
*                                                                      *
****** NOSC ************************************************************
*                                                                      *
*     Change default for the prescreening.
*
 965  PreScr  = .False.
      Go To 998
*                                                                      *
****** ONEO ************************************************************
*                                                                      *
*     Do not compute two electron integrals.
*
 990  Onenly = .TRUE.
      Go To 998
*                                                                      *
****** TEST ************************************************************
*                                                                      *
*     Process only the input.
*
 991  Test = .TRUE.
      Go To 998
*                                                                      *
****** SHOW ************************************************************
*                                                                      *
*-----Raise the printlevel to show gradient contributions
*
 992  Continue
      Show=.true.
      Go To 998
*                                                                      *
****** EXTR ************************************************************
*                                                                      *
*     Put the program name and the time stamp onto the extract file
*
971   Write (6,*) 'InputH: EXTRACT option is redundant and is ignored!'
      Go To 998
*                                                                      *
****** VERB ************************************************************
*                                                                      *
*     Verbose output
*
 912  nPrint( 1)=6
      nPrint(99)=6
      Go To 998
*                                                                      *
****** NONA ************************************************************
*                                                                      *
*     Compute the anti-symmetric overlap gradient only.
*
972   Nona=.true.
      Run_MCLR=.False.
      Go To 998
*                                                                      *
****** NOMC ************************************************************
*                                                                      *
*     Request no automatic run of MCLR
*
973   Run_MCLR=.False.
      Go To 998
************************************************************************
*                                                                      *
*                          End of input section.                       *
*                                                                      *
************************************************************************
 997  Continue
*
      iPrint=nPrint(iRout)
*
      iOpt = 1
      if (onenly) iopt=0
      iRC = -1
      Lu_Mck=35
      Call OpnMck(irc,iOpt,'MCKINT',Lu_Mck)
      If (iRC.ne.0) Then
         Write (6,*) 'InputH: Error opening MCKINT'
         Call QTrace
         Call Abend()
      End If
      If (ipert.eq.1) Then
        Label2='Geometry'
        LabelOp='PERT    '
        irc=-1
        iopt=0
        Call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
        sIrrep=.true.
        iCntrl=1
      Else If (ipert.eq.2) Then
        Label2='Hessian'
        LabelOp='PERT    '
        Call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
        iCntrl=1
      Else If (ipert.eq.3) Then
        LabelOp='PERT    '
        Label2='Magnetic'
        Call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
        Write (6,*) 'InputH: Illegal perturbation option'
        Write (6,*) 'iPert=',iPert
        Call QTrace
        Call Abend()
      Else If (ipert.eq.4) Then
        LabelOp='PERT    '
        Label2='Relativistic'
        Call cWrMck(iRC,iOpt,LabelOp,1,Label2,iDummer)
        Write (6,*) 'InputH: Illegal perturbation option'
        Write (6,*) 'iPert=',iPert
        Call QTrace
        Call Abend()
      Else
        Write (6,*) 'InputH: Illegal perturbation option'
        Write (6,*) 'iPert=',iPert
        Call QTrace
        Call Abend()
      End If

*     If (lEq)  TRSymm=.False.
*     If (Slct) TRSymm=.False.
*
      mDisp = 0
      mdc = 0
      Do 10 iCnttp = 1, nCnttp
         Do 20 iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            mDisp = mDisp + 3*(nIrrep/nStab(mdc))
 20      Continue
 10   Continue
*
      Write (6,*)
      Write (6,'(20X,A,E10.3)')
     &  ' Threshold for contributions to the gradient or Hessian:',
     &   CutInt
      Write (6,*)
*
      If (Nona) Then
         Write (6,*)
         Write (6,'(20X,A)')
     &   ' McKinley only is computing the antisymmetric gradient '//
     &   ' of the overlap integrals for the NonAdiabatic Coupling.'
         Write (6,*)
      End If
*
      If (iCntrl.eq.1) Then
*
*
*     Generate symmetry adapted cartesian displacements
*
      If (iPrint.ge.6) Then
      Write (6,*)
      Write (6,'(20X,A)')
     &           '********************************************'
      Write (6,'(20X,A)')
     &           '* Symmetry Adapted Cartesian Displacements *'
      Write (6,'(20X,A)')
     &           '********************************************'
      Write (6,*)
      End If
      Call ICopy(mxdc*8,[0],0,IndDsp,1)
      Call ICopy(mxdc*3,[0],0,InxDsp,1)
      Call GetMem('ATDISP','ALLO','INTE',ipad,mdisp)
      Call GetMem('DEGDISP','ALLO','INTE',ipdd,mdisp)
      nDisp = 0
      Do 100 iIrrep = 0, nIrrep-1
         lDisp(iIrrep) = 0
         Type = .True.
*        Loop over basis function definitions
         mdc = 0
         mc = 1
         Do 110 iCnttp = 1, nCnttp
*           Loop over unique centers associated with this basis set.
            Do 120 iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               IndDsp(mdc,iIrrep) = nDisp
*              Loop over the cartesian components
               Do 130 iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                nIrrep/nStab(mdc),iChTbl,iIrrep,
     &                iComp,nStab(mdc)) ) Then
                      nDisp = nDisp + 1
                      If (nDisp.gt.mDisp) Then
                         Write (6,*) 'nDisp.gt.mDisp'
                         Call Abend
                      End If
                      If (iIrrep.eq.0) InxDsp(mdc,iCar+1) = nDisp
                      lDisp(iIrrep) = lDisp(iIrrep) + 1
                      If (Type) Then
                         If (iPrint.ge.6) Then
                         Write (6,*)
                         Write (6,'(10X,A,A)')
     &                    ' Irreducible representation : ',
     &                      lIrrep(iIrrep)
                         Write (6,'(10X,2A)')
     &                      ' Basis function(s) of irrep: ',
     &                       lBsFnc(iIrrep)
                         Write (6,*)
                         Write (6,'(A)')
     &                   ' Basis Label        Type   Center Phase'
                         End If
                         Type = .False.
                      End If
                      If (iPrint.ge.6)
     &                Write (6,'(I4,3X,A8,5X,A1,7X,8(I3,4X,I2,4X))')
     &                      nDisp,LblCnt(mdc),xyz(iCar),
     &                      (mc+iCo,iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                      iOper,nIrrep),iComp)*
     &                      iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                      iOper,nIrrep)),
     &                      iCo=0,nIrrep/nStab(mdc)-1 )
                      Write (ChDisp(nDisp),'(A,1X,A1)')
     &                       LblCnt(mdc),xyz(iCar)
                      iwork(ipad+ndisp-1)=icnttp
                      iwork(ipdd+ndisp-1)=nIrrep/nstab(mdc)
                  End If
*
 130           Continue
               mc = mc + nIrrep/nStab(mdc)
 120        Continue
 110     Continue
*
 100  Continue
*
      If (nDisp.ne.mDisp) Then
         Write (6,*) 'InputH: nDisp.ne.mDisp'
         Write (6,*) 'nDisp,mDisp=',nDisp,mDisp
         Call QTrace
         Call Abend()
      End If
      If (sIrrep) Then
         ndisp=ldisp(0)
         Do i= 1,nIrrep-1
            lDisp(i)=0
         End Do
      End If
      Call GetMem('TDISP','Allo','INTE',ipTD,ndisp)
      Call ICopy(nDisp,[30],0,iWork(ipTD),1)
      iOpt = 0
      iRC = -1
      labelOp='ndisp   '
      Call iWrMck(iRC,iOpt,labelop,1,[ndisp],iDummer)
      If (iRC.ne.0) Then
         Write (6,*) 'InputH: Error writing to MCKINT'
         Write (6,'(A,A)') 'labelOp=',labelOp
         Call QTrace
         Call Abend()
      End If
      LABEL='DEGDISP'
      iRc=-1
      iOpt=0
      Call iWrMck(iRC,iOpt,Label,idum,iwork(ipdd),idum)
      If (iRC.ne.0) Then
         Write (6,*) 'InputH: Error writing to MCKINT'
         Write (6,'(A,A)') 'LABEL=',LABEL
         Call QTrace
         Call Abend()
      End If
      Call GetMem('DEGDISP','Free','INTE',ipdd,mDisp)
      LABEL='NRCTDISP'
      iRc=-1
      iOpt=0
      Call iWrMck(iRC,iOpt,Label,idum,iwork(ipad),idum)
      If (iRC.ne.0) Then
         Write (6,*) 'InputH: Error writing to MCKINT'
         Write (6,'(A,A)') 'LABEL=',LABEL
         Call QTrace
         Call Abend()
      End If
      Call GetMem('ATDISP','Free','INTE',ipad,mDisp)
      LABEL='TDISP'
      iRc=-1
      iOpt=0
      Call iWrMck(iRC,iOpt,Label,idum,iWork(ipTD),idum)
      If (iRC.ne.0) Then
         Write (6,*) 'InputH: Error writing to MCKINT'
         Write (6,'(A,A)') 'LABEL=',LABEL
         Call QTrace
         Call Abend()
      End If
      Call GetMem('TDISP','FREE','INTE',ipTD,ndisp)
*
      Else If (iCntrl.eq.2) Then
          Write(6,*) 'Svaret aer 48 '
      Else If (iCntrl.eq.3) Then
          Write(6,*) 'Svaret aer 48'
      End If
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
*     Set up data for the utilization of the translational
*     and rotational invariance of the energy.
*
      If (TRSymm) Then
         iSym(1) = 0
         iSym(2) = 0
         iSym(3) = 0
         Do 15 i = 1, Min(nIrrep-1,5)
            j = i
            If (i.eq.3) j = 4
            Do 16 k = 1, 3
               If (iAnd(iOper(j),2**(k-1)).ne.0) iSym(k) = 2**(k-1)
   16       Continue
   15    Continue
         nTR = 0
*--------Translational equations
         Do 150 i = 1, 3
            If (iSym(i).eq.0) nTR = nTR + 1
 150     Continue
         If (iPrint.ge.99) Write (6,*) ' nTR=',nTR
*--------Rotational equations
         Do 160 i = 1,3
            j = i+1
            If (j.gt.3) j = j-3
            k = i+2
            If (k.gt.3) k = k-3
            ijSym = iEor(iSym(j),iSym(k))
            If (ijSym.eq.0) nTR = nTR + 1
 160     Continue
         If (nTR.eq.0) Then
            TRSymm = .False.
            Go To 9876
         End If
         If (iPrint.ge.99) Write (6,*) ' nTR=',nTR
         Call GetMem('Amtrx','Allo','Real',ipAm,lDisp(0)**2)
         Call GetMem('Temp ','Allo','Real',ipTmp,nTR**2)
         Call GetMem('Coor ','Allo','Real',ipC,lDisp(0)*4)
         Call GetMem('Car  ','Allo','Inte',ipCar,lDisp(0))
*
         call dcopy_(nTR*lDisp(0),[Zero],0,Work(ipAm),1)
         call dcopy_(4*lDisp(0),[Zero],0,Work(ipC),1)
*
*        Generate temporary information of the symmetrical
*        displacements.
*
         ldsp = 0
         mdc = 0
         iIrrep = 0
         Do 2100 iCnttp = 1, nCnttp
            jxyz = dbsc(iCnttp)%ipCntr
            Do 2200 iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
*              Call RecPrt(' Coordinates',' ',Work(jxyz),1,3)
               Fact = Zero
               iComp = 0
               If (Work(jxyz  ).ne.Zero) iComp = iOr(iComp,1)
               If (Work(jxyz+1).ne.Zero) iComp = iOr(iComp,2)
               If (Work(jxyz+2).ne.Zero) iComp = iOr(iComp,4)
               Do 2250 jIrrep = 0, nIrrep-1
                  If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                  nIrrep/nStab(mdc),iChTbl,jIrrep,
     &                  iComp,nStab(mdc)) ) Then
                     Fact = Fact + One
                  End If
 2250          Continue
               Do 2300 iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                  nIrrep/nStab(mdc),iChTbl,iIrrep,
     &                  iComp,nStab(mdc)) ) Then
                     ldsp = ldsp + 1
*--------------------Transfer the coordinates
                     ip = 4*(ldsp-1) + ipC
                     call dcopy_(3,Work(jxyz),1,Work(ip),1)
*--------------------Transfer the multiplicity factor
                     Work(ip+3) = Fact
                     iWork(ipCar-1+ldsp) = iCar + 1
                  End If
 2300          Continue
               jxyz = jxyz + 3
 2200       Continue
 2100    Continue
         If (iPrint.ge.99) Then
            Call RecPrt(' Information',' ',Work(ipC),4,lDisp(0))
            Write (6,*) (iWork(i),i=ipCar,ipCar+lDisp(0)-1)
         End If
*
*--------Set up coefficient for the translational equations
*
         iTR = 0
         Do 1110 i = 1,3
            If (iSym(i).ne.0) Go To 1110
            iTR = iTR + 1
            Do 1120 ldsp = 1, lDisp(0)
               If (iWork(ipCar+ldsp-1).eq.i) Then
                  ipOut= ipC + 4*(ldsp-1) + 3
                  ipIn = nTR*(ldsp-1) + iTR + ipAm - 1
                  Work(ipIn) = Work(ipOut)
               End If
 1120       Continue
 1110    Continue
*
*--------Set up coefficient for the rotational invariance
*
         Do 1210 i = 1, 3
            j = i + 1
            If (j.gt.3) j = j - 3
            k = i + 2
            If (k.gt.3) k = k - 3
            ijSym = iEor(iSym(j),iSym(k))
            If (ijSym.ne.0) Go To 1210
            iTR = iTR + 1
            Do 1220 ldsp = 1, lDisp(0)
               ipIn = nTR*(ldsp-1) + iTR + ipAm - 1
               ipOut = ipC + 4*(ldsp-1)
               If (iWork(ipCar+ldsp-1).eq.j) Then
                  Fact = Work(ipOut+3) * Work(ipOut+k-1)
                  Work(ipIn) = Fact
               Else If (iWork(ipCar+ldsp-1).eq.k) Then
                  Fact = Work(ipOut+3) * Work(ipOut+j-1)
                  Work(ipIn) = -Fact
               End If
 1220       Continue
 1210    Continue
         If (iPrint.ge.99)
     &      Call RecPrt(' The A matrix',' ',Work(ipAm),nTR,lDisp(0))
*
*--------Now, transfer the coefficient of those gradients which will
*        not be computed directly.
*        The matrix to compute the inverse of is determined via
*        a Gram-Schmidt procedure.
*
*--------Pick up the other vectors
         Do 1230 iTR = 1, nTR
*           Write (*,*) ' Looking for vector #',iTR
            ovlp = Zero
            kTR = 0
*-----------Check all the remaining vectors
            Do 1231 ldsp = 1, lDisp(0)
               Do 1235 jTR = 1, iTR-1
                  If (iTemp(jTR).eq.ldsp) Go To 1231
 1235          Continue
*              Write (*,*) ' Checking vector #', ldsp
               ipNew = ipAm + nTR*(ldsp-1)
               ipIn = ipTmp + nTR*(iTR-1)
               call dcopy_(nTR,Work(ipNew),1,Work(ipIn),1)
*              Call RecPrt(' Vector',' ',Work(ipIn),nTR,1)
*--------------Gram-Schmidt orthonormalize against accepted vectors
               Do 1232 lTR = 1, iTR-1
                  ipOld = ipTmp + nTR*(lTR-1)
                  alpha = DDot_(nTR,Work(ipIn),1,Work(ipOld),1)
*                 Write (*,*) ' <x|y> =', alpha
                  Call DaXpY_(nTR,-alpha,Work(ipOld),1,Work(ipIn),1)
 1232          Continue
*              Call RecPrt(' Remainings',' ',Work(ipIn),nTR,1)
               alpha = DDot_(nTR,Work(ipIn),1,Work(ipIn),1)
*              Write (*,*) ' Remaining overlap =', alpha
*--------------Check the remaining magnitude of vector after Gram-Schmidt
               If (alpha.gt.ovlp) Then
                  kTR = ldsp
                  ovlp = alpha
               End If
 1231       Continue
            If (kTR.eq.0) Then
               Write (6,*) ' No Vector found!'
               Call Abend
            End If
*           Write (*,*) ' Selecting vector #', kTR
*-----------Pick up the "best" vector
            ipNew = ipAm + nTR*(kTR-1)
            ipIn = ipTmp + nTR*(iTR-1)
            call dcopy_(nTR,Work(ipNew),1,Work(ipIn),1)
            Do 1233 lTR = 1, iTR-1
               ipOld = ipTmp + nTR*(lTR-1)
               alpha = DDot_(nTR,Work(ipIn),1,Work(ipOld),1)
               Call DaXpY_(nTR,-alpha,Work(ipOld),1,Work(ipIn),1)
 1233       Continue
            alpha = DDot_(nTR,Work(ipIn),1,Work(ipIn),1)
            Call DScal_(nTR,One/Sqrt(alpha),Work(ipIn),1)
            iTemp(iTR) = kTR
 1230    Continue
         Do 1234 iTR = 1, nTR
            ipNew = ipAm + nTR*(iTemp(iTR)-1)
            ipIn  = ipTmp + nTR*(iTR-1)
            call dcopy_(nTR,Work(ipNew),1,Work(ipIn),1)
            call dcopy_(nTR,[Zero],0,Work(ipNew),1)
 1234    Continue
         If (iPrint.ge.99) Then
            Call RecPrt(' The A matrix',' ',Work(ipAm),nTR,lDisp(0))
            Call RecPrt(' The T matrix',' ',Work(ipTmp),nTR,nTR)
            Write (6,*) (iTemp(iTR),iTR=1,nTR)
         End If
*
*        Compute the inverse of the T matrix
*
         Call MatInvert(Work(ipTmp),nTR)
         If (IPrint.ge.99)
     &      Call RecPrt(' The T-1 matrix',' ',Work(ipTmp),nTR,nTR)
         Call DScal_(nTR**2,-One,Work(ipTmp),1)
*
*        Generate the complete matrix
*
         Call GetMem(' Temp2','Allo','Real',ipScr,nTR*lDisp(0))
         Call DGEMM_('N','N',
     &               nTR,lDisp(0),nTR,
     &               1.0d0,Work(ipTmp),nTR,
     &               Work(ipAm),nTR,
     &               0.0d0,Work(ipScr),nTR)
         If (IPrint.ge.99)
     &      Call RecPrt(' A-1*A',' ',Work(ipScr),nTR,lDisp(0))
         call dcopy_(lDisp(0)**2,[Zero],0,Work(ipAm),1)
         call dcopy_(lDisp(0),[One],0,Work(ipAm),lDisp(0)+1)
         Do 1250 iTR = 1, nTR
            ldsp = iTemp(iTR)
            ipOut = ipScr + iTR - 1
            ipIn  = ipAM  + ldsp - 1
            call dcopy_(lDisp(0),Work(ipOut),nTR,Work(ipIn),lDisp(0))
 1250    Continue
         If (iPrint.ge.99)
     &      Call RecPrt('Final A matrix',' ',
     &                  Work(ipAm),lDisp(0),lDisp(0))
*
*
         Call GetMem(' Temp2','Free','Real',ipScr,nTR*lDisp(0))
         Call GetMem('Car  ','Free','Inte',ipCar,lDisp(0))
         Call GetMem('Coor ','Free','Real',ipC,lDisp(0)*4)
         Call GetMem('Temp ','Free','Real',ipTmp,nTR**2)
         Do 1501 iTR = 1, nTR
            ldsp = iTemp(iTR)
            LPert(ldsp)=.False.
 1501    Continue
*
         Write (6,*)
         Write (6,'(20X,A,A)')
     &      ' Automatic utilization of translational and',
     &      ' rotational invariance of the energy is employed.'
         Write (6,*)
         Do 7000 i = 1, lDisp(0)
            If (lpert(i)) Then
               Write (6,'(25X,A,A)') Chdisp(i), ' is independent'
            Else
               Write (6,'(25X,A,A)') Chdisp(i), ' is dependent'
            End If
 7000    Continue
         Write (6,*)
*
      Else
         nTR = 0
         If (iPrint.ge.6) Then
         Write (6,*)
         Write (6,'(20X,A,A)')
     &      ' No automatic utilization of translational and',
     &      ' rotational invariance of the energy is employed.'
         Write (6,*)
         End If
      End If
*
      If (Slct) Then
         Write (6,*)
         Write (6,'(20X,A)') ' The Selection option is used'
         Write (6,*)
         Do 7100 i = 1, lDisp(0)
            If (lpert(i)) Then
               Write (6,'(25X,A,A)') Chdisp(i), ' is computed'
            Else
               Write (6,'(25X,A,A)') Chdisp(i), ' is set to zero'
            End If
 7100    Continue
         Write (6,*)
      End If
*
 9876 Continue
      Call Datimx(KWord)
      goto 888
        Call qTrace
        Write(6,*) ' *** Error in subroutine INPUTG ***'
        Write(6,*) '     Abend in subroutine WrOne'
        Call Abend

 888  Continue
      Call ICopy(nIrrep,[0],0,nFck,1)
      Do iIrrep=0,nIrrep-1
        If (iIrrep.ne.0) Then
          Do jIrrep=0,nIrrep-1
           kIrrep=NrOpr(iEOR(ioper(jIrrep),ioper(iIrrep)),
     &                  iOper,nIrrep)
           If (kIrrep.lt.jIrrep)
     &     nFck(iIrrep)=nFck(iIrrep)+nBas(jIrrep)*nBas(kIrrep)
          End Do
        Else
           Do jIrrep=0,nIrrep-1
             nFck(0)=nFck(0)+nBas(jIrrep)*(nBas(jIrrep)+1)/2
           End Do
        End If
      End Do
*
      Call QExit('InputH')
      Return
      End
