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
* Copyright (C) Yannick Carissan                                       *
*               2005, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Localisation(iReturn)
c
c     Author: Yannick Carissan
c
c     Modifications:
c        - October 10, 2005 (Thomas Bondo Pedersen):
c          completely restructed; introduce Boys and Cholesky
c          localisations.
c        - December 2005 / January 2006 (Thomas Bondo Pedersen):
c          Edmiston-Ruedenberg, PAO, and pair domain analysis included.
c
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "inflocal.fh"
#include "debug.fh"
#include "real.fh"
      Character*80 Title, Txt
      Character NameFile*20, Filename*6
      Character*180 Line, Get_Ln
      External Get_Ln

cvv      Integer  LocUtil_Models
cvv      External LocUtil_Models

      Character*7  matname
      Character*2  PreFix
      Character*4  Model
      Character*12 SecNam
      Character*15 AddInfoString
      Parameter (SecNam = 'Localisation')

      Real*8 ERFun(2), xnr0(8)
      Integer IndT(7,8)

      Call qEnter(SecNam)

C     Start timing.
C     -------------

      Call CWTime(C1,W1)

C     Dummy allocation used to flush memory at the end.
C     -------------------------------------------------

      l_Dum = 1
      Call GetMem('LocDum','Allo','Real',ip_Dum,l_Dum)

C     Print banner.
C     -------------

c     Call Banner_Localisation()

C     Read basic info from runfile and INPORB.
C     ----------------------------------------

C     Quick and dirty read of the FileOrb name before INPORB is opened
      LuSpool=17
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'LOCALISATION')
      LC_FileOrb = ' '
100   Line = Get_Ln(LuSpool)
      Call UpCase(Line)
      If (Line(1:4).eq.'FILE') Goto 200
      If (Line(1:4).eq.'END ') Goto 999
      Goto 100
200   Line = Get_Ln(LuSpool)
      Call FileOrb(Line,LC_FileOrb)
999   Call Close_LuSpool(LuSpool)

      Call GetInfo_Localisation_0()

C     Read and process input.
C     -----------------------

      Call ReadInp_Localisation()

C     Check that all is OK.
C     ---------------------

      irc=0
      Call Chk_Input(irc)
      If (irc .gt. 0) Then
         Call SysAbendMsg(SecNam,'Inconsistent input',' ')
      Else If (irc .lt. 0) Then
         iReturn = 0
         Go To 1 ! nothing to do; exit...
      End If

C     If test option was specified, we need to keep a copy of the
C     original MOs.
C     -----------------------------------------------------------

      If (Test_Localisation.or.Analysis.or.AnaPAO.or.AnaPAO_Save
     &    .or.LocNatOrb .or.LocCanOrb) Then
         lMOrig = nCMO
         Call GetMem('MOrig','Allo','Real',ipMOrig,lMOrig)
         Call dCopy_(lMOrig,Work(ipCMO),1,Work(ipMOrig),1)
      Else
         ipMOrig = -99999999
         lMOrig  = 0
      End If

C     Print initial orbitals.
C     -----------------------

      If (.not.Silent .and. PrintMOs) Then
         Write(Title,'(80x)')
         Write(Title,'(a)') 'Initial MO''s'
         iPrWay=2 ! short
         Call PriMO_Localisation(Title,.True.,.True.,
     &                           -1.0d0,1.0d5,
     &                           nSym,nBas,nOrb,Name,Work(ipEor),
     &                           Work(ipOcc),Work(ipCMO),
     &                           iPrWay,iWork(ipInd))
      End If

      If (Debug) Then
         Write(6,'(A,A,I2)') SecNam,': debug info at start:'
         Write(6,'(A,8I8)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
         Write(6,'(A,8I8)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
         Write(6,'(A,8I8)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
         Write(6,'(A,8I8)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
      End If

C     Evaluate ER functional for initial orbitals.
C     --------------------------------------------

      If (EvalER) Then
         ERFun(1) = 0.0d0
         Call ComputeFuncER(ERFun(1),Work(ipCMO),nBas,nOrb2Loc,nFro,
     &                      nSym,Timing)
      End If

C     Localise orbitals (if the user did not request us to skip it).
C     --------------------------------------------------------------

      If (Skip) Then
         Write(6,'(//,5X,A,//)')
     &   'NOTICE: LOCALISATION PROCEDURE SKIPPED BY USER REQUEST!'
         AddInfoString='SKIPPED LOCALI '
         AddInfoVal=0.0d0
         iTol=15
      Else
         If (Timing) Call CWTime(C1_Loc,W1_Loc)
         If (LocModel.eq.1 .or. LocModel.eq.2 .or. LocModel.eq.4) Then
            If (LocModel .eq. 1) Then
               Model='Pipe'
               AddInfoString='PIPEKFUNCTIONAL'
            Else If (LocModel .eq. 2) Then
               Model='Boys'
               AddInfoString='BOYSFUNCTIONAL '
            Else If (LocModel .eq. 4) Then
               Model='Edmi'
               AddInfoString='ERFUNCTIONAL   '
            End If
            irc = 0
            Call Localise_Iterative(irc,Model,Functional)
            If (irc .ne. 0) Then
               Write(Txt,'(A,I3)')
     &         'Return code from Localise_Iterative:',irc
               Call SysAbendMsg(SecNam,'Localisation failed!',Txt)
            End If
            AddInfoVal=Functional
            iTol = 4
         Else If (LocModel .eq. 3) Then
            If (LocPAO) Then
               Model='PAO '
               AddInfoString='CHOLESKY PAO   '
            Else
               Model='Chol'
               AddInfoString='CHOLESKY NORM  '
            End If
            irc = 0
            Call Localise_Noniterative(irc,Model,xNrm)
            If (irc .ne. 0) Then
               Write(Txt,'(A,I3)')
     &         'Return code from Localise_Noniterative:',irc
               Call SysAbendMsg(SecNam,'Localisation failed!',Txt)
            End If
            AddInfoVal=xNrm
            iTol = 4
         Else If (Wave) Then ! wavelet transform
            irc = 0
            Call Wavelet_Transform(irc,ipCMO,nSym,nBas,nFro,nOrb2Loc,
     &                                 iWave,.false.,xNrm)
            If (irc .ne. 0) Then
               Write(Txt,'(A,I3)')
     &         'Return code from Wavelet_Transform:',irc
               Call SysAbendMsg(SecNam,'Localisation failed!',Txt)
            End If
            AddInfoVal=xNrm
            iTol = 4
         Else If (DoCNOs) Then ! Constrained NOs (analysis)
            irc = 0
            Call Get_CNOs(irc,nFro,nOrb2Loc,xNrm)
            If (irc .ne. 0) Then
               Write(Txt,'(A,I3)')
     &         'Return code from Get_CNOs:',irc
               Call SysAbendMsg(SecNam,'Localisation failed!',Txt)
            End If
            AddInfoVal=xNrm
            iTol = 4
         Else
            Call SysAbendMsg(SecNam,'Logical bug','(LocModel)')
            AddInfoString='?!?!?!?!?!?!?!?'
            AddInfoVal=-9.9d15
            iTol = 15
         End If
         If (Timing) Then
            Call CWTime(C2_Loc,W2_Loc)
            Write(6,'(/,1X,A,F10.2,A)')
     &      'CPU  time for localisation procedure:',
     &      C2_Loc-C1_Loc,' seconds'
            Write(6,'(1X,A,F10.2,A,/)')
     &      'Wall time for localisation procedure:',
     &      W2_Loc-W1_Loc,' seconds'
         End If
      End If

C     Info for check system.
C     ----------------------

      Call Add_Info(AddInfoString,AddInfoVal,1,iTol)

C     Order local orbitals according to Cholesky ordering.
C     ----------------------------------------------------

      If (Order) Then
         Write(6,'(/,1X,A,A)')
     &   'Sorting local orbitals according to Cholesky ordering.',
     &   ' (Based on overlap U=X^TSC.)'
         Call Sort_Localisation(Work(ipCMO),nBas,nOrb2Loc,nFro,nSym)
      End If

C     Evaluate ER functional for local orbitals.
C     ------------------------------------------

      If (EvalER) Then
         ERFun(2) = 0.0d0
         Call ComputeFuncER(ERFun(2),Work(ipCMO),nBas,nOrb2Loc,nFro,
     &                      nSym,Timing)
         Write(6,'(/,1X,A,1P,D15.8,/,1X,A,D15.8,/)')
     &    'ER functional for initial orbitals: ',ERFun(1),
     &    'ER functional for local   orbitals: ',ERFun(2)
      End If

C     Test section.
C     -------------

      If (Test_Localisation) Then
         Write(6,'(//,1X,A)')
     &   'Testing orbital localisation...'
         irc = -1
         Call TestLoc(irc)
         If (irc .eq. 0) Then
            Write(6,*) '...OK!'
         Else
            Write(6,*) SecNam,': localisation error detected!'
            Write(6,*) ' TestLoc returned ',irc
            iReturn = 1
            Go To 1 ! exit
         End If
      End If

C     Analysis.
C     ---------

      If (Analysis) Then
         PreFix = 'L_'
         If (AnaAtom) Then
            Call BitMap_Localisation_Atom(PreFix)
         Else
            Call BitMap_Localisation(PreFix)
         End If
      End If

C     Orbital domains.
C     ----------------

      If (DoDomain) Then
         irc = 0
         Call Domain_Localisation(irc)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': Domain error!'
            Write(6,*) ' Domain_Localisation returned ',irc
            Write(6,*) ' Program continues nevertheless...'
         End If
      End If

C     Local Natural Orbitals
C     ----------------------

      If (LocNatOrb .or. LocCanOrb) Then
         nbo=0
         Do iSym=1,nSym
            nbo=nbo+nBas(iSym)*nOrb2Loc(iSym)
         End Do
         Call GetMem('XCMO','Allo','Real',jCMO,2*nbo)
         jXMO=jCMO+nbo
*
         ibo=0
         jbo=0
         Do iSym=1,nSym
            kCMO=ibo+nBas(iSym)*nFro(iSym)
            call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(ipMOrig+kCMO),1,
     &                                           Work(jbo+jCMO),1)
            call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(ipCMO+kCMO),1,
     &                                           Work(jbo+jXMO),1)
            ibo=ibo+nBas(iSym)**2
            jbo=jbo+nBas(iSym)*nOrb2Loc(iSym)
         End Do
*
         If (LocNatOrb) Then
            matname='Density'
            jXarray=ipOcc
         Else
            matname='Fock'
            jXarray=ipEor
         EndIf
         lOff=0
         Do iSym=1,nSym
            xnr0(iSym) = ddot_(nOrb2Loc(iSym),1.0d0,0,
     &                             Work(jXarray+lOff+nFro(iSym)),1)
            lOff=lOff+nBas(iSym)
         End Do
*
         Call Loc_Nat_Orb(irc,Work(jCMO),Work(jXMO),Work(jXarray),
     &                        nOrb2Loc)
         If (irc.ne.0) Then
            Write(6,*) SecNam,': localisation error detected!'
            Write(6,*) ' Loc_Nat_Orb returned ',irc
            iReturn = 1
            Go To 1 ! exit
         EndIf
*
         write(6,*)
         write(6,*)' ------------------------------------------------- '
         write(6,*)' Sum of the eigenvalues of the partial ',matname,' '
         write(6,*)' matrix of the orbitals that have been localised   '
         write(6,*)' ------------------------------------------------- '
         write(6,*)'    Symm.        before     / after localisation   '
         write(6,*)' ------------------------------------------------- '
         lOff=0
         Do iSym=1,nSym
            xnr1 = ddot_(nOrb2Loc(iSym),1.0d0,0,
     &                                 Work(jXarray+lOff+nFro(iSym)),1)
            lOff=lOff+nBas(iSym)
            write(6,'(3X,I4,8X,F11.5,4X,F11.5)') iSym, xnr0(iSym), xnr1
         End Do
         write(6,*)' ------------------------------------------------- '
         write(6,*)
*
         ibo=0
         jbo=0
         Do iSym=1,nSym
            kCMO=ibo+nBas(iSym)*nFro(iSym)
            call dcopy_(nBas(iSym)*nOrb2Loc(iSym),Work(jbo+jXMO),1,
     &                                           Work(ipCMO+kCMO),1)
            ibo=ibo+nBas(iSym)**2
            jbo=jbo+nBas(iSym)*nOrb2Loc(iSym)
         End Do
         Call GetMem('XCMO','Free','Real',jCMO,2*nbo)
      EndIf

C-TBP, July 2010 (in connection with fixing deleted orbitals bug, patch
C 7.7.073_Localisation):
C Moved the following section outside print section, which is only
C executed if PrintMOs=.True., implying that the info on LOCORB
C would differ depending on this print flag!!

C     Print a warning if localisation is done in a subset of orbitals
C     belonging to more than one subset (e.g. mixing occupied and
C     virtual orbitals and thus breaking the variational principle).
C     Set orbital energies to zero (unless localisation was skipped).
C     Also zero occupations if localisation mixed different subsets.
C     --------------------------------------------------------------

      iCheck=0
      iOff=0
      iSym=1
      jTyp=Min(6,Max(2,iWork(ipInd+nFro(1)))) ! Fro=Ina and Del=Vir
      Do While (iSym.le.nSym .and. iCheck.eq.0)
        jInd=ipInd+iOff+nFro(iSym)
        j=1
        Do while ( j.lt.nOrb2Loc(iSym) .and. iCheck.eq.0)
           iCheck=Min(6,Max(2,iWork(jInd+j)))-jTyp
           j=j+1
        End Do
        iOff=iOff+nBas(iSym)
        iSym=iSym+1
      End Do
      If (.not.LocCanOrb .and. .not.Skip) Then
         iOff=0
         Do iSym=1,nSym
            kEor=ipEor+iOff+nFro(iSym)
            Call FZero(Work(kEor),nOrb2Loc(iSym))
            iOff=iOff+nBas(iSym)
         End Do
      End If
      If (iCheck.eq.0) Then
         jPrt=0
         If (jTyp.gt.3 .and. jTyp.lt.6) jPrt=1
      Else
         write(6,*) '****  WARNING  ***  WARNING  ***  WARNING  ****'
         write(6,*) ' The orbitals for which localisation is reque'//
     &              'sted belong to more than one of the subspaces:'
         write(6,*) ' (Froz+Inac|RAS1|RAS2|RAS3|Sec+Del). '
         write(6,*) '***********************************************'
         jPrt=1
      End If
      If (.not.LocNatOrb .and. .not.Skip .and. .not.DoCNOs) Then
         If(jPrt.eq.1) Then
            iOff=0
            Do iSym=1,nSym
               kOcc=ipOcc+iOff+nFro(iSym)
               Call FZero(Work(kOcc),nOrb2Loc(iSym))
               iOff=iOff+nBas(iSym)
            End Do
         End If
      End If

C     Print final localised orbitals.
C     -------------------------------

      If (PrintMOs) Then
         Write(Title,'(80x)')
         Write(Title,'(a)') 'Final localised MO''s'
         iPrWay=2 ! short
         Call PriMO_Localisation(Title,.True.,.True.,
     &                          -1.0d0,1.0d5,
     &                          nSym,nBas,nOrb,Name,Work(ipEor),
     &                          Work(ipOcc),Work(ipCMO),
     &                          iPrWay,iWork(ipInd))
      End If

C     Write LOCORB file.
C     ------------------

      Write(Namefile,'(A)')  'LOCORB'
      Write(Title,'(80X)')
      Write(Title,'(A)') 'Localised orbitals'
      LU_=isFreeUnit(11)
      j=ipInd-1
      Call iZero(IndT,56)
      Do iSym=1,nSym
         Do k=1,nOrb(iSym)
            kIndT=iWork(j+k)
            If (kIndT.gt.0 .and. kIndT.le.7) Then
               IndT(kIndT,iSym)=IndT(kIndT,iSym)+1
            Else
               Call WarningMessage(2,
     &                             'Localisation: Illegal orbital type')
               Write(6,'(A,I6,A,I2,A,I9)')
     &         'Orbital',k,' of sym.',iSym,' has illegal type:',kIndT
               Call Abend()
            End If
         End Do
         j=j+nBas(iSym)
      End Do
      Call WrVec_Localisation(Namefile,LU_,'COEI',nSym,nBas,nBas,
     &                        Work(ipCMO),Work(ipOcc),Work(ipEor),IndT,
     &                        Title)
      If (.not.Silent) Then
         Write(6,'(1X,A)') 'The LOCORB file has been written.'
      End If

C     Write MOLDEN file.
C     ------------------

      iUHF=0
      Filename='MD_LOC'
      Call Molden_Interface(iUHF,Namefile,Filename,.False.)
      If (.not.Silent) Then
         Write(6,'(1X,A)') 'The MOLDEN file has been written.'
      End If

C     Set return code.
C     ----------------

      iReturn = 0

C     Free memory (by flushing).
C     --------------------------

    1 Call GetMem('LocDum','Flus','Real',ip_Dum,l_Dum)
      Call GetMem('LocDum','Free','Real',ip_Dum,l_Dum)

C     Print timing.
C     -------------

      If (.not.Silent) Then
         Call CWTime(C2,W2)
         CPUtot = C2 - C1
         WLLtot = W2 - W1
         Call Cho_CnvTim(CPUtot,icHour,icMin,cSec)
         Call Cho_CnvTim(WLLtot,iwHour,iwMin,wSec)
         Write(6,'(/,1X,A,I8,A,I2,A,F6.2,A)')
     &   '*** Total localisation time (CPU) : ',
     &   icHour,' hours ',icMin,' minutes ',cSec,' seconds ***'
         Write(6,'(1X,A,I8,A,I2,A,F6.2,A,/)')
     &   '*** Total localisation time (Wall): ',
     &   iwHour,' hours ',iwMin,' minutes ',wSec,' seconds ***'
      End If

C     That's it!
C     ----------

      Call qExit(SecNam)
      End
