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
      Subroutine RunTinker(natom,ipCord,ipMltp,ipIsMM,MltOrd,
     &                     DynExtPot,iQMChg,nAtMM,StandAlone,
     &                     DoDirect)
      Implicit Real*8 (a-h,o-z)
#include "para_info.fh"
*
#include "espf.fh"
*
      Character*180 Line
      Character*180 Get_Ln
      Character*12 ExtPotFormat
      Character*256 TkLine
      Logical DynExtPot,lFirst,StandAlone,DoDirect
      Integer RC
1000  Format('Molcas  ',i2,2x,i2)
1010  Format(3F15.8)
*
      Call qEnter('runtinker')
      iPL = iPL_espf()
      Write(ExtPotFormat,'(a4,i2,a6)') '(I4,',MxExtPotComp,'F13.8)'
*
* Always update the coordinates of the tinker xyz file
* WARNING: coordinates are converted to Angstroms
* This is done through a communication file: Project.qmmm
*
      lFirst = (ipMltp .eq. ip_Dummy)
*
* Only call Tinker on the master node
*
      ITkQMMM = 1
      If (MyRank .eq. 0) Then
        ITkQMMM = IsFreeUnit(ITkQMMM)
        Call Molcas_Open (ITkQMMM,'QMMM')
*
*     The MM subsystem can relax (microiterations, MD, ...) unless:
*     1) there are no QM multipoles
*     2) this is a call to retrieve MM energy/gradient/electrostatic potential only
*
        iRelax = 1
        If (lFirst .or. .not. StandAlone) iRelax = 0
        If (DoDirect) Then
          Write(ITkQMMM,1000) iRelax,-1
        Else
          Write(ITkQMMM,1000) iRelax,MltOrd/4
        End If
        Do iAtom = 1, natom
          Write(ITkQMMM,1010)(Work(ipCord+(iAtom-1)*3+K)*Angstrom,K=0,2)
        End Do
        If (.not.lFirst) Then
          Write(ITkQMMM,'(A)') 'Multipoles'
          If (iQMChg .eq. 0) Then
            If(iPL.ge.3) Write(6,'(A)') ' Multipoles passed to Tinker'
            iMlt = 0
            Do iAtom = 1, natom
               If (iWork(ipIsMM+iAtom-1).eq.0) Then
                  If (MltOrd.eq.1) Then
                     Write(ITkQMMM,'(I6,4F15.8)') iAtom,
     &                                  Work(ipMltp+iMlt),Zero,Zero,Zero
                  Else
                     Write(ITkQMMM,'(I6,4F15.8)') iAtom,
     &                                       (Work(ipMltp+iMlt+j),j=0,3)
                  End If
                  iMlt = iMlt + MltOrd
               Else
                  Write(ITkQMMM,'(I6,4F15.8)') iAtom,Zero,Zero,Zero,Zero
               End If
            End Do
          Else If (iQMChg .eq. 1) Then
            If (          (StandAlone.and.iPL.ge.2)
     &           .or.(.not.StandAlone.and.iPL.ge.3))
     &               Write(6,'(A)') ' ESPF multipoles passed to Tinker'
            iMlt = 0
            Do iAtom = 1, natom
               If (iWork(ipIsMM+iAtom-1).eq.0) Then
                  If (MltOrd.eq.1) Then
                     Write(ITkQMMM,'(I6,4F15.8)') iAtom,
     &                                  Work(ipMltp+iMlt),Zero,Zero,Zero
                  Else
                     Write(ITkQMMM,'(I6,4F15.8)') iAtom,
     &                                       (Work(ipMltp+iMlt+j),j=0,3)
                  End If
                  iMlt = iMlt + MltOrd
               Else
                  Write(ITkQMMM,'(I6,4F15.8)') iAtom,Zero,Zero,Zero,Zero
               End If
            End Do
          Else If (iQMChg.eq.2) Then
            If (          (StandAlone.and.iPL.ge.2)
     &           .or.(.not.StandAlone.and.iPL.ge.3))
     &              Write(6,'(A)') ' Mulliken charges passed to Tinker'
            Call GetMem('Mulliken','Allo','Real',ipMull,natom)
            Call Get_dArray('Mulliken Charge',Work(ipMull),natom)
            Do iAtom = 0, natom-1
               If (iWork(ipIsMM+iAtom).eq.0)
     &            Write(ITkQMMM,'(I6,4F15.8)') iAtom+1,
     &                                 Work(ipMull+iAtom),Zero,Zero,Zero
            End Do
            Call GetMem('Mulliken','Free','Real',ipMull,natom)
          Else If (iQMChg.eq.3) Then
            If (          (StandAlone.and.iPL.ge.2)
     &           .or.(.not.StandAlone.and.iPL.ge.3))
     &                Write(6,'(A)') ' LoProp charges passed to Tinker'
            Call GetMem('LoProp','Allo','Real',ipLPC,natom)
            Call Get_dArray('LoProp Charge',Work(ipLPC),natom)
            Do iAtom = 0, natom-1
               If (iWork(ipIsMM+iAtom).eq.0)
     &            Write(ITkQMMM,'(I6,4F15.8)') iAtom+1,
     &                                  Work(ipLPC+iAtom),Zero,Zero,Zero
            End Do
            Call GetMem('LoProp','Free','Real',ipLPC,natom)
          End If
        End If
        Close (ITkQMMM)
*
* Tinker is running
*
        Call Getenvf('Project ',Line)
        mLine = Len(Line)
        iLast = iCLast(Line,mLine)
        Line = Line(1:iLast)//'.xyz'
        Line = Line(1:iLast)//'.key'
        Line = '/tkr2qm_s ${Project}.xyz>${Project}.Tinker.log'
        Call Getenvf('TINKER ',TkLine)
        mLine = Len(TkLine)
        iLast = iCLast(TkLine,mLine)
        If (iLast.eq.0) Then
          Call Getenvf('MOLCAS',TkLine)
          mLine = Len(TkLine)
          iLast = iCLast(TkLine,mLine)
          TkLine = TkLine(1:iLast)//'/tinker/bin'
        End If
        iLast = iCLast(TkLine,mLine)
        nLine = Len(Line)
        jLast = iCLast(Line,nLine)
        Line = TkLine(1:iLast)//Line(1:jLast)
        Call StatusLine(' espf:',' Calling Tinker')
        RC=0
        Call Systemf(Line(1:iLast+jLast),RC)
      End If
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        Call GA_Sync()
        Call PFGET_ASCII('TINKER.LOG')
        Call PFGET_ASCII('QMMM')
        Call PFGET_ASCII('ESPF.EXTPOT')
        Call GA_Sync()
      End If
#endif
      If (           (StandAlone.and.iPL.ge.2)
     &     .or. (.not.StandAlone.and.iPL.ge.3)) Then
         iSomething = 0
         Lu=55
         Lu=IsFreeUnit(Lu)
         Call Molcas_Open(Lu,'TINKER.LOG')
 666     Continue
         Read(Lu,'(A)',End=667) Line
         iSomething = iSomething + 1
         nLine=LEN(Line)
         iLast = iCLast(Line,nLine)
         Write(6,*) Line(1:iLast)
         Go To 666
 667     Close (Lu)
         If (iSomething.eq.0) Then
            Write(6,*) ' Something bad with Tinker: no output !'
            Call Quit_OnUserError()
         End If
      End If
*
* Tinker post-processing
* WARNING: all Tinker results must be converted to atomic units !!!
* Convert the ESPF external potential and derivatives to something
* understandable by molcas, stored in the ESPF.EXTPOT file
*
      If (iPL.ge.2) Then
         Write(6,*) 'Back from Tinker'
         Write(6,*)
      End If
      Call GetMem('PotInESPF','Allo','Real',iESPF,natom*MxExtPotComp)
      Call DCopy_(MxExtPotComp*natom,Zero,0,Work(iESPF),1)
      ITkQMMM = IsFreeUnit(ITkQMMM)
      Call Molcas_Open (ITkQMMM,'QMMM')
      Line = Get_Ln(ITkQMMM)
      If (Index(Line,'MMisOK') .eq. 0) Then
         Write(6,*) 'Something wrong happend with Tinker'
         Call Abend()
      End If
      Do While (Index(Line,'TheEnd ') .eq. 0)
         Line=Get_Ln(ITkQMMM)
         If (Index(Line,'NMM ').ne.0) Then
            Call Get_I(2,nAtMM,1)
         Else If (Index(Line,'ESPF1 ').ne.0) Then
            Call Get_I(2,iAtom,1)
            Call Get_F(3,Work(iESPF+(iAtom-1)*MxExtPotComp),4)
         Else If (Index(Line,'ESPF21 ').ne.0) Then
            Call Get_I(2,iAtom,1)
            Call Get_F(3,Work(iESPF+(iAtom-1)*MxExtPotComp+4),3)
         Else If (Index(Line,'ESPF22 ').ne.0) Then
            Call Get_I(2,iAtom,1)
            Call Get_F(3,Work(iESPF+(iAtom-1)*MxExtPotComp+7),3)
         Else If (Index(Line,'FullCoupling ').ne.0) Then
            DynExtPot = .True.
         Else If (Index(Line,'MMq ').ne.0) Then
            Call Get_I(2,nMMq,1)
            Call GetMem('MMq_coord','Allo','Real',iMMqx,4*nMMq)
            Do iq = 1, nMMq
               Line=Get_Ln(ITkQMMM)
               Call Get_F(1,Work(iMMqx+(iq-1)*4),4)
            End Do
         End If
      End Do
      Close (ITkQMMM)
      ITkPot = IsFreeUnit(ITkQMMM)
      Call Molcas_Open (ITkPot,'ESPF.EXTPOT')
      If (DoDirect) Then
         Write(ITkPot,'(I10,1X,I2)') nMMq,0
         Do iq = 1, nMMq
            ibla = iMMqx+(iq-1)*4
            Work(ibla  ) = Work(ibla  ) * Angstrom
            Work(ibla+1) = Work(ibla+1) * Angstrom
            Work(ibla+2) = Work(ibla+2) * Angstrom
            Write(ITkPot,'(4F15.8)') (Work(ibla+J),J=0,3)
         End Do
         Call GetMem('MMq_coord','Free','Real',iMMqx,4*nMMq)
      Else
         Write(ITkPot,'(I1)') 0
         Do iAtom = 1, natom
            ibla = iESPF+(iAtom-1)*MxExtPotComp
            Work(ibla  ) = Work(ibla  ) * Angstrom
            Work(ibla+1) = Work(ibla+1) * Angstrom2
            Work(ibla+2) = Work(ibla+2) * Angstrom2
            Work(ibla+3) = Work(ibla+3) * Angstrom2
            Work(ibla+4) = Work(ibla+4) * Angstrom3
            Work(ibla+5) = Work(ibla+5) * Angstrom3
            Work(ibla+6) = Work(ibla+6) * Angstrom3
            Work(ibla+7) = Work(ibla+7) * Angstrom3
            Work(ibla+8) = Work(ibla+8) * Angstrom3
            Work(ibla+9) = Work(ibla+9) * Angstrom3
            Write(ITkPot,ExtPotFormat) iAtom,(Work(ibla+J),
     &                                      J=0,MxExtPotComp-1)
         End Do
      End If
      Close (ITkPot)
      Call GetMem('PotInESPF','Free','Real',iESPF,natom*MxExtPotComp)
*
      Call QExit('runtinker')
      Return
      End
