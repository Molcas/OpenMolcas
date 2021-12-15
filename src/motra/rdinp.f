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
* Copyright (C) 1991, Markus P. Fuelscher                              *
************************************************************************
      Subroutine RdInp_Motra
************************************************************************
*                                                                      *
*     Subroutine RdInp                                                 *
*                                                                      *
*     Purpose: Read input from inputstream                             *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)

      Parameter ( nCmd   = 16    )
      Parameter ( lCmd   = 4     )
      Character*4 CmdTab(nCmd)
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
      Character*180 Line,Blank
      Integer nDel2(8)
      character*(180) Get_Ln
      External Get_Ln
      Data CmdTab/'TITL','FROZ','DELE','PRIN','MOLO',
     &            'LUMO','JOBI','ONEL','FILE','AUTO',
     &            'EXTR','RFPE','CTON','DIAG','HDF5',
     &            'END '/
*
      COMMON / CHO_Minp / iCTonly, iDoInt
      Character*3 tv2disk
      COMMON / CHOTRAW / tv2disk
*


      iCTonly=0
      iDoInt =0
      ihdf5  =0
      tv2disk='PQK'
*----------------------------------------------------------------------*
*     Initialize some arrays                                           *
*----------------------------------------------------------------------*
      Do iSym=1,mxSym
        nDel(iSym)=0
        nOrb(iSym)=0
        CutThrs(iSym)=0.0D0
      End Do
      Call Get_iArray('Non valence orbitals',nFro,nSym)
      Do i=1,72
        Blank(i:i)=' '
      End Do
*----------------------------------------------------------------------*
*     Locate "start of input"                                          *
*----------------------------------------------------------------------*
      LuSpool = 17
      Call SpoolInp(LuSpool)
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'MOTRA')
*----------------------------------------------------------------------*
*     Read the input stream line by line and identify key command      *
*----------------------------------------------------------------------*
100   Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 100
      Call UpCase(Line)
110   jCmd=0
      Do iCmd=1,nCmd
        If ( Line(1:lCmd).eq.CmdTab(iCmd)(1:lCmd) ) jCmd=iCmd
      End Do
      If ( jCmd.eq.0 ) Then
         Write (6,*) 'RdInp: Unknown command at line: ', trim(Line)
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     Branch to the processing of the command sections                 *
*----------------------------------------------------------------------*
      Goto (1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
     &      1110,1120,1130,1140,1150,2000),jCmd
*---  Process the "TITLe" command -------------------------------------*
1010  nTit=0
15    Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 15
      Call UpCase(Line)
      if(nTit.gt.0) Then
       Do iCmd=1,nCmd
         If ( Line(1:lCmd).eq.CmdTab(iCmd)(1:lCmd) ) Goto 110
       End Do
      End If
      nTit=nTit+1
      If ( nTit.le.mxTit ) Title(nTit)=Trim(Line)
      Goto 15
*---  Process the "FROZen orbitals" command ---------------------------*
1020  Continue
*
      If (iPrint.GE.0) then
         Write(6,*)
         Write(6,'(6X,A)')'*** WARNING: Default frozen orbitals is '//
     &                    'overwritten by user input.'
         Write(6,'(6X,A,8I4)')
     &                    '*** Default values:',(nFro(iSym),iSym=1,nSym)
      EndIf
*
25    Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 25
      Read(Line,*,Err=994) (nFro(iSym),iSym=1,nSym)
      Goto 100
*---  Process the "DELEted orbitals" command --------------------------*
1030  Continue
35    Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 35
      Read(Line,*,Err=994) (nDel(iSym),iSym=1,nSym)
      Goto 100
*---  Process the "PRINt level" command -------------------------------*
1040  Continue
45    Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 45
      Read(Line,*,Err=994) iPrint
      Goto 100
*---  Process the "MOLOrb" command ------------------------------------*
1050  Continue
      iVecTyp=1
      Goto 100
*---  Process the "LUMOrb" command ------------------------------------*
1060  Continue
      iVecTyp=2
      Goto 100
*---  Process the "JOBIph" command ------------------------------------*
1070  Continue
      iVecTyp=3
      Goto 100
*---  Process the "ONEL only" command ---------------------------------*
1080  Continue
      iOneOnly=1
      Goto 100
*---  Process the "FILEORB" command------------------------------------*
1090  Continue
      iVecTyp=2
      Line=Get_Ln(LuSpool)
      write(6,*)' RdInp_Motra before calling fileorb.'
      write(6,*)' Line:'//line(1:60)
      write(6,*)'   Calling fileorb now...'
      call fileorb(Line,FnInpOrb)
      write(6,*)'   Back from fileorb.'
      write(6,*)' Line:'//line(1:60)
      write(6,*) FnInpOrb
      Goto 100
*---  Process the "AUTO delete" command--------------------------------*
1100  Continue
      iAutoCut=1
105   Read(LuSpool,'(A)') Line
      If ( Line(1:72).eq.Blank(1:72) .or. Line(1:1).eq.'*' ) Goto 105
      Read(Line,*,Err=994) (CutThrs(iSym),iSym=1,nSym)
      Goto 100
*---  Process the "EXTRact" command------------------------------------*
1110  Write (6,*) 'The EXTRACT option is redundant and is ignored!'
      Goto 100
*---  Process the "RFperturbation" command ----------------------------*
1120  Continue
      iRFpert=1
      Goto 100
*---  Process the "CTonly" to perform exclusively CD vectors transform-*
1130  Continue
      Line=Get_Ln(LuSpool)
      Call UpCase(Line)
      Call LeftAd(Line)
      tv2disk=Line(1:3)
      If (tv2disk.ne.'PQK' .and. tv2disk.ne.'KPQ') tv2disk='PQK' !def
      iCTonly=1
      Goto 100
*---  Process the "DIAGonal ERI evaluation" command -------------------*
1140  Continue
      iDoInt=1
      Goto 100
*---  Process the "HDF5 output file" command --------------------------*
1150  Continue
      ihdf5 = 1
      Goto 100
*---  Process the "END of input" command ------------------------------*
2000  Continue

* New rules for title lines...warning needed?
      If ( nTit.gt.mxTit) Then
       Write(6,*)' NOTE: New input specifications says TITLE keyword'
       Write(6,*)' must be followed by exactly one title line.'
       nTit=mxTit
      End If
*
*---- Check for deleted orbitals
*
      Call Get_iArray('nDel',nDel2,nSym)
      Do iSym = 1, nSym
         If (nDel2(iSym).gt.nDel(iSym)) Then
            Write (6,*)
            Write (6,*) 'Orbitals deleted at an earlier stage.'
            Write (6,*) 'iSym=',iSym
            Write (6,*) 'Input value  :',nDel(iSym)
            Write (6,*) 'Earlier value:',nDel2(iSym)
            Write (6,*)
            Write (6,*) 'Input value is now updated to earlier value!'
            Write (6,*)
            nDel(iSym)=nDel2(iSym)
         End If
      End Do
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      nOrbt=0
      nOrbtt=0
      Do iSym=1,nSym
        nOrb(iSym)=nBas(iSym)-nFro(iSym)-nDel(iSym)
        nOrbt=nOrbt+nOrb(iSym)
        nOrbtt=nOrbtt+nOrb(iSym)*(nOrb(iSym)+1)/2
      End Do
      Call Put_iArray('nFro',nFro,nSym)
      close(LuSpool)
      Return
*----------------------------------------------------------------------*
*     Error Exit                                                       *
*----------------------------------------------------------------------*
994   Write (6,*) 'RdInp: error readin input file!'
      Write (6,*) 'Command=',CmdTab(jCmd)
      Call Abend()
      End
