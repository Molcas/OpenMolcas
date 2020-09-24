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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_RdInp()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Parse and process RPA input.
C
      Implicit None
#include "WrkSpc.fh"
#include "rpa_config.fh"
#include "rpa_data.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_RdInp')
      Integer mLine
      Parameter (mLine=100000)

      Logical Reduce_Prt
      External Reduce_Prt
      Character*180 Get_Ln
      External Get_Ln
      Integer  iPrintLevel, RPA_iUHF
      External iPrintLevel, RPA_iUHF

      Character*4   Keyword
      Character*180 Line

      Logical doTitle
      Logical EndOfInput
      Logical DebugPrint

      Integer ip_Integer, l_Integer
      Integer Lu
      Integer iLine
      Integer iUHF
      Integer i

      ! register entry
      Call qEnter(SecNam)

      ! set default print level
      iPrint=max(iPrintLevel(-1),0)
      If (iPrint.gt.0 .and. iPrint.lt.3) Then
         If (Reduce_Prt()) Then
           iPrint=0
         End If
      End If

      ! set debug print
#if defined (_DEBUGPRINT_)
      DebugPrint=iPrint.gt.4
#else
      DebugPrint=.false.
#endif
      If (DebugPrint) Then
         Write(6,'(A,A,I3)') SecNam,': default print level is',iPrint
      End If

      ! Set restricted (1) or unrestricted (2)
      iUHF=RPA_iUHF()

      ! set defaults
      dRPA=.true.
      SOSEX=.false.
      LumOrb=.false.
      nTitle=0
      nFreeze(1)=0
      nFreeze(2)=0

      ! open input file and find RPA section
      Lu=17
      Call SpoolInp(Lu)
      Rewind(Lu)
      Call RdNLst(Lu,'RPA')

      ! allocate dummy arrays for reading input
      ! (dimension should be the longest needed for reading)
      l_Integer=2
      Call GetMem('InpRdI','Allo','Inte',ip_Integer,l_Integer)

      ! parse input
      doTitle=.false.
      EndOfInput=.false.
      iLine=0
      Do While (.not.EndOfInput .and. iLine.lt.mLine)
         iLine=iLine+1 ! eliminate risk of an infinite loop
         Line=Get_Ln(Lu)
         Call StdFmt(Line,Keyword)
         If (DebugPrint) Then
            Write(6,'(A,A)') SecNam,' processing line:'
            Write(6,'(A)') Line
            Write(6,'(A,A)') 'Key=',Keyword
         End If
         !****************************
         If (Keyword(1:1).eq.' ') Then
         !****************************
            ! blank line: pass
            Continue
         !*********************************
         Else If (Keyword(1:1).eq.'*') Then
         !*********************************
            ! comment line: pass
            Continue
         !***********************************
         Else If (Keyword(1:3).eq.'END') Then
         !***********************************
            ! end of input
            EndOfInput=.true.
         !*******************************
         Else If (Keyword.eq.'TITL') Then
         !*******************************
            ! title
            doTitle=.true.
            nTitle=nTitle+1
            If (nTitle.gt.mTitle) Then
               Write(6,'(A,I8)') 'Maximum number of title lines is',
     *                           mTitle
               Write(6,'(A)') 'Current input line:'
               Write(6,'(A)') Line
               Call RPA_Warn(2,'Too many title lines in RPA input')
            Else
               Line=Get_Ln(Lu)
               Title(nTitle)=Line(1:80)
            End If
         !*******************************
         Else If (Keyword.eq.'PRIN') Then
         !*******************************
            ! print level
            doTitle=.false.
            Call RPA_ReadIntegerInput('PRIN',1,Lu,
     *                                iWork(ip_Integer),l_Integer)
            iPrint=max(iWork(ip_Integer),0)
#if defined (_DEBUGPRINT_)
            DebugPrint=DebugPrint.or.iPrint.gt.4
#endif
         !*******************************
         Else If (Keyword.eq.'LUMO') Then
         !*******************************
            ! use orbitals from InpOrb file
            doTitle=.false.
            LumOrb=.true.
         !*******************************
         Else If (Keyword.eq.'RUNO') Then
         !*******************************
            ! use orbitals from Runfile
            doTitle=.false.
            LumOrb=.false.
         !*******************************
         Else If (Keyword.eq.'DRPA') Then
         !*******************************
            ! direct RPA
            doTitle=.false.
            dRPA=.true.
            SOSEX=.false.
         !*******************************
         Else If (Keyword.eq.'SOSE') Then
         !*******************************
            ! direct RPA + SOSEX
            doTitle=.false.
            dRPA=.true.
            SOSEX=.true.
         !*******************************
         Else If (Keyword.eq.'ALLE') Then
         !*******************************
            ! All electrons correlated (no frozen)
            doTitle=.false.
            Call iZero(nFro,16)
         !*******************************
         Else If (Keyword.eq.'FREE') Then
         !*******************************
            ! Freeze occupied orbitals
            doTitle=.false.
            Call iZero(nFro,16)
            Call RPA_ReadIntegerInput('FREE',iUHF,Lu,
     *                                iWork(ip_Integer),l_Integer)
            Do i=1,iUHF
               nFreeze(i)=iWork(ip_Integer-1+i)
            End Do
         !*******************************
         Else If (Keyword.eq.'DELE') Then
         !*******************************
            ! Delete virtual orbitals
            doTitle=.false.
            Call RPA_Warn(2,
     *                  'Virtual orbital deletion not implemented yet!')
         Else
            ! no keyword match
            ! => either a title line or an unknown keyword
            If (doTitle) Then
               ! process as title line
               nTitle=nTitle+1
               If (nTitle.gt.mTitle) Then
                  Write(6,'(A,I8)') 'Maximum number of title lines is',
     *                              mTitle
                  Write(6,'(A)') 'Current input line:'
                  Write(6,'(A)') Line
                  Call RPA_Warn(2,'Too many title lines in RPA input')
               Else
                  Title(nTitle)=Line(1:80)
               End If
            Else
               ! unknown keyword
               Write(6,'(A)') 'Offending input line:'
               Write(6,'(A)') Line
               Write(6,'(A,A,A)') 'Equivalent keyword input "',Keyword,
     *                            '" not recognized!'
               Call RPA_Warn(2,'RPA input keyword not recognized')
            End If
         End If
      End Do

      ! deallocation
      Call GetMem('InpRdI','Free','Inte',ip_Integer,l_Integer)

      ! close input file
      Call Close_LuSpool(Lu)

      ! register exit
      Call qExit(SecNam)

      End
************************************************************************
      Subroutine RPA_ReadIntegerInput(Key,nInp,Lu,iVal,n)
      Implicit None
      Character*4 Key
      Integer nInp
      Integer Lu
      Integer n
      Integer iVal(n)
      Character*180 Line

      Character*180 Get_Ln
      External Get_Ln

      If (n.ge.nInp) Then
         Line=Get_Ln(Lu)
         Call Get_I(1,iVal,nInp)
      Else
         ! insufficent memory for reading (fix in calling routine)
         Call RPA_Warn(3,'Integer read problem for keyword '//Key)
      End If

      End
