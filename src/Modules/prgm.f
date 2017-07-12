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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************

* Module replacing the code in prgminit.c in molcas-extra

! from getenvc.c
#define MAXSTR 256

!define _DEBUG_

      Module prgm
      Implicit None
      Private
      Type FileEntry
        Character (Len=MAXSTR) :: Filename
        Character (Len=MAXSTR) :: Shortname
        Character (Len=16) :: Attributes
      End Type FileEntry
      Type(FileEntry), Dimension(:), Allocatable :: FileTable
      Character (Len=MAXSTR) :: WorkDir='', FastDir='', Project='Noname'
      Character (Len=16) :: SlaveDir='', SubDir=''

      Public :: PrgmInitC, PrgmFree, PrgmTranslateC, SetSubDir
#ifndef _GA_
      Public :: IsInMem
#endif

      Contains

      Subroutine PrgmInitC(ModName, l)
      Character (Len=*), Intent(In) :: ModName
      Integer, Intent(In) :: l
      Call PrgmCache()
      Call ReadPrgmFile(ModName)
      Call ReadPrgmFile('global')
      Return
      ! Avoid unused argument warnings
      If (.False.) Call Unused_integer(l)
      End Subroutine PrgmInitC

      Subroutine PrgmFree()
      If (Allocated(FileTable)) Deallocate(FileTable)
      End Subroutine PrgmFree

      Subroutine PrgmTranslateC(InStr,l1,OutStr,l2,Par)
      Character (Len=*), Intent(In) :: InStr
      Character (Len=*), Intent(Out) :: OutStr
      Integer, Intent(In) :: Par, l1
      Integer, Intent(Out) :: l2
      Character (Len=Len(InStr)) :: Input
      Character (Len=MAXSTR) :: WD, Ext
      Character (Len=16) :: Attr
      Integer :: Num, Loc
      Logical :: Found, Lustre=.False.
      Input = Strip(InStr, Char(0))
#ifdef _DEBUG_
      Write(6,*) 'Translating ', Trim(Input)
#endif
      Inquire(File=Input,Exist=Found)
      If (Found) Then
#ifdef _DEBUG_
        Write(6,*) 'File ', Trim(Input), ' exists'
#endif
        OutStr = Input
      Else
        WD = WorkDir
        Num = FindFile(Input, FileTable)
        If (Num .gt. 0) Then
#ifdef _DEBUG_
          Write(6,*) 'File ', Trim(Input), ' found in table: ',
     &               Trim(FileTable(Num)%Filename)
#endif
          ! the value of $WorkDir will depend on some attributes,
          ! pass it to the ExpandVars function
          Attr = FileTable(Num)%Attributes
          If (Index(Attr, 'f') .gt. 0) WD = FastDir
#ifdef MOLCAS_LUSTRE
          Lustre = (Index(Attr, 'l') .gt. 0)
#endif
          If ((Par .eq. 1) .and. (.not. Lustre)) WD = Trim(WD)//SlaveDir
          OutStr = FileTable(Num)%Filename
          OutStr = ExpandVars(OutStr, Trim(WD)//SubDir)
          ! add "extension" to multi files
          If (Index(Attr, '*') .gt. 0) Then
            Ext = Input(Len_Trim(FileTable(Num)%Shortname)+1:)
            OutStr = Trim(OutStr)//Ext
          Else If (Index(Attr, '.') .gt. 0) Then
            Ext = Input(Len_Trim(FileTable(Num)%Shortname)+1:)
            Loc = Index(OutStr, '.', Back=.True.)
            OutStr = ReplaceSubstr(OutStr, Loc, Loc, Trim(Ext)//'.')
          End If
        Else
#ifdef _DEBUG_
          Write(6,*) 'Assuming $WorkDir/'//Trim(Input)
#endif
          OutStr = ExpandVars('$WorkDir/'//Input, Trim(WD)//SubDir)
        End If
      End If
      l2 = Len_Trim(OutStr)
#ifdef _DEBUG_
      Write(6,*) 'Output: ', Trim(OutStr)
#endif
      Return
      ! Avoid unused argument warnings
      If (.False.) Call Unused_integer(l1)
      End Subroutine PrgmTranslateC

      Subroutine SetSubDir(Dir)
      Character (Len=*), Intent(In) :: Dir
      SubDir = Trim(Dir)
      EndSubroutine SetSubDir

#ifndef _GA_
      Function IsInMem(Filename)
      Character (Len=*), Intent(In) :: Filename
      Integer :: IsInMem
      ! since FiM is not in OpenMolcas, always return 0
      IsInMem = 0
      End Function IsInMem
#endif

! Private procedures follow

! Subroutine to read the contents of a .prgm file and append them to FileTable
      Subroutine ReadPrgmFile(ModName)
      Character (Len=*), Intent(In) :: ModName
      Character (Len=MAXSTR) :: Dir, Line
      Character (Len=2*MAXSTR) :: FileName
      Logical :: Found
      Integer :: Lu, Error, Num, i, j, k
      Integer, External :: isFreeUnit
      Type(FileEntry), Dimension(:), Allocatable :: TempTable, NewTable

      If (.Not. Allocated(FileTable)) Allocate(FileTable(0))

#ifdef _DEBUG_
      Write(6,*) 'ModName: '//Trim(ModName)
#endif

      ! If other locations are enabled, pymolcas should be changed too
      ! Or this should be revamped so that pymolcas writes a pre-parsed file
      Call GetEnvF('MOLCAS', Dir)
      Dir = Trim(Dir)//'/data'
      FileName = Trim(Dir)//'/'//Trim(ModName)//'.prgm'
      Inquire(File=FileName,Exist=Found)
      If (Found) Then
#ifdef _DEBUG_
        Write(6,*) 'FileName: '//Trim(FileName)
#endif
        Lu = isFreeUnit(30)
        Call molcas_open(Lu, Trim(FileName))
        Num = 0
        Do
          Read(Lu,*,IOStat=Error)
          If (Error .ne. 0) Exit
          Num = Num+1
        End Do
        Allocate(TempTable(Num))
        Rewind(Lu)
        Num = 0
        Do
          Read(Lu,'(A)',IOStat=Error) Line
          If (Error .ne. 0) Exit
          Line = AdjustL(Line)
          If (Line(1:1) .eq. '#') Cycle
          If (Index(Line, '(prgm)') .ne. 0) Cycle
          If (Index(Line, '(file)') .eq. 0) Cycle
          Num = Num+1
          !Strip quotes and tabs
          Line = Strip(Line, '"'//Char(9))
          Line = AdjustL(Line(Index(Line,' '):))
          TempTable(Num)%Shortname = Line(:Index(Line,' '))
          Line = AdjustL(Line(Index(Line,' '):))
          TempTable(Num)%Filename = Line(:Index(Line,' '))
          Line = AdjustL(Line(Index(Line,' '):))
          TempTable(Num)%Attributes = Line(:Index(Line,' '))
        End Do
        ! Make sure leftover entries are blanked
        Do i=Num+1,Size(TempTable)
          TempTable(i)%Shortname = ''
        End Do
        ! Count new unique entries
        j = 0
        Do i=1,Num
          If (FindFile(TempTable(i)%Shortname, FileTable) .gt. 0) Cycle
          If (FindFile(TempTable(i)%Shortname, TempTable(1:i-1)) .gt. 0)
     &       Cycle
          j = j+1
        End Do
        Num = j
        ! Merge and update FileTable
        Allocate(NewTable(Size(FileTable)+Num))
        NewTable(1:Size(FileTable)) = FileTable
        k = Size(FileTable)
        Do i=1,Size(TempTable)
          If (TempTable(i)%Shortname .eq. '') Exit
          Num = k+1
          j = FindFile(TempTable(i)%Shortname, NewTable(1:k))
          If (j .gt. 0) Num=j
          NewTable(Num) = TempTable(i)
          k = Max(k, Num)
        End Do
        Deallocate(FileTable)
        Call Move_Alloc(NewTable, FileTable)
        Deallocate(TempTable)
        Close(Lu)
      End If
#ifdef _DEBUG_
      Call ReportPrgm()
#endif
      End Subroutine ReadPrgmFile

! Subroutine to print the contents of FileTable
      Subroutine ReportPrgm()
      Character (Len=MAXSTR) :: FmtStr
      Integer, Dimension(2) :: MaxLen
      Integer :: i
      MaxLen = 0
      Do i=1,Size(FileTable)
        MaxLen(1) = Max(MaxLen(1), Len_Trim(FileTable(i)%Shortname))
        MaxLen(2) = Max(MaxLen(2), Len_Trim(FileTable(i)%Filename))
      End Do
      Write(FmtStr,*) '(A',MaxLen(1),',1X,A',MaxLen(2),',1X,"<",A,">")'
      Write(6,*) '===================================='
      Write(6,*) 'Number of entries: ', Size(FileTable)
      Do i=1,Size(FileTable)
      Write(6,FmtStr) FileTable(i)%Shortname, FileTable(i)%Filename,
     &                Trim(FileTable(i)%Attributes)
      End Do
      Write(6,*) '===================================='
      End Subroutine ReportPrgm

! Function to strip all the characters in Chars from String (in any position)
      Function Strip(String, Chars)
      Character (Len=*), Intent(In) :: String
      Character (Len=*), Intent(In) :: Chars
      Character (Len=:), Allocatable :: Strip
      Character (Len=Len(String)) :: Tmp
      Integer :: i, j
      j=0
      Do i=1,Len_Trim(String)
        If (Index(Chars, String(i:i)) .ne. 0) Cycle
        j = j+1
        Tmp(j:j) = String(i:i)
      End Do
      Strip = Trim(Tmp(1:j))
      End Function Strip

! Function to find a file, given its Shortname, in a file table
! Returns the index in the table (or 0 if not found)
      Function FindFile(Short, Table)
      Character (Len=*), Intent(In) :: Short
      Type(FileEntry), Dimension(:), Intent(In) :: Table
      Integer :: FindFile, i
      FindFile=0
      Do i=1,Size(Table)
        ! an entry matches not only if it's equal, also if
        ! the beginning matches and it is a "multi" file
        If (Index(Short, Trim(Table(i)%Shortname)) .eq. 1) Then
          If ((Trim(Short) .eq. Trim(Table(i)%Shortname)) .or.
     &        (Index(Table(i)%Attributes, '*') .gt. 0) .or.
     &        (Index(Table(i)%Attributes, '.') .gt. 0)) Then
            FindFile=i
            Exit
          End If
        End If
      End Do
      End Function

! Function to replace environment variables in a string
! A variable starts with '$' and ends with [ $/.] or the end of the string
      Function ExpandVars(String, WD)
      Character (Len=*), Intent(In) :: String, WD
      Character (Len=:), Allocatable :: ExpandVars
      Character (Len=MAXSTR) :: Var, Val
      Integer :: i, j, Ini, Fin
      ! start with a copy of the string, this copy will be modified on the fly
      ExpandVars = String
      Ini = 0
      i = 0
      ! scan the string, not using a fixed loop because the length and index
      ! will change as replacements are done
      Do
        i = i+1
        If (i .gt. Len(ExpandVars)) Exit
        ! the $ marks the start of the variable, find the end
        ! (by default the end of the string)
        If (ExpandVars(i:i) .eq. '$') Then
          Ini = i
          Fin = Len(ExpandVars)
          Do j = i+1, Len(ExpandVars)
            If (Index(' $/.', ExpandVars(j:j)) .ne. 0) Then
              Fin = j-1
              Exit
            End If
          End Do
          ! get the value of the variable and replace it in the string
          Var = ExpandVars(Ini+1:Fin)
          Select Case (Var)
            Case ('Project')
              Val = Project
            Case ('WorkDir')
              Val = WD
            Case Default
              Call GetEnvF(Var, Val)
              If (Trim(Val) .eq. '') Then
                If (Var .ne. 'SubProject') Val = 'UNK_VAR'
              End If
          End Select
          ExpandVars = ReplaceSubstr(ExpandVars, Ini, Fin, Trim(Val))
          ! update the index to continue
          i = Ini + Len_Trim(Val) - 1
        End If
      End Do
      End Function ExpandVars

! Function to replace a substring (between the Ini and Fin positions) with Repl
      Function ReplaceSubstr(String,Ini,Fin,Repl)
      Character (Len=*), Intent(In) :: String, Repl
      Integer, Intent(In) :: Ini,Fin
      Character (Len=:), Allocatable :: ReplaceSubstr
      Integer :: i,j
      ! make sure the indices are within limits
      i = Min(Max(1,Ini),Len(String))
      j = Min(Max(1,Fin),Len(String))
      j = Max(i,j)
      ReplaceSubstr = Trim(String(:i-1) // Repl // String(j+1:))
      End Function ReplaceSubstr

! Save some often used variables as module variables, for faster access
      Subroutine PrgmCache
#include "para_info.fh"
      Call GetEnvF('WorkDir', WorkDir)
      Call GetEnvF('FastDir', FastDir)
      Call GetEnvF('Project', Project)
      If (Trim(Project) .eq. '') Project = 'Noname'
      If (MyRank .gt. 0) Write(SlaveDir,'(A,I0)') '/tmp_', MyRank
      End Subroutine

      End Module prgm
