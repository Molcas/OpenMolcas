!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

! Module replacing the code in prgminit.c in molcas-extra

#include "compiler_features.h"
#include "macros.fh"
! from getenvc.c
#define MAXSTR 256

module prgm

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
private

type FileEntry
  character(len=MAXSTR) :: Filename, Shortname
  character(len=16) :: Attr
end type FileEntry
type(FileEntry), allocatable :: FileTable(:)
character(len=MAXSTR) :: WorkDir = '', FastDir = '', Project = 'Noname', StatusFile = ''
character(len=16) :: SlaveDir = '', SubDir = ''

public :: PrgmInitC, PrgmFree, PrgmTranslate_Mod, SetSubDir
#ifndef _GA_
public :: IsInMem
#endif

! Private extensions to mma interfaces

interface cptr2loff
  module procedure :: fe_cptr2loff
end interface
interface mma_allocate
  module procedure :: fe_mma_allo_1D, fe_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: fe_mma_free_1D
end interface

contains

subroutine PrgmInitC(ModName,l)
  character(len=*), intent(in) :: ModName
  integer(kind=iwp), intent(in) :: l
  integer(kind=iwp) :: n
  unused_var(l)
  call PrgmCache()
  call ReadPrgmFile(ModName)
  call ReadPrgmFile('global')
  ! Save location of the "status" file, since it will be needed after deallocation
  call PrgmTranslate_Mod('status',5,StatusFile,n,0)
  return
end subroutine PrgmInitC

subroutine PrgmFree()
  if (allocated(FileTable)) call mma_deallocate(FileTable)
  return
end subroutine PrgmFree

subroutine PrgmTranslate_Mod(InStr,l1,OutStr,l2,Par)
  character(len=*), intent(in) :: InStr
  integer(kind=iwp), intent(in) :: l1, Par
  character(len=*), intent(out) :: OutStr
  integer(kind=iwp), intent(out) :: l2
  integer(kind=iwp) :: Num, Loc
  logical(kind=iwp) :: Found, Lustre
  character(len=len(InStr)) :: Input
  character(len=MAXSTR) :: WD, Ext
  character(len=16) :: Attr
  unused_var(l1)
  Lustre = .false.
  Input = InStr
  Loc = index(Input,char(0))
  if (Loc > 0) Input(Loc:) = ''
# ifdef _DEBUGPRINT_
  write(u6,*) 'Translating ',trim(Input)
# endif
  inquire(file=Input,exist=Found)
  if (Found) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'File ',trim(Input),' exists'
#   endif
    OutStr = Input
  else
    WD = WorkDir
    if (trim(WD) == '') WD = '.'
    if (allocated(FileTable)) then
      Num = FindFile(Input,FileTable)
    else
      Num = -1
    end if
    if (Num > 0) then
#     ifdef _DEBUGPRINT_
      write(u6,*) 'File ',trim(Input),' found in table: ',trim(FileTable(Num)%Filename)
#     endif
      ! the value of $WorkDir will depend on some attributes,
      ! pass it to the ExpandVars function
      Attr = FileTable(Num)%Attr
      if (index(Attr,'f') > 0) WD = FastDir
#     ifdef MOLCAS_LUSTRE
      Lustre = (index(Attr,'l') > 0)
#     endif
      if ((Par == 1) .and. (.not. Lustre)) WD = trim(WD)//SlaveDir
      OutStr = FileTable(Num)%Filename
      OutStr = ExpandVars(OutStr,trim(WD)//SubDir)
      ! add "extension" to multi files
      if (index(Attr,'*') > 0) then
        Ext = Input(len_trim(FileTable(Num)%Shortname)+1:)
        OutStr = trim(OutStr)//Ext
      else if (index(Attr,'.') > 0) then
        Ext = Input(len_trim(FileTable(Num)%Shortname)+1:)
        Loc = index(OutStr,'.',back=.true.)
        OutStr = ReplaceSubstr(OutStr,Loc,Loc,trim(Ext)//'.')
      end if
    else if ((Num < 0) .and. (Input == 'status')) then
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Using cached "status" entry'
#     endif
      OutStr = trim(StatusFile)
    else
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Assuming $WorkDir/'//trim(Input)
#     endif
      if ((Par == 1)) WD = trim(WD)//SlaveDir
      OutStr = ExpandVars('$WorkDir/'//Input,trim(WD)//SubDir)
    end if
  end if
  l2 = len_trim(OutStr)
# ifdef _DEBUGPRINT_
  write(u6,*) 'Output: ',trim(OutStr)
# endif
  return
end subroutine PrgmTranslate_Mod

subroutine SetSubDir(Dir)
  character(len=*), intent(in) :: Dir
  SubDir = trim(Dir)
  return
end subroutine SetSubDir

#ifndef _GA_
function IsInMem(Filename)
  integer(kind=iwp) :: IsInMem
  character(len=*), intent(in) :: Filename
  unused_var(len(Filename)) ! GCC 4.8 doesn't like just "unused_var(Filename)" here
  ! since FiM is not in OpenMolcas, always return 0
  IsInMem = 0
  return
end function IsInMem
#endif

! Private procedures follow

! Subroutine to read the contents of a .prgm file and append them to FileTable
subroutine ReadPrgmFile(ModName)
  character(len=*), intent(in) :: ModName
  character(len=MAXSTR) :: Dir, Line
  character(len=2*MAXSTR) :: FileName
  logical(kind=iwp) :: Found
  integer(kind=iwp) :: Lu, Error, Num, i, j, k
  integer(kind=iwp), external :: isFreeUnit
  type(FileEntry), allocatable :: TempTable(:), NewTable(:)

# include "macros.fh"
  unused_proc(mma_allocate(FileTable,[0,0]))
  if (.not. allocated(FileTable)) call mma_allocate(FileTable,0,label='FileTable')

# ifdef _DEBUGPRINT_
  write(u6,*) 'ModName: '//trim(ModName)
# endif

  ! If other locations are enabled, pymolcas should be changed too
  ! Or this should be revamped so that pymolcas writes a pre-parsed file
  call GetEnvF('MOLCAS',Dir)
  Dir = trim(Dir)//'/data'
  FileName = trim(Dir)//'/'//trim(ModName)//'.prgm'
  inquire(file=FileName,exist=Found)
  if (Found) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'FileName: '//trim(FileName)
#   endif
    Lu = isFreeUnit(30)
    call molcas_open(Lu,trim(FileName))
    Num = 0
    do
      read(Lu,*,iostat=Error)
      if (Error /= 0) exit
      Num = Num+1
    end do
    call mma_allocate(TempTable,Num,label='TempTable')
    rewind(Lu)
    Num = 0
    do
      read(Lu,'(A)',iostat=Error) Line
      if (Error /= 0) exit
      Line = adjustl(Line)
      if (Line(1:1) == '#') cycle
      if (index(Line,'(prgm)') /= 0) cycle
      if (index(Line,'(file)') == 0) cycle
      Num = Num+1
      ! Strip quotes and tabs
      Line = Strip(Line,'"'//char(9))
      Line = adjustl(Line(index(Line,' '):))
      TempTable(Num)%Shortname = Line(:index(Line,' '))
      Line = adjustl(Line(index(Line,' '):))
      TempTable(Num)%Filename = Line(:index(Line,' '))
      Line = adjustl(Line(index(Line,' '):))
      TempTable(Num)%Attr = Line(:index(Line,' '))
    end do
    ! Make sure leftover entries are blanked
    do i=Num+1,size(TempTable)
      TempTable(i)%Shortname = ''
    end do
    ! Count new unique entries
    j = 0
    do i=1,Num
      if (FindFile(TempTable(i)%Shortname,FileTable,Exact=.true.) > 0) cycle
      if (FindFile(TempTable(i)%Shortname,TempTable(1:i-1),Exact=.true.) > 0) cycle
      j = j+1
    end do
    Num = j
    ! Merge and update FileTable
    call mma_allocate(NewTable,size(FileTable)+Num,label='FileTable')
    NewTable(1:size(FileTable)) = FileTable
    k = size(FileTable)
    do i=1,size(TempTable)
      if (TempTable(i)%Shortname == '') exit
      Num = k+1
      j = FindFile(TempTable(i)%Shortname,NewTable(1:k),Exact=.true.)
      if (j > 0) Num = j
      NewTable(Num) = TempTable(i)
      k = max(k,Num)
    end do
    call mma_deallocate(FileTable)
    call move_alloc(NewTable,FileTable)
    call mma_deallocate(TempTable)
    close(Lu)
  end if
# ifdef _DEBUGPRINT_
  call ReportPrgm()
# else
  unused_proc(ReportPrgm())
# endif
  return
end subroutine ReadPrgmFile

! Subroutine to print the contents of FileTable
subroutine ReportPrgm()
  character(len=MAXSTR) :: FmtStr
  integer(kind=iwp) :: MaxLen(2)
  integer(kind=iwp) :: i
  MaxLen = 0
  do i=1,size(FileTable)
    MaxLen(1) = max(MaxLen(1),len_trim(FileTable(i)%Shortname))
    MaxLen(2) = max(MaxLen(2),len_trim(FileTable(i)%Filename))
  end do
  write(FmtStr,*) '(A',MaxLen(1),',1X,A',MaxLen(2),',1X,"<",A,">")'
  write(u6,*) '===================================='
  write(u6,*) 'Number of entries: ',size(FileTable)
  do i=1,size(FileTable)
    write(u6,FmtStr) FileTable(i)%Shortname,FileTable(i)%Filename,trim(FileTable(i)%Attr)
  end do
  write(u6,*) '===================================='
  return
end subroutine ReportPrgm

! Function to strip all the characters in Chars from String (in any position)
function Strip(String,Chars)
# ifdef ALLOC_ASSIGN
  character(len=:), allocatable :: Strip
# else
  character(len=len(String)) :: Strip
# endif
  character(len=*), intent(in) :: String, Chars
  character(len=len(String)) :: Tmp
  integer(kind=iwp) :: i, j
  j = 0
  do i=1,len_trim(String)
    if (index(Chars,String(i:i)) /= 0) cycle
    j = j+1
    Tmp(j:j) = String(i:i)
  end do
  Strip = trim(Tmp(1:j))
  return
end function Strip

! Function to find a file, given its Shortname, in a file table
! Returns the index in the table (or 0 if not found)
function FindFile(Short,Table,Exact)
  integer(kind=iwp) :: FindFile
  character(len=*), intent(in) :: Short
  type(FileEntry), intent(in) :: Table(:)
  logical(kind=iwp), intent(in), optional :: Exact
  logical(kind=iwp) :: FindExact
  integer(kind=iwp) :: i
  FindFile = 0
  if (present(Exact)) then
    FindExact = Exact
  else
    FindExact = .false.
  end if
  do i=1,size(Table)
    ! an entry matches not only if it's equal, also if
    ! the beginning matches and it is a "multi" file
    if (FindExact) then
      if (trim(Short) == trim(Table(i)%Shortname)) then
        FindFile = i
        exit
      end if
    else
      if (index(Short,trim(Table(i)%Shortname)) == 1) then
        if ((trim(Short) == trim(Table(i)%Shortname)) .or. (index(Table(i)%Attr,'*') > 0) .or. &
            (index(Table(i)%Attr,'.') > 0)) then
          FindFile = i
          exit
        end if
      end if
    end if
  end do
  return
end function FindFile

! Function to replace environment variables in a string
! A variable starts with '$' and ends with [ $/.] or the end of the string
function ExpandVars(String,WD)
# ifdef ALLOC_ASSIGN
  character(len=:), allocatable :: ExpandVars
# define LEV Len(ExpandVars)
# else
  character(len=MAXSTR) :: ExpandVars
# define LEV Len_Trim(ExpandVars)
# endif
  character(len=*), intent(in) :: String, WD
  character(len=MAXSTR) :: Var, Val
  integer(kind=iwp) :: i, j, Ini, Fin
  ! start with a copy of the string, this copy will be modified on the fly
  ExpandVars = String
  Ini = 0
  i = 0
  ! scan the string, not using a fixed loop because the length and index
  ! will change as replacements are done
  do
    i = i+1
    if (i > LEV) exit
    ! the $ marks the start of the variable, find the end
    ! (by default the end of the string)
    if (ExpandVars(i:i) == '$') then
      Ini = i
      Fin = LEV
      do j=i+1,LEV
        if (index(' $/.',ExpandVars(j:j)) /= 0) then
          Fin = j-1
          exit
        end if
      end do
      ! get the value of the variable and replace it in the string
      Var = ExpandVars(Ini+1:Fin)
      select case (Var)
        case ('Project')
          Val = Project
        case ('WorkDir')
          Val = WD
        case Default
          call GetEnvF(Var,Val)
          if (trim(Val) == '') then
            if (Var /= 'SubProject') Val = 'UNK_VAR'
          end if
      end select
      ExpandVars = ReplaceSubstr(ExpandVars,Ini,Fin,trim(Val))
      ! update the index to continue
      i = Ini+len_trim(Val)-1
    end if
  end do
  return
end function ExpandVars

! Function to replace a substring (between the Ini and Fin positions) with Repl
function ReplaceSubstr(String,Ini,Fin,Repl)
# ifdef ALLOC_ASSIGN
  character(len=:), allocatable :: ReplaceSubstr
# else
  character(len=MAXSTR) :: ReplaceSubstr
# endif
  character(len=*), intent(in) :: String, Repl
  integer(kind=iwp), intent(in) :: Ini, Fin
  integer(kind=iwp) :: i, j
  ! make sure the indices are within limits
  i = min(max(1,Ini),len(String))
  j = min(max(1,Fin),len(String))
  j = max(i,j)
  ReplaceSubstr = trim(String(:i-1)//Repl//String(j+1:))
  return
end function ReplaceSubstr

! Save some often used variables as module variables, for faster access
subroutine PrgmCache()
  use Para_Info, only: mpp_id
  call GetEnvF('WorkDir',WorkDir)
  call GetEnvF('FastDir',FastDir)
  call GetEnvF('Project',Project)
  if (trim(Project) == '') Project = 'Noname'
  if (mpp_id() > 0) write(SlaveDir,'(A,I0)') '/tmp_',mpp_id()
  return
end subroutine PrgmCache

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define fe_cptr2loff, fe_mma_allo_1D, fe_mma_allo_1D_lim, fe_mma_free_1D
#define _TYPE_ type(FileEntry)
#  define _FUNC_NAME_ fe_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ fe_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'fe_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module prgm
