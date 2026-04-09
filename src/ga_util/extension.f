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
      Subroutine Extension(FileName)
      implicit none
      Character(LEN=*), intent(inout):: FileName

      Call Ext_PID(FileName)
      End Subroutine Extension

      Subroutine Ext_PID(FileName)
#ifdef _MOLCAS_MPP_
      use definitions, only: iwp
      use UnixInfo, only: PID
#endif
      Implicit None
#ifdef _MOLCAS_MPP_
      Integer(kind=iwp) Length,NameLength
      Integer(kind=iwp), External ::StrnLn
#endif
      Character(LEN=*), intent(inout):: FileName

#ifdef _MOLCAS_MPP_
      Length=Len(FileName)
      NameLength=StrnLn(FileName)
      Filename(NameLength+1:NameLength+1)='_'
      Call WrNumber(Filename(NameLength+2:Length),PID)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_character(FileName)
#endif
      End Subroutine Ext_PID

      Subroutine Ext_Rank(FileName)
#ifdef _MOLCAS_MPP_
      use definitions, only: iwp
      Use Para_Info, Only: MyRank
#endif
      Implicit None
#ifdef _MOLCAS_MPP_
      Integer(kind=iwp) Length,NameLength
      Integer(kind=iwp), External ::StrnLn
#endif
      Character(LEN=*), intent(inout):: FileName

#ifdef _MOLCAS_MPP_
      Length=Len(FileName)
      NameLength=StrnLn(FileName)
      Filename(NameLength+1:NameLength+1)='_'
      Call WrNumber(FileName(NameLength+2:Length),MyRank)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_character(FileName)
#endif
      End Subroutine Ext_Rank

      Subroutine WrNumber(Name,Number)
      use definitions, only: iwp, wp
      Implicit None
      integer(kind=iwp), intent(in):: Number
      Character(LEN=*), intent(inout):: Name

      Character(LEN=10) Format
      integer(kind=iwp) Limit,iTens
      real(kind=wp) ANumber

      Format=' '
      If(Number.Ge.0)Then
        Limit=0
        Do iTens=0,99
        Limit=Limit+9*(10**iTens)
        If(Number.Le.Limit)Then
          Write(Format,'(A,I4.4,A)')'(I',iTens+1,')'
          Write(Name,Format)Number
          Return
        EndIf
        EndDo
      Else
        ANumber=-Dble(Number)
        Limit=0
        Do iTens=0,99
        Limit=Limit+9*(10**iTens)
        If(ANumber.Le.Dble(Limit))Then
          Write(Format,'(A,I4.4,A)')'(A,I',iTens+1,')'
          Write(Name,Format)'-',ANumber
          Return
        EndIf
        EndDo
      EndIf
      Call SysAbendMsg('wrnumber', 'Number too large',' ' )
      End Subroutine WrNumber
