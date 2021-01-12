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
      Character*(*) FileName

      Call Ext_PID(FileName)
      Return
      End
      Subroutine Ext_PID(FileName)
      Implicit Real*8 (a-h,o-z)
#ifdef _MOLCAS_MPP_
      External StrnLn
      Integer StrnLn
#endif
#include "unixinfo.fh"
      Character*(*) FileName

#ifdef _MOLCAS_MPP_
      Length=Len(FileName)
      NameLength=StrnLn(FileName)
      Filename(NameLength+1:NameLength+1)='_'
      Call WrNumber(Filename(NameLength+2:Length),PID)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_character(FileName)
#endif
      Return
      End
      Subroutine Ext_Rank(FileName)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: MyRank
#endif
      Implicit Real*8 (a-h,o-z)
#ifdef _MOLCAS_MPP_
      External StrnLn
      Integer StrnLn
#endif
      Character*(*) FileName

#ifdef _MOLCAS_MPP_
      Length=Len(FileName)
      NameLength=StrnLn(FileName)
      Filename(NameLength+1:NameLength+1)='_'
      Call WrNumber(FileName(NameLength+2:Length),MyRank)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_character(FileName)
#endif
      Return
      End
      Subroutine WrNumber(Name,Number)
      Implicit Real*8 (a-h,o-z)
      Character*(*) Name
      Character*10 Format

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
      End
