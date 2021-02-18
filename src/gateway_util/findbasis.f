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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
          Subroutine Find_Basis_Set(DirName,ExtBasDir,type)
************************************************************************
* Object: make a guess for a file                                      *
* Author: Valera Veryazov, Theoretical Chemistry, Chemical Centre      *
*         University of Lund, Lund, Sweden                             *
************************************************************************

          Character*(*) DirName,ExtBasDir,type
          Character*512 tmp
          character Molcas*256, CurrDir*256
          integer iAbsName
          logical Exist
          if(ExtBasDir.ne.' ') then
             i=index(ExtBasDir,' ')
             iAbsName=0
             if(ExtBasDir(1:1).ne.'/') then
              iAbsName=1
              CurrDir=' '
              call getenvf('CurrDir',CurrDir)
             endif
             if(iAbsName.eq.0) then
             tmp=ExtBasDir(1:i-1)//'/'//type
             else
             tmp=CurrDir(1:index(CurrDir,' ')-1)//
     &             '/'//ExtBasDir(1:i-1)//'/'//type
             endif
c             print *,tmp
          Call f_Inquire(tmp,Exist)
          if (Exist) then
             if(iAbsName.eq.0) then
             tmp=ExtBasDir(1:i-1)
             else
             tmp=CurrDir(1:index(CurrDir,' ')-1)//
     &             '/'//ExtBasDir(1:i-1)
             endif
             i=index(tmp,' ')
             DirName=tmp(1:i-1)
c          Print *,'>',DirName(1:i-1),'<'
             return
          endif
          endif
          if(DirName(1:13).ne.'basis_library') return

          Molcas=' '
          call getenvf('MOLCAS_BASIS',Molcas)
          If (Molcas.eq.' ') Then
             call getenvf('MOLCAS',Molcas)
             i=index(Molcas,' ')
             DirName=Molcas(1:i-1)//'/basis_library'
          Else
             i=index(Molcas,' ')
             DirName=Molcas(1:i-1)
          End If
          i=index(DirName,' ')
          if(i.eq.0) then
             Call WarningMessage(2,'Too long path to Molcas')
             Call Abend()
          endif
          Return
          End
