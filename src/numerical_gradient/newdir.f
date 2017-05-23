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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
* Interfaces to C functions needed to create and remove
* a subdirectory in WorkDir
      Module NewDir
      Use iso_c_binding
      Implicit None
#include "molcastypes.fh"
      Interface
        Function c_getcwd(path,size) bind(C,name="getcwd_for_f")
          Use iso_c_binding
          Type(c_ptr) :: c_getcwd
          Character(Kind=c_char) :: path(*)
          Integer(Kind=MOLCAS_C_INT) :: size
        End Function c_getcwd
        Function c_chdir(path) bind(C,name="chdir_for_f")
          Use iso_c_binding
          Integer(Kind=MOLCAS_C_INT) :: c_chdir
          Character(Kind=c_char) :: path(*)
        End Function c_chdir
        Function c_mkdir(path,mode) bind(C,name="mkdir_for_f")
          Use iso_c_binding
          Integer(Kind=MOLCAS_C_INT) :: c_mkdir
          Character(Kind=c_char) :: path(*)
          Integer(Kind=MOLCAS_C_INT) :: mode
        End Function c_mkdir
        Function c_rmrf(path) bind(C,name="rmrf")
          Use iso_c_binding
          Integer(Kind=MOLCAS_C_INT) :: c_rmrf
          Character(Kind=c_char) :: path(*)
        End Function c_rmrf
        Subroutine c_setsubdir(sub) bind(C,name="setsubdir")
          Use iso_c_binding
          Character(Kind=c_char) :: sub(*)
        End Subroutine c_setsubdir
      End Interface
      Character(Len=1024) :: Sub, OldWorkDir, NewWorkDir

      Contains

      Subroutine f_getcwd(path)
        Use iso_c_binding
        Implicit None
        Character(Kind=c_char,Len=*) :: path
        Type(c_ptr) :: ret
        Integer :: i
        ret = c_getcwd(path, Len(path))
        Do i=Len(path),1,-1
          If (path(i:i).eq.c_null_char) path(i:i)=' '
        End Do
      End Subroutine f_getcwd

      Subroutine f_chdir(path, err)
        Use iso_c_binding
        Implicit None
        Character(*) :: path
        Integer, Optional, Intent(Out) :: err
        Integer :: loc_err
        loc_err = c_chdir(Trim(path)//c_null_char)
        If (Present(err)) err = loc_err
      End Subroutine f_chdir

      Subroutine f_mkdir(path, err)
        Use iso_c_binding
        Implicit None
        Character(*) :: path
        Integer, Optional, Intent(Out) :: err
        Integer :: loc_err
        loc_err = c_mkdir(Trim(path)//c_null_char, Int(o'772'))
        If (Present(err)) err = loc_err
      End Subroutine f_mkdir

      Subroutine f_rmrf(path, err)
        Use iso_c_binding
        Implicit None
        Character(*) :: path
        Integer, Optional, Intent(Out) :: err
        Integer :: loc_err
        loc_err = c_rmrf(Trim(path)//c_null_char)
        If (Present(err)) err = loc_err
      End Subroutine f_rmrf

      Subroutine f_setsubdir(sub)
        Use iso_c_binding
        Implicit None
        Character(*) :: sub
        If (Trim(sub).eq.'') Then
          Call c_setsubdir(''//c_null_char)
        Else
          Call c_setsubdir('/'//Trim(sub)//c_null_char)
        End If
      End Subroutine f_setsubdir

      End Module NewDir
