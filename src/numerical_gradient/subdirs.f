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
      Module subdirs
      Use iso_c_binding
      implicit none
      private
      public :: f_setsubdir, Sub, OldWorkDir, NewWorkDir
#include "molcastypes.fh"
      Interface
#ifdef _HAVE_EXTRA_
        Subroutine c_setsubdir(sub) bind(C, name="setsubdir")
          Use iso_c_binding
          Character(Kind=c_char) :: sub(*)
        End Subroutine c_setsubdir
#endif
      End Interface
      Character(Len=1024), save :: Sub, OldWorkDir, NewWorkDir

      Contains

#ifdef _HAVE_EXTRA_
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
#else
      Subroutine f_setsubdir(sub)
        Use Prgm
        Implicit None
        Character(*) :: sub
        If (Trim(sub).eq.'') Then
          Call SetSubDir('')
        Else
          Call SetSubDir('/'//Trim(sub))
        End If
      End Subroutine f_setsubdir
#endif

      End Module subdirs
