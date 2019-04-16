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
      Subroutine ParentWorkDir
      use subdirs, only : f_setsubdir, Sub, OldWorkDir, NewWorkDir
      use filesystem, only : remove_, chdir_
      Implicit None
      Integer :: i

      Sub=''
      Call f_setsubdir(Sub)
      Call chdir_(OldWorkDir)

*     Remove the subdirectory,
*     try to make sure it is indeed a subdirectory
      If (Index(NewWorkDir,Trim(OldWorkDir)) .eq. 1) Then
        i = Len_Trim(OldWorkDir)
        If ((Len_Trim(NewWorkDir) .ge. i+2) .and.
     &      (NewWorkDir(i+1:i+1) .eq. '/') .and.
     &      (NewWorkDir(i+2:i+2) .ne. '/')) Then
          Call remove_(NewWorkDir)
        End If
      End If

      End Subroutine ParentWorkDir
