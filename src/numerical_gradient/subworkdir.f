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
      Subroutine SubWorkDir
      use subdirs, only : f_setsubdir, Sub, OldWorkDir, NewWorkDir
      use filesystem, only : getcwd_, chdir_, mkdir_
      Implicit None
      Integer :: i,Length,iErr
      Integer, Parameter :: nFiles=21
      Character(Len=1024) :: Names(nFiles),
     &                       OldFile(nFiles),NewFile(nFiles)
      Logical :: Found

*     Define name of subdirectory and files that must be copied over
      Sub='NG'
      Names( 1)='RUNFILE'
      Names( 2)='SEWARINP'
      Names( 3)='SCFINP'
      Names( 4)='RASSCINP'
      Names( 5)='CASPTINP'
      Names( 6)='MBPT2INP'
      Names( 7)='RASSIINP'
      Names( 8)='MOTRAINP'
      Names( 9)='CCSDTINP'
      Names(10)='CHCCINP'
      Names(11)='CHT3INP'
      Names(12)='ESPFINP'
      Names(13)='JOBIPH'
      Names(14)='ESPF.SAV'
      Names(15)='TINKER.XYZ'
      Names(16)='TINKER.KEY'
      Names(17)='MCPDFINP'
      Names(18)='CHEMNATFIE'
      Names(19)='CHEMCANFIE'
      Names(20)='CHEMNATMPS0'
      Names(21)='CHEMCANMPS0'

*     Get real filenames to copy
      Do i=1,nFiles
        Call prgmtranslate(Names(i),OldFile(i),Length)
      End Do

*     Create the new directory and switch to it
      Call getcwd_(OldWorkDir)
      NewWorkDir=Trim(OldWorkDir)//'/'//Trim(Sub)
      Call mkdir_(NewWorkDir)
      Call chdir_(NewWorkDir)
      Call f_setsubdir(Sub)

*     Get real target filenames
      Do i=1,nFiles
        Call prgmtranslate(Names(i),NewFile(i),Length)
*       ESPF.SAV is copied to ESPF.DATA
        If (Names(i).eq.'ESPF.SAV') Then
          Call prgmtranslate('ESPF.DATA',NewFile(i),Length)
        End If
      End Do

*     Copy the files from the old directory to the new
      Do i=1,nFiles
        Call f_inquire(OldFile(i),Found)
        If (Found) Call fcopy(OldFile(i),NewFile(i),iErr)
      End Do

*     The INPORB file is special...
      Call f_inquire('../INPORB',Found)
      If (Found) Call fcopy('../INPORB','./INPORB',iErr)

      End Subroutine SubWorkDir
