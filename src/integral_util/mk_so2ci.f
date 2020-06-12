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
      Subroutine Mk_SO2cI(SO2cI,iSO2Shell,nSOs)
      Implicit Real*8 (a-h,o-z)
      Integer SO2cI(2,nSOs), iSO2Shell(nSOs)
*                                                                      *
************************************************************************
*                                                                      *
*-----Set up table SO to contigues index over the shell
*
*     Write (*,*) 'Enter SO2cI'
      Call ICopy(2*nSOs,[0],0,SO2cI,1)
      Call Nr_Shells(nShell)
      Do iShell = 1, nShell
*
*------- Generate contigues index for this shell
*
         Index=0
         Do iSO = 1, nSOs
            If (iSO2Shell(iSO).eq.iShell) Then
               Index = Index + 1
               SO2cI(1,iSO) = Index
            End If
         End Do
*
*------- Store dimension for this shell
*
         Do iSO = 1, nSOs
            If (iSO2Shell(iSO).eq.iShell) Then
               SO2cI(2,iSO) = Index
            End If
         End Do
*
      End Do ! iShell
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
