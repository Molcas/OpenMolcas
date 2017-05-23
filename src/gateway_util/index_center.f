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
      Integer Function Index_Center(iCnt,iR,Index,iAtoms,nAtoms)
      Integer Index(2,nAtoms)
*
 98   Continue
      Index_Center=0
      Do i = 1, iAtoms
         If (Index(1,i).eq.iCnt .and.
     &       Index(2,i).eq.iR   ) Then
            Index_Center=i
            Go To 99
         End If
      End Do
*
      iAtoms=iAtoms+1
      Index(1,iAtoms)=iCnt
      Index(2,iAtoms)=iR
      Go To 98
*
 99   Return
      End
