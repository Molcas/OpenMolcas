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
      Integer Function Index_NoSym(iCntr,iCmp,iCnt,iAng,iR,Index,iBas,
     &                             nBas)
      Integer Index(5,nBas)
*
 98   Continue
      Index_NoSym=0
      Do i = 1, iBas
         If (Index(1,i).eq.iCntr .and.
     &       Index(2,i).eq.iCmp  .and.
     &       Index(3,i).eq.iCnt  .and.
     &       Index(4,i).eq.iAng  .and.
     &       Index(5,i).eq.iR   ) Then
            Index_NoSym=i
            Go To 99
         End If
      End Do
*
      iBas=iBas+1
      Index(1,iBas)=iCntr
      Index(2,iBas)=iCmp
      Index(3,iBas)=iCnt
      Index(4,iBas)=iAng
      Index(5,iBas)=iR
      Go To 98
*
 99   Return
      End
