!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine GS_Order(T,nInter,nVec)
      Implicit Real*8 (a-h,o-z)
      Real*8 T(nInter,nVec)
!
#ifdef _DEBUGPRINT_
      Call RecPrt('GS_Order: T(orig)','(12F6.2)',T,nInter,nVec)
#endif
      Do j = 1, nVec-1
         DiagMax=DDot_(nInter,T(1,j),1,T(1,j),1)
         iDiag=j
         Do i = j+1, nVec
            Diag=DDot_(nInter,T(1,i),1,T(1,i),1)
            If (T(i,i).gt.DiagMax) Then
               DiagMax=Diag
               iDiag=i
            End If
         End Do
         If (iDiag.ne.j) call dswap_(nInter,T(1,iDiag),1,T(1,j),1)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('GS_Order: T(ordered)','(12F6.2)',T,nInter,nVec)
#endif
      Return
      End
