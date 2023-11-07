!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************
      SubRoutine ExpArr(Array,Ind,nArray,lArray)
!***********************************************************************
!                                                                      *
! Object: to expand arrays according to an index array.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             Augusti '91                                              *
!***********************************************************************
      use Constants, only: Zero
      Implicit None
      Integer lArray, nArray
      Real*8 Array(lArray,nArray)
      Integer Ind(nArray)

      Integer iArray, jArray
!
      Do 100 iArray = nArray, 1, -1
         jArray = Ind(iArray)
         If (jArray.le.0) Then
!           Set column iArray to zero
            Array(:,iArray)=Zero
         Else If (jArray.lt.iArray) Then
!           Copy row jArray to position iArray
            Array(:,iArray)=Array(:,jArray)
         End If
 100  Continue
      Return
      End SubRoutine ExpArr
