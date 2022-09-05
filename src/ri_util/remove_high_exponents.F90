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
! Copyright (C) 2007,2008, Roland Lindh                                *
!***********************************************************************
      Subroutine Remove_High_Exponents(iD,nD,List2,mData,nTheta_All)
      Use Basis_Info, only: Shells
      Implicit Real*8 (a-h,o-z)
!***********************************************************************
!                                                                      *
!     Experimental code to be used with care.                          *
!                                                                      *
!***********************************************************************
      Integer iD(nD), List2(mData,nTheta_All)
      Logical Skip
!
      Call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)
      mD = nD
      i = 1
 100  Continue
         iTheta_All=iD(i)
         Skip=.False.
         kAng  = List2(1,iTheta_All)
         lAng  = List2(2,iTheta_All)
         k     = List2(5,iTheta_All)
         l     = List2(6,iTheta_All)
         kShll = List2(7,iTheta_All)
         lShll = List2(8,iTheta_All)
         If (kAng.eq.lAng) Then
            l     = List2(6,iTheta_All)
            Skip = (k.eq.1.and.l.eq.1).and.Shells(kShll)%nExp.ne.1
         Else
            Skip=l.eq.1.and.Shells(lShll)%nExp.ne.1
         End If
         If (Skip) Then
            If (mD.eq.i) Then
               mD = mD -1
               Go To 200
            End If
            Do j = i+1, mD
               iD(j-1) = iD(j)
            End Do
            mD = mD -1
            Go To 100
         End If
         i = i + 1
         If (i.le.mD) Go To 100
 200  Continue
      nD = mD
      Call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)
!
      Return
      End
