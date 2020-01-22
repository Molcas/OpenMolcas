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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************
Subroutine SortDiag(HH,EigVec,nVec,nDim)

Implicit None
Integer, Intent(In) :: nVec, nDim
Real*8, Intent(InOut) :: HH(*), EigVec(nDim,nVec)
Integer :: i, iMax, ii, jj
Integer, External :: idAMax_

Do i=1,nVec-1
  iMax = idAMax_(nVec-i+1,EigVec(i,i),nDim)
  If (iMax > 1) Then
    iMax = iMax+i-1
    ii = (i*(i+1))/2
    jj = (iMax*(iMax+1))/2
    Call dSwap_(1,HH(ii),1,HH(jj),1)
    Call dSwap_(nDim,EigVec(1,i),1,EigVec(1,iMax),1)
  End If
End Do

End Subroutine SortDiag
