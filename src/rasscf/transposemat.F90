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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

      Subroutine TransposeMat(Matout,Matin,nElem,nRow_in,nCol_in)
      Implicit None
      INTEGER nElem,nRow_in,nCol_in,iRow,iCol,iOff1,iOff2
      Real*8 Matin(nElem),Matout(nElem)

      IF(nRow_in*nCol_in.ne.nElem) THEN
       write(6,*) 'Error in TransposeMat()'
       write(6,*) 'nRow_in*nCol_in != nElem'
      END IF

      DO iCol=1,nCol_in
       iOff1=(iCol-1)*nRow_in
       Do iRow=1,nRow_in
        iOff2=(iRow-1)*nCol_in
        Matout(iOff2+iCol)=Matin(iOff1+iRow)
       End Do
      END DO

      RETURN
      End Subroutine
