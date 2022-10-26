!**********************************************************************
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
! Jie J. Bao, on May. 19, 2022, created this file.               *
!*****************************************************************
!***********************************************************************
      Subroutine ReadMat2(FileName,MatInfo,Matrix,NRow,NCol,            &
     &LenName,LenInfo,Trans)

!     This subroutine is to replace ReadMat in the long run.
      INTEGER NRow,NCol,LenName
      CHARACTER(Len=LenName)::FileName
      CHARACTER(Len=LenInfo)::MatInfo
      CHARACTER(Len=1)::Trans
      Real*8,DIMENSION(NRow*NCol)::Matrix

      INTEGER LU,IsFreeUnit,IRow,ICol,iOff
      External IsFreeUnit

      IF(LenName.gt.0) THEN
       LU=100
       LU=IsFreeUnit(LU)
       CALL Molcas_Open(LU,FileName)
      ELSE
       LU=6
      END IF
      IF(Trans.eq.'T') THEN
       DO IRow=1,NRow
        iOff=(IRow-1)*nCol
        read(LU,*) (Matrix(iOff+ICol),ICol=1,NCol)
       END DO
      ELSE
       DO ICol=1,NCol
        read(LU,*) (Matrix((iRow-1)*nCol+iCol),IRow=1,NRow)
       END DO
      END IF
      Read(LU,*)MatInfo
      IF(LenName.gt.0) THEN
       Close(LU)
      END IF
      RETURN
      End Subroutine
