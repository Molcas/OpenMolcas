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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
      Subroutine PrintMat(FileName,MatInfo,Matrix,NRow,NCol,            &
     &                    LenName,LenInfo,Trans)
      Implicit None

      INTEGER NRow,NCol,LenName,LenInfo
      CHARACTER(Len=LenName)::FileName
      CHARACTER(Len=LenInfo)::MatInfo
      CHARACTER(Len=1)::Trans
      CHARACTER(Len=80)::PrtFmt
      Real*8,DIMENSION(NRow,NCol)::Matrix

      INTEGER LU,IsFreeUnit,IRow,ICol
      External IsFreeUnit


      IF(LenName.gt.0) THEN
      LU=100
      LU=IsFreeUnit(LU)
      CALL Molcas_Open(LU,FileName)
      ELSE
      LU=6
      END IF

      IF(Trans.eq.'N') THEN
       WRITE(PrtFmt,'(A,I5,A)')                                         &
     & '(',NCol,'(ES24.14E4,1X))'
       DO IRow=1,NRow
        write(LU,PrtFmt)                                                &
     &  (Matrix(IRow,ICol),ICol=1,NCol)
       END DO
      ELSE
       WRITE(PrtFmt,'(A,I5,A)')                                         &
     & '(',NRow,'(ES24.14E4,1X))'
       DO ICol=1,NCol
        write(LU,PrtFmt)                                                &
     &  (Matrix(IRow,ICol),IRow=1,NRow)
       END DO
      END IF
      WRITE(LU,*)MatInfo
      IF(LenName.gt.0) THEN
       Close(LU)
      END IF
      End Subroutine PrintMat
