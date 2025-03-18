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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
      Subroutine CMSRdMat(Mat,NRow,NCol,FileName,NameLen)
      Implicit None
      INTEGER NRow,NCol,NameLen
      Real*8,DIMENSION(NRow*NCol)::Mat
      CHARACTER(Len=NameLen)::FileName
      INTEGER I,J,LU
      Integer, External :: IsFreeUnit

      LU=233
      LU=IsFreeUnit(LU)
      CALL Molcas_Open(LU,FileName)
      DO I=1,NRow
       read(LU,*)(Mat((I-1)*NCol+J),J=1,NCol)
      END DO
      CLOSE(LU)

      END Subroutine CMSRdMat
