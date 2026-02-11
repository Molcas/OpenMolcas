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
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************
      Subroutine InitRotMat(RotMat,lRoots,CMSSFile,LenCMSS)
      INTEGER LenCMSS,lRoots
      CHARACTER(Len=LenCMSS)::CMSSFile
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      CHARACTER(Len=16)::ScrChar

      IF(CMSSFile.eq.'XMS') THEN
        CALL ReadMat('ROT_VEC',ScrChar,RotMat,lroots,lroots,7,16,'T')
      ELSE
        CALL ReadMat(CMSSFile ,ScrChar,RotMat,lroots,lroots,LenCMSS,16, &
     &               'T')
      END IF
      RETURN
      End Subroutine
