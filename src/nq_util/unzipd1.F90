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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
      Subroutine UnzipD1(D1Unzip,D1MO,nD1MO)
      use nq_Info

!*****Input
      INTEGER nD1MO
      Real*8,DIMENSION(nD1MO)::D1MO
!*****Output
      Real*8,DIMENSION(NASHT**2)::D1Unzip
!*****Intermediate
      INTEGER iv,ix,iLoc1,iLoc2,iLoc3

      CALL FZero(D1Unzip,NASHT**2)
      DO iv=1,NASHT
       Do ix=1,iv-1
        iLoc1=(iv-1)*NASHT+ix
        iLoc2=(ix-1)*NASHT+iv
        iLoc3=(iv-1)*iv/2+ix
        D1Unzip(iLoc1)=0.5d0*D1MO(iLoc3)
        D1Unzip(iLoc2)=D1Unzip(iLoc1)
       End Do
       ix=iv
       iLoc1=(iv-1)*NASHT+ix
       iLoc3=(iv+1)*iv/2
       D1Unzip(iLoc1)=0.5d0*D1MO(iLoc3)
      END DO

      RETURN
      End Subroutine
