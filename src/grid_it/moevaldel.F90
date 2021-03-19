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
! Copyright (C) 1995, Roland Lindh                                     *
!               2000, Valera Veryazov                                  *
!               2014, Thomas Dresselhaus                               *
!***********************************************************************
      Subroutine MOEvalDel(MOValueD,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt)

      Implicit Real*8 (A-H,O-Z)
      Real*8 Ccoor(3,nCoor),MOValueD(4*nCoor*nMOs),CMOs(nCMO)
      Integer DoIt(nMOs),mAO
      integer nDrv

      mAO  = 4
      nDrv = 1

      Call MOEval(MOValueD,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)

!        IJ1=1+(I-1)*4
!        IJ2=2+(I-1)*4
!        IJ3=3+(I-1)*4
!        IJ4=4+(I-1)*4
!
!        MOValue(I)=Work(MOTmp-1+IJ1)
!        MOValueDX(I)=Work(MOTmp-1+IJ2)
!        MOValueDY(I)=Work(MOTmp-1+IJ3)
!        MOValueDZ(I)=Work(MOTmp-1+IJ4)
!      END DO

      Return
      End
