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
      Subroutine CalcFckS(FckO,GDMat,FckS)
      use rasscf_global, only: lRoots, nAc
      Implicit None


!*****Input
      Real*8,DIMENSION(NAC,NAC)::FckO
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat
!*****Output
      Real*8,DIMENSION(lRoots,lRoots)::FckS
!*****Auxiliary variables
      INTEGER IState,JState, iOrb, jOrb

      FckS(:,:)=0.0d0

      DO IState=1,lRoots
       Do JState=1,IState
        dO IOrb=1,NAC
         do JOrb=1,NAC
          FckS(IState,JState)=FckS(IState,JState)+FckO(IOrb,JOrb)*      &
     &GDMat(IState*(IState-1)/2+JState,IOrb,JOrb)
         end do
        eND DO
        FckS(JState,IState)=FckS(IState,JState)
        End Do
      END DO

!      CALL PrintMat('XMS_Mat','test',FckS,LRoots,LRoots,0,4,'N')

      END Subroutine CalcFckS
