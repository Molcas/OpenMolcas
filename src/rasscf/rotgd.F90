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

! This file contains subroutines relating to generalized 1-e
! density matrix (GD) called in CMSNewton, including
! CalcGD:       calculating GD with lucia.
! CalcDg:       calculating Dg matrix, namely sum_{vx}{GD^KL_vx * g_tuvx}
! RotGD:        GD^KL_tu = sum_{MN}{U^KM * U^LN * GD^MN_tu}
! TransposeMat: transform GD^KL_tu from leading with state indices
!               to leading with orbital indices in mode 1, and vice
!               versa in mode 2.



      Subroutine RotGD(GD,R,nGD,lRoots,NAC2)
      use CMS, only: RGD
      Implicit None
      INTEGER nGD,lRoots,NAC2
      Real*8 GD(nGD),R(lRoots**2)

      INTEGER iNAC2,iLoc,lRoots2
!      Real*8 RGD(lRoots**2),RGDR(lRoots**2)

      lRoots2=lRoots**2


!      write(6,*) 'rotation matrix in RotGD'
!      CALL RecPrt(' ',' ',R,lRoots,lRoots)
!
!      write(6,*) 'GD matrix after rotation'
!      CALL RecPrt(' ',' ',GD,lRoots2,NAC2)


      DO iNAC2=1,NAC2
       iLoc=(iNAC2-1)*lRoots2+1
       CALL DGEMM_('T','N',lRoots,lRoots,lRoots,                        &
     &               1.0d0,R       ,lRoots,GD(iLoc),lRoots,             &
     &               0.0d0,RGD     ,lRoots)
       CALL DGEMM_('N','N',lRoots,lRoots,lRoots,                        &
     &               1.0d0,RGD     ,lRoots,R       ,lRoots,             &
     &               0.0d0,GD(iLoc),lRoots)
      END DO

!      write(6,*) 'GD matrix after rotation'
!      CALL RecPrt(' ',' ',GD,lRoots2,NAC2)

      RETURN
      End Subroutine
