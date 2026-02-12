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
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
!*****************************************************************

subroutine CalcDDg(DDg,GD,Dg,nDDg,nGD,lRoots2,NAC2)
!*****************************************************************
!     Gtuvx : two-electron integral, g_tuvx                      *
!     GD    : "generalized 1-e density matrix"                   *
!              GD^KL: transition density matrix from L to K      *
!              GD^KK: density matrix for state K                 *
!     Dg    :  sum_{vx}{GD^KL_vx * g_tuvx}                       *
!     In GDorbit and Dgorbit, the leading index is orbital index;*
!     In GDstate and Dgstate, the leading index is state index.  *
!                                                                *
!     DDg   :  sum_{tuvx}{GD^KL_tu * GD^MN_vx * g_tuvx}          *
!      namely, sum_{tu}{GD^KL_tu * Dg^MN_tu}                     *
!*****************************************************************

integer nDDg, nGD, lRoots2, NAC2
real*8 DDg(nDDg), GD(nGD), Dg(nGD)

call DGEMM_('T','N',lRoots2,lRoots2,NAC2,1.0d0,Dg,NAC2,GD,NAC2,0.0d0,DDg,lRoots2)

return

end subroutine CalcDDg
