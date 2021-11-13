!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module mula_global

use Constants, only: Two, Pi, cLight
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: maxMax_n = 64, MaxNumAt = 50
real(kind=wp), parameter :: hbarcm = 1.0e-2/(Two*Pi*cLight)

integer(kind=iwp) :: inpUnit, mdim1, mdim2, ndata, ndim1, ndim2, ngdim, nPolyTerm = 0, nvar
real(kind=wp) :: cmstart, cmend, energy1, energy2, LifeTime, TranDip(3)
logical(kind=iwp) :: broadplot, Huge_Print, OscStr, plotwindow, Use_cm, Use_nm, VibModPlot, WriteVibLevels
integer(kind=iwp), allocatable :: ipow(:,:), m_plot(:), n_plot(:), NormModes(:)
real(kind=wp), allocatable :: AtCoord1(:,:), AtCoord2(:,:), Hess1(:,:), Hess2(:,:), t_dipin1(:), t_dipin2(:), TranDipGrad(:,:), &
                              var(:,:), yin1(:), yin2(:)

public :: AtCoord1, AtCoord2, broadplot, cmstart, cmend, energy1, energy2, hbarcm, Hess1, Hess2, Huge_Print, inpUnit, ipow, &
          LifeTime, m_plot, maxMax_n, MaxNumAt, mdim1, mdim2, n_plot, ndata, ndim1, ndim2, ngdim, NormModes, nPolyTerm, nvar, &
          OscStr, plotwindow, t_dipin1, t_dipin2, TranDip, TranDipGrad, Use_cm, Use_nm, var, VibModPlot, WriteVibLevels, yin1, yin2

end module mula_global
