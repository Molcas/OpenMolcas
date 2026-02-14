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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 01, 2022, created this file.               *
! ****************************************************************

module CMS

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iCMSOpt, nCMSScale, nPosHess
real(kind=wp) :: CMSThres, LargestQaaGrad
logical(kind=iwp) :: BigQaaGrad, CMSGiveOpt, CMSNotConverged, NeedMoreStep, PosHess
character(len=128) :: cmsguessfile
real(kind=wp), allocatable :: RGD(:)

public :: BigQaaGrad, CMSGiveOpt, cmsguessfile, CMSNotConverged, CMSThres, iCMSOpt, LargestQaaGrad, nCMSScale, NeedMoreStep, &
          nPosHess, PosHess, RGD

end module CMS
