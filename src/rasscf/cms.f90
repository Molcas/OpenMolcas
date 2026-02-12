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

implicit none

logical CMSNotConverged
logical CMSGiveOpt
real*8 CMSThres
real*8, dimension(:), allocatable :: RGD
integer iCMSOpt
logical PosHess, BigQaaGrad, NeedMoreStep
integer nPosHess, nCMSScale
real*8 LargestQaaGrad
character*128 cmsguessfile

end module CMS
