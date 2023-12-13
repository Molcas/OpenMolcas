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

Module CMS
logical CMSNotConverged
Logical CMSGiveOpt
Real*8  CMSThres
Real*8,DIMENSION(:),Allocatable:: RGD
INTEGER iCMSOpt
Logical PosHess,BigQaaGrad,NeedMoreStep
INTEGER nPosHess,nCMSScale
Real*8  LargestQaaGrad
CHARACTER*128 cmsguessfile
End Module CMS

