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

subroutine CMSHeader(CMSSFile,LenCMSS)

use CMS, only: iCMSOpt, CMSGuessFile
use rasscf_global, only: CMSThreshold, iCMSIterMin, iCMSIterMax, lRoots

implicit none
integer LenCMSS
character(len=LenCMSS) :: CMSSFile
#include "warnings.h"

write(6,*)
write(6,*)
write(6,'(4X,A35)') 'CMS INTERMEDIATE-STATE OPTIMIZATION'
if (CMSSFile == 'XMS') then
  write(6,'(5X,A12,8X,A25)') 'START MATRIX','XMS INTERMEDIATE STATES'
else
  write(6,'(5X,A12,8X,A25)') 'START MATRIX',CMSGuessFile
end if
if (iCMSOpt == 1) then
  write(6,'(5X,A8,12X,A25)') 'OPT ALGO','NEWTON'
else if (iCMSOpt == 2) then
  write(6,'(5X,A8,12X,A25)') 'OPT ALGO','JACOBI'
end if
write(6,'(5X,A15,5X,16X,ES9.2E2)') 'Q_a-a THRESHOLD',CMSThreshold
if (iCMSOpt == 1) write(6,'(5X,A15,5X,16X,ES9.2E2)') 'GRAD  THRESHOLD',CMSThreshold*1.0d-2
write(6,'(5X,A10,10X,I25)') 'MAX CYCLES',ICMSIterMax
write(6,'(5X,A10,10X,I25)') 'MIN CYCLES',ICMSIterMin
write(6,*) repeat('=',71)
if (iCMSOpt == 2) then
  if (lRoots > 2) then
    write(6,'(4X,A8,2X,2(A16,11X))') 'Cycle','Q_a-a','Difference'
  else
    write(6,'(4X,A8,2X,A18,6X,A8,12X,A12)') 'Cycle','Rot. Angle (deg.)','Q_a-a','Q_a-a Diff.'
  end if
else
  write(6,'(6X,A5,7X,A5,8X,A10,2X,A6,5X,A7,4X,A4)') 'Cycle','Q_a-a','Difference','# Pos.','Largest','Step'
  write(6,'(43X,A7,4X,A8,3X,A6)') 'Hessian','Gradient','Scaled'
end if
write(6,*) repeat('-',71)

return

end subroutine CMSHeader
