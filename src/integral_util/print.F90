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
Module Print
private
#include "print.fh"
!     Integer, parameter:: nRout=1024
!     Integer :: nPrint(nRout) = [5,i=1,nRout]
!     Logical :: Show=.False.
!     Integer :: icolorize=0
public nPrint, Show, iColorize
End Module Print
