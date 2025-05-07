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
! Copyright (C) 2012, Victor P. Vysotskiy                              *
!               2025, Ignacio Fdez. Galvan                             *
!***********************************************************************

integer function kind2goff(var)

#include "mama.fh"
character(len=4) var

kind2goff = 0
if (var == 'INTE') kind2goff = iofint
if (var == 'REAL') kind2goff = iofdbl
if (var == 'CHAR') kind2goff = iofchr

end function kind2goff
