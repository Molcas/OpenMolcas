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
Module Disp_MCLR
#include "LenIn.fh"
Integer, Parameter:: LENIN6=LENIN+6
Integer lDisp(8)
Integer, Parameter,Private:: Mxdccc=500
Character ChDisp(Mxdccc*3)*(LENIN+6)
Character(LEN=8) SwLbl(Mxdccc)
integer   dspvec(mxdccc),nhess
Private LenIN,LENIN6
END Module Disp_MCLR
