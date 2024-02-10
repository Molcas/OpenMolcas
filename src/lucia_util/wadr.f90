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
Module wadr
Private
Real*8, Allocatable, Public:: TUVX(:), FockOcc(:), DSPN(:), DMAT(:), PMAT(:), PA(:)
! These arrays are used for the SXCTL part of the code.
Real*8, Allocatable, Public:: DIA(:), SXN(:), BM(:), F1(:), F2(:), SXG(:), SXHD(:), SXH(:)
INTEGER, Public:: NLX
! These arrays are used for the TRACTL2 part of the code.
Real*8, Allocatable, Public:: FI(:), FA(:), D1I(:), D1A(:), OccN(:)
INTEGER, Public:: nPWXY
End Module wadr
