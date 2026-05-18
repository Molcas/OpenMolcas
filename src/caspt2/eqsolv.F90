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

module EQSOLV

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXCASE = 13
integer(kind=iwp) :: IDBMAT(8,MXCASE), IDSMAT(8,MXCASE), IDSTMAT(8,MXCASE), IDTMAT(8,MXCASE), IRHS, IVECC, IVECC2, IVECR, IVECW, &
                     IVECX, LLIST(8,8,17), MODVEC(8,MXCASE), MXSCT, NLIST(8,8,17), NLSTOT
! This is not a parameter, because the G1 correction modifies it
integer(kind=iwp) :: IFCOUP(MXCASE,MXCASE) = reshape( &
                     [0, 1, 2, 0, 3, 4, 5, 0, 0, 0, 0, 0, 0, &
                      0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, &
                      0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, &
                      0, 0, 0, 0, 8, 0, 0, 9,10,11,12, 0, 0, &
                      0, 0, 0, 0, 0,13,14, 0, 0,15,16,23,24, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,17, 0, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,18, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0,19, 0, 0, 0, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,20, 0, 0, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,21, 0, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,22, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], &
                     [MXCASE,MXCASE])

public :: IDBMAT, IDSMAT, IDSTMat, IDTMAT, IFCOUP, iRHS, iVecC, iVecC2, iVecR, iVecW, iVecX, LLIST, MODVEC, MXSCT, NLIST, NLSTOT

end module EQSOLV
