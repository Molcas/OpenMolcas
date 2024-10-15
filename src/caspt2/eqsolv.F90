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
Module EQSOLV
INTEGER, PARAMETER, Private:: MXCASE=13
INTEGER, PARAMETER :: MXVEC=6, MXBLK=40*256*256

INTEGER IDSMAT(8,MXCASE),IDBMAT(8,MXCASE),IDTMAT(8,MXCASE),           &
        IDSTMAT(8,MXCASE),MODVEC(8,MXCASE),NLIST(8,8,17),             &
        LLIST(8,8,17),NLSTOT, MXSCT

INTEGER IRHS,IVECX,IVECR,IVECC,IVECC2,IVECW
END Module EQSOLV
