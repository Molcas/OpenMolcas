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
      Subroutine MINV(ARRAY,ARRINV,ISING,DET,NDIM)
      Implicit Real*8 (a-h,o-z)
      Real*8 ARRAY(NDIM,NDIM), ARRINV(NDIM,NDIM)
#include "WrkSpc.fh"
!
      Call Allocate_Work(ipA,NDIM**2)
      Call Allocate_Work(ipB,NDIM**2)
      Call Allocate_Work(ipBUF,NDIM)
      Call Allocate_iWork(IPIV,NDIM)
      Call Allocate_iWork(JPIV,NDIM)
!
      Call MINV_INNER(ARRAY,ARRINV,ISING,DET,NDIM,Work(ipA),            &
     &                Work(ipBUF),Work(ipB),iWork(IPIV),                &
     &                iWork(JPIV))
!
      Call Free_iWork(JPIV)
      Call Free_iWork(IPIV)
      Call Free_Work(ipBUF)
      Call Free_Work(ipB)
      Call Free_Work(ipA)
!
      Return
      End
