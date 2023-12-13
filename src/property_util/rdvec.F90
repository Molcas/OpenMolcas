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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  RDVEC
!
!> @brief
!>   A routine to read MO coefficients, occupation numbers, one-electron energies and type index information from ``INPORB`` file
!> @author V. Veryazov
!>
!> @details
!> New version of ::rdvec routine.
!> ::RDVEC is a wrapper to ::RDVEC_, which read UHF
!> information from ``INPORB`` file.
!>
!> \p Label defines the type of information to read from ``INPORB`` file
!> Valid targets are: ``C``---CMO, ``O``---OCC, ``E``---EORB, ``I``---INDT
!> ``A``---alpha values, ``B``---beta values
!>
!> ::RdVec checks that \p NBAS / \p NORB information is consistent,
!> and reacts according to \p iWarn. ``0``: No checks for \p NBAS / \p NORB;
!> ``1``: Print error message; ``2``: ::Abend.
!>
!> Example: Get CMO coeff. and OCC for RHF:
!>
!> \code
!> call RdVec('INPORB',Lu,'CO',NSYM,NBAS,NBAS,CMO,OCC,Dummy,iDummy,Title,0,iErr)
!> \endcode
!>
!> @param[in]  FName File name
!> @param[in]  LU_   Unit number
!> @param[in]  LABEL Task
!> @param[in]  NSYM  N symmetries
!> @param[in]  NBAS  N basis functions
!> @param[in]  NORB  N orbitals
!> @param[out] CMO   MO coefficients
!> @param[out] OCC   Occupations
!> @param[out] EORB  One electron energies
!> @param[out] INDT  Type Index information
!> @param[out] TITLE Title of orbitals
!> @param[in]  IWARN Warning level
!> @param[out] IERR  Return code
!***********************************************************************

subroutine RDVEC(FName,LU_,LABEL,NSYM,NBAS,NORB,CMO,OCC,EORB,INDT,TITLE,iWarn,iErr)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
character(len=*), intent(in) :: FName, LABEL
integer(kind=iwp), intent(in) :: LU_, NSYM, NBAS(NSYM), NORB(NSYM), iWarn
real(kind=wp), intent(_OUT_) :: CMO(*), OCC(*), EORB(*)
integer(kind=iwp), intent(_OUT_) :: INDT(*)
character(len=*), intent(out) :: TITLE
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: iWFType
real(kind=wp) :: vDum(2)

call RdVec_(FName,LU_,LABEL,0,NSYM,NBAS,NORB,CMO,vDum,OCC,vDum,EORB,vDum,INDT,TITLE,iWarn,iErr,iWFtype)

return

end subroutine RDVEC
