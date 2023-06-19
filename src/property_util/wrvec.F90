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
!  WRVEC
!
!> @brief
!>   A routine to write MO coefficients, occupation numbers, one-electron energies and type index information to ``INPORB`` file
!> @author V. Veryazov
!>
!> @details
!> New version of ::wrvec routine.
!> ::WRVEC is a wrapper to ::WRVEC_, which writes UHF
!> information to ``INPORB`` file.
!>
!> \p Label defines the type of information to write to ``INPORB`` file.
!> Valid targets are:
!> ``C``---CMO, ``O``---OCC, ``E``---EORB, ``I``---INDT, ``A``---Append Index, ``K``---Coordinates, ``B``---Basis section
!>
!> Example: Write CMO coeff. for RHF:
!>
!> \code
!> call WrVec('INPORB',Lu,'C',NSYM,NBAS,NBAS,CMO,Dummy,Dummy,iDummy,Title)
!> \endcode
!>
!> @param[in] FName File name
!> @param[in] LU_   Unit number
!> @param[in] LABEL Task
!> @param[in] NSYM  N symmetries
!> @param[in] NBAS  N basis functions
!> @param[in] NORB  N orbitals
!> @param[in] CMO   MO coefficients
!> @param[in] OCC   Occupations
!> @param[in] EORB  One electron energies
!> @param[in] INDT  Type Index information
!> @param[in] TITLE Title of orbitals
!***********************************************************************

subroutine WRVEC(FName,LU_,LABEL,NSYM,NBAS,NORB,CMO,OCC,EORB,INDT,TITLE)

use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: FName, LABEL, TITLE
integer(kind=iwp), intent(in) :: LU_, NSYM, NBAS(NSYM), NORB(NSYM), INDT(7,8)
real(kind=wp), intent(in) :: CMO(*), OCC(*), EORB(*)
real(kind=wp) :: vDum(2)

call WrVec_(FName,LU_,LABEL,0,NSYM,NBAS,NORB,CMO,vDum,OCC,vDum,EORB,vDum,INDT,TITLE,0)

return

end subroutine WRVEC
