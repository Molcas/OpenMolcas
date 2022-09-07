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

subroutine Mk_nacCD_Shells(kCnttp,lCnttp)
!***********************************************************************
!                                                                      *
!     Flowchart                                                        *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: kCnttp, lCnttp

!***********************************************************************
!                                                                      *
!     Produce product basis                                            *
!                                                                      *
!***********************************************************************

! Loop over iAng and jAng
! Pick up exponents and coefficients

!***********************************************************************
!                                                                      *
!     aCD: Transform to auxiliary basis set                            *
!                                                                      *
!***********************************************************************

! aCD subroutine:
! Still looping over yang and jAng
! Cholesky decomposition (only contracted + with a threshold)
! End loop

!***********************************************************************
!                                                                      *
!     naCD: Make the set nodeless                                      *
!***********************************************************************

! Loop until yAng=nTest*2 (naCD shells)
! Loop over aCD shells (iAng and jAng)
! Put exponents at the "right place" i.e. grouping the similar shell's exponents
! End aCD loop
! Fitting of coefficients to their own original ones
! Projection of the space spanned by the auxiliary basis set formed
!    by aCD methods to the space spanned by the new naCD
!    (calculating AA and AB matrix elements)

!***********************************************************************
!                                                                      *
!     nacCD: Contract the set                                          *
!                                                                      *
!***********************************************************************

! Still looping over yAng
! Fit the primitives
! Reduce the exponents: eliminate contraction functions below the CD threshold
! Refit the contraction coefficients
! End loop

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(kCnttp)
  call Unused_integer(lCnttp)
end if

end subroutine Mk_nacCD_Shells
