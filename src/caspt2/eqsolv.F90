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

use caspt2_module, only: FockType, G1SecIn
use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXCASE = 13
integer(kind=iwp) :: IDBMAT(8,MXCASE), IDSMAT(8,MXCASE), IDSTMAT(8,MXCASE), IDTMAT(8,MXCASE), IRHS, IVECC, IVECC2, IVECR, IVECW, &
                     IVECX, LLIST(8,8,17), MODVEC(8,MXCASE), MXSCT, NLIST(8,8,17), NLSTOT
! Coupling cases of the off-diagonal blocks of H0 (KOD of SGM and OFFDNS).
! Use Cpl_Code to decide whether a coupling enters a given quantity.
integer(kind=iwp), parameter :: IFCOUP(MXCASE,MXCASE) = reshape( &
                     [ 0, 1, 2, 0, 3, 4, 5, 0, 0, 0, 0, 0, 0, &
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

public :: Cpl_Code, IDBMAT, IDSMAT, IDSTMat, IDTMAT, IFCOUP, iRHS, iVecC, iVecC2, iVecR, iVecW, iVecX, LLIST, MODVEC, MXSCT, &
          NLIST, NLSTOT

contains

function Cpl_Code(iCase2,iCase1,H0) result(KOD)
! Coupling code of the cases iCase2 > iCase1, 0 if they do not couple.
! H0 = .true.: in H0 and its derivatives (sigma, gradients). With
! FOCKtype G1 the inactive-secondary block of the Fock matrix is left out
! of H0, as in K. Andersson, Theor. Chim. Acta 91, 31 (1995); it enters
! only the couplings A-E (KOD 4, 5), C-G (11, 12) and D-H (23, 24).
! The keyword G1SEcin keeps it (H0 = full projected Fock operator).
! H0 = .false.: all couplings (e.g. the expectation value <PSI1|E_pq|PSI1>).

  integer(kind=iwp) :: KOD
  integer(kind=iwp), intent(in) :: iCase2, iCase1
  logical(kind=iwp), intent(in) :: H0

  KOD = IFCOUP(iCase2,iCase1)
  if (H0 .and. (FockType == 'G1') .and. (.not. G1SecIn)) then
    select case (KOD)
      case (4,5,11,12,23,24)
        KOD = 0
    end select
  end if

end function Cpl_Code

end module EQSOLV
