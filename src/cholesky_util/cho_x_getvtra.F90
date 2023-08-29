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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!  Cho_X_getVtra
!
!> @brief
!>    This routine performs a half-MO-transformation of Cholesky vectors stored in reduced storage
!> @author F. Aquilante
!>
!> @details
!> This routine performs a half-MO-transformation of Cholesky vectors stored in reduced
!> storage. For \p DoRead = ``.true.`` the vectors are read from
!> disk using array \p RedVec as scratch space, whereas for
!> \p DoRead = ``.false.`` the reduced vectors must be supplied in
!> array \p RedVec.
!>
!> Given (\p ChoT),the target arrays of type SBA_Type,
!> the routine performs a half-MO-transformation of \p NUMV Cholesky
!> vectors of compound symmetry \p ISYM starting with
!> vector \p IVEC1 and returns them in the target arrays.
!>
!> - \p iSwap = ``0``: \f$ L(k,b,J) \f$ is returned
!> - \p iSwap = ``1``: \f$ L(a,k,J) \f$ is returned
!> - \p iSwap = ``2``: \f$ L(k,J,b) \f$ is returned
!> - \p iSwap = ``3``: \f$ L(a,J,k) \f$ is returned
!>
!>
!> - \p IREDC: reduced set in core at the moment of the call to the routine.
!>             Can be set to ``-1`` (= unknown or undefined) by the calling routine.
!>
!> @param[out]    irc     return code
!> @param[in,out] RedVec  Vectors stored in reduced set(s) [\p DoRead option off] or scratch space for reading reduced vectors [\p DoRead option on]
!> @param[in]     lRedVec size of the \p RedVec
!> @param[in]     IVEC1   first vector to read
!> @param[in]     NUMV    number of vectors to transform starting from \p IVEC1
!> @param[in]     ISYM    compound symmetry of the Cholesky vectors
!> @param[in]     iSwap   type of the full storage for the half transformed Cholesky vectors
!> @param[in,out] IREDC   reduced set in core
!> @param[in]     nDen    total number of densities to which MOs refer
!> @param[in]     kDen    first density for which the MO transformation has to be performed
!> @param[in]     MOs     the MOs coefficients stored in the data type DSBA_Type, i.e. symmetry blocked.
!> @param[out]    ChoT    the half transformed vectors, symmetry blocked as type SBA_Type
!> @param[in]     DoRead  flag for reading the reduced vectors
!***********************************************************************

subroutine Cho_X_getVtra(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,iSwap,IREDC,nDen,kDen,MOs,ChoT,DoRead)

use Data_Structures, only: DSBA_Type, SBA_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lRedVec, IVEC1, NUMV, ISYM, iSwap, nDen, kDen
real(kind=wp), intent(inout) :: RedVec(lRedVec)
integer(kind=iwp), intent(inout) :: IREDC
type(DSBA_Type), intent(in) :: MOs(nDen)
type(SBA_Type), intent(_OUT_) :: ChoT(nDen)
logical(kind=iwp), intent(in) :: DoRead
integer(kind=iwp) :: IVEC2, jDen, JNUM, JVEC1, JVREF, MUSED, MXUSD

MXUSD = 0
MUSED = 0
! zeroing the target arrays
!--------------------------
do jDen=kDen,nDen
  ChoT(jDen)%A0(:) = Zero
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoRead) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  JVEC1 = IVEC1         ! Absolute starting index
  IVEC2 = JVEC1+NUMV-1  ! Absolute ending index

  do while (jVec1 <= iVec2)
    call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)
    MXUSD = max(MXUSD,MUSED)

    if ((JNUM <= 0) .or. (JNUM > IVEC2-JVEC1+1)) then
      irc = 77
      return
    end if

    JVREF = JVEC1-IVEC1+1 ! Relative index
    call cho_vTra(irc,RedVec,lRedVec,JVREF,JVEC1,JNUM,NUMV,ISYM,IREDC,iSwap,nDen,kDen,MOs,ChoT)
    if (irc /= 0) return

    JVEC1 = jVec1+JNUM

  end do  ! end the while loop
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else ! only MO transformation
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  JNUM = NUMV
  JVREF = 1
  call cho_vTra(irc,RedVec,lRedVec,JVREF,IVEC1,JNUM,NUMV,ISYM,IREDC,iSwap,nDen,kDen,MOs,ChoT)
  if (irc /= 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
irc = 0

return

end subroutine Cho_X_getVtra
