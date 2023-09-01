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
!  Cho_X_getVfull
!
!> @brief
!>   Reorder Cholesky vectors from reduced to full storage
!> @author F. Aquilante
!>
!> @details
!> This routine reorders Cholesky vectors from reduced to
!> full storage. For \p DoRead = ``.true.`` the vectors are read from
!> disk using array \p RedVec as scratch space, whereas for
!> \p DoRead = ``.false.`` the reduced vectors must be supplied in
!> array \p RedVec.
!>
!> Given a set of target arrays (\p Wab),
!> the routine reorders \p NUMV Cholesky
!> vectors of compound symmetry \p ISYM starting with
!> vector \p JVEC1 and returns them in the target arrays.
!> Each pointer should thereby point to a
!> location where the corresponding Cholesky
!> vector of a given unique symmetry pair
!> of indices has to be stored.
!>
!> - \p iSwap = ``0``: \f$ L(a,b,J) \f$ is returned (in LT-storage if \f$ \mathrm{sym}(a)=\mathrm{sym}(b) \f$)
!> - \p iSwap = ``1``: \f$ L(a,J,b) \f$ is returned
!> - \p iSwap = ``2``: \f$ L(a,J,b) \f$ is returned (in SQ-storage if \f$ \mathrm{sym}(a)=\mathrm{sym}(b) \f$)
!>
!> - \p iSkip(syma) = ``0``: skip the symmetry block \f$ a \f$. Any vector \f$ L_{ab} \f$ or \f$ L_{ba} \f$
!>                           with \c syma &times; \c symb = \p ISYM won't be returned in the target array
!>
!> - \p IREDC: reduced set in core at the moment of the call to the routine.
!>             Can be set to ``-1`` (= unknown or undefined) by the calling routine.
!>
!> @param[out]    irc     return code
!> @param[in,out] RedVec  vectors stored in reduced set(s) [\p DoRead option off]
!>                        or scratch space for reading reduced vectors [\p DoRead option on]
!> @param[in]     lRedVec size of the \p RedVec
!> @param[in]     IVEC1   first vector to read
!> @param[in]     NUMV    number of vectors to read starting from \p IVEC1
!> @param[in]     ISYM    compound symmetry of the Cholesky vectors
!> @param[in]     iSwap   type of the full storage for the returned Cholesky vectors
!> @param[in,out] IREDC   current reduced set in core (location ``3``)
!> @param[in,out] Wab     target arrays
!> @param[in]     iSkip   skipping parameters for each symmetry block \f$ (ab) \f$ of compound symmetry \p ISYM
!> @param[in]     DoRead  flag for reading reduced vectors from disk
!***********************************************************************

subroutine Cho_X_getVfull(irc,RedVec,lRedVec,IVEC1,NUMV,ISYM,iSwap,IREDC,Wab,iSkip,DoRead)

use Symmetry_Info, only: Mul
use Cholesky, only: nBas, nSym
use Data_structures, only: Map_to_SBA, SBA_Type
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lRedVec, IVEC1, NUMV, ISYM, iSwap, iSkip(*)
real(kind=wp), intent(inout) :: RedVec(lRedVec)
integer(kind=iwp), intent(inout) :: IREDC
type(SBA_Type), intent(inout) :: Wab
logical(kind=iwp), intent(in) :: DoRead
integer(kind=iwp) :: i, iOff, ipChoV(8), ipVec(8), iSymp, iSymq, IVEC2, j, JNUM, JVEC1, jVref, MUSED, MXUSD, n2BSF(8,8), nnBSF(8,8)
integer(kind=iwp), external :: ip_of_Work

MXUSD = 0
MUSED = 0

ipChoV(1:nSym) = -6666
call Map_to_SBA(Wab,ipChoV)
! Get pointers relative to Wab%A0
iOff = ip_of_Work(Wab%A0(1))
ipChoV(1:nSym) = ipChoV(1:nSym)-iOff+1

ipVec(1:nSym) = ipChoV(1:nSym)

call set_nnBSF(nSym,nBas,nnBSF,n2BSF)

if (iSwap == 0) then

  do iSymq=1,nSym
    iSymp = mul(iSym,iSymq)
    if ((nnBSF(iSymp,iSymq) > 0) .and. ((iSymp >= iSymq) .and. (iSkip(iSymp) /= 0))) Wab%SB(iSymp)%A2(:,:) = Zero
  end do

else if ((iSwap == 1) .or. (iSwap == 2)) then

  do iSymq=1,nSym
    iSymp = mul(iSym,iSymq)
    if ((n2BSF(iSymp,iSymq) > 0) .and. ((iSymp >= iSymq) .and. (iSkip(iSymp) /= 0))) Wab%SB(iSymp)%A2(:,:) = Zero
  end do

else

  write(u6,*) 'Wrong parameter! iSwap= ',iSwap
  irc = 66
  return

end if

if (DoRead) then

  JVEC1 = IVEC1
  IVEC2 = JVEC1+NUMV-1

  do while (jVec1 <= iVec2)
    call CHO_VECRD(RedVec,lRedVec,JVEC1,IVEC2,ISYM,JNUM,IREDC,MUSED)
    MXUSD = max(MXUSD,MUSED)

    if ((JNUM <= 0) .or. (JNUM > IVEC2-JVEC1+1)) then
      irc = 77
      return
    end if

    jVref = JVEC1-IVEC1+1
    call cho_Reordr(irc,RedVec,lRedVec,jVref,JVEC1,JNUM,NUMV,ISYM,IREDC,iSwap,ipVec,Wab%A0,iSkip)

    if (irc /= 0) return

    jVec1 = jVec1+JNUM

    do i=1,nSym
      j = mul(i,ISYM)
      if ((j >= i) .and. (iSkip(j) /= 0)) then
        if (iSwap == 0) then
          ipVec(j) = ipVec(j)+nnBSF(j,i)*JNUM
        else if (iSwap == 1) then
          ipVec(j) = ipChoV(j)
        else if (iSwap == 2) then
          ipVec(j) = ipVec(j)+n2BSF(j,i)*JNUM
        end if
      end if
    end do

  end do  ! end the while loop

else ! only reorder

  JNUM = NUMV

  call cho_Reordr(irc,RedVec,lRedVec,1,IVEC1,JNUM,NUMV,ISYM,IREDC,iSwap,ipVec,Wab%A0,iSkip)

  if (irc /= 0) return

end if

irc = 0

return

end subroutine Cho_X_getVfull
