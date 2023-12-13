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
! Copyright (C) 2004,2005, Giovanni Ghigo                              *
!***********************************************************************
!  Def_TCVx
!
!> @brief
!>   The routine defines which Transformed Cholesky (TCVx) to generate
!>   setting ``.True.`` the logical matrix \c TCVXist(iType,iSym,jSym)
!> @author Giovanni Ghigo
!>
!> @details
!> - \c iType = ``1``: TCVA
!> - \c iType = ``2``: TCVB
!> - \c iType = ``3``: TCVC
!> - \c iType = ``4``: TCVD
!> - \c iType = ``5``: TCVE
!> - \c iType = ``6``: TCVF
!> - \c iType = ``7``: TCVG
!>
!> @param[in] iSym Symmetry(``i``) of the Cholesky Full Vector
!> @param[in] jSym Symmetry(``j``) of the Cholesky Full Vector
!***********************************************************************

subroutine Def_TCVx(iSym,jSym)
!***********************************************************************
! Author  :  Giovanni Ghigo                                            *
!            Lund University, Sweden                                   *
! Written :  October 2004                                              *
! Modified:  January 2005                                              *
!----------------------------------------------------------------------*
! Define which Transformed Cholesky Full Vectors (TCVx) to generate.   *
! TCVXist(iType,iSym,jSym) is .True. if the TCVx must be generated.    *
! iType(x):  1=A, 2=B, 3=C, 4=D, 5=E, 6=F, 7=G                         *
!***********************************************************************

use Cho_Tra, only: DoTCVA, nAsh, nIsh, nSsh, TCVXist
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, jSym

if (nIsh(jSym) > 0) then
  if (DoTCVA) then
    if (nIsh(iSym) > 0) then
      TCVXist(1,iSym,jSym) = .true.
      TCVXist(1,jSym,iSym) = .true.  ! Aji = T(Aij) !
    end if
    if (nAsh(iSym) > 0) then
      TCVXist(2,iSym,jSym) = .true.
      TCVXist(7,jSym,iSym) = .true.
    end if
  end if
  if (nSsh(iSym) > 0) then
    TCVXist(3,iSym,jSym) = .true.
  end if
end if

if ((nAsh(jSym) > 0) .and. DoTCVA) then
  if ((nIsh(iSym) > 0) .and. (iSym /= jSym)) then
    TCVXist(2,jSym,iSym) = .true.
    TCVXist(7,iSym,jSym) = .true.
  end if
  if (nAsh(iSym) > 0) then
    TCVXist(4,iSym,jSym) = .true.
    TCVXist(4,jSym,iSym) = .true.  ! Dji = T(Dij) !
  end if
  if (nSsh(iSym) > 0) then
    TCVXist(5,iSym,jSym) = .true.
  end if
end if

if ((nSsh(jSym) > 0) .and. (iSym /= jSym)) then
  if (nIsh(iSym) > 0) then
    TCVXist(3,jSym,iSym) = .true.
  end if
  if ((nAsh(iSym) > 0) .and. DoTCVA) then
    TCVXist(5,jSym,iSym) = .true.
  end if
end if

if ((nSsh(jSym) > 0) .and. (nSsh(iSym) > 0) .and. DoTCVA) then
  TCVXist(6,iSym,jSym) = .true.
end if

return

end subroutine Def_TCVx
