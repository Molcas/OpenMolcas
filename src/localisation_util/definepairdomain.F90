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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  DefinePairDomain
!
!> @brief
!>   Define pair domains
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Pair domains are defined by the union of individual orbital
!> domains, \p iDomain(0:nAtom,nOcc), see SubRoutine ::DefineDomain.
!>
!> On exit, the contents of \p iPairDomain array are:
!>
!> - \p iPairDomain(0,ij): number of atoms in pair domain \c ij.
!> - \p iPairDomain(n,ij): id of atom \c n (\c n = ``1``, ``2``, ..., \p iPairDomain(0,ij)) in pair domain \c ij.
!>
!> The contents of \p Rmin array are:
!>
!> - \p Rmin(ij): minimum distance between atoms in pair \p ij.
!>
!> The contents of \p iClass array are (for \c i = ``1``, ``2`` , ``3``, ..., \p nRThr):
!>
!> - \p iClass(ij) = \c i-1 implies \p Rmin(ij) &le; \p RThr(i) and
!> - \p iClass(ij) = \p nRThr implies \p Rmin(ij) > \p RThr(nRThr)
!>
!> Note that if classification is not wanted, simply specify
!> \p nRThr = ``0`` (in which case arrays \p iClass and \p RThr are not referenced,
!> may be dummy arguments in the call to this routine).
!> The \p iDomain array should be set by SubRoutine ::DefineDomain
!> before calling this routine.
!>
!> Return codes:
!>
!> - \p irc = ``0``: all OK.
!>
!> (at the moment, this is the only possible return code, but that
!> might change.)
!>
!> @param[out] irc         Return code
!> @param[out] iPairDomain Pair domain definition
!> @param[out] iClass      Classification
!> @param[out] Rmin        Minimum interatomic distance in each pair
!> @param[in]  iDomain     Orbital domain definition
!> @param[in]  RThr        Thresholds for classification
!> @param[in]  Coord       Nuclear coordinates
!> @param[in]  nAtom       Number of atoms
!> @param[in]  nOcc        Number of orbitals for which domains are defined
!> @param[in]  nRThr       Number of thresholds for classification
!***********************************************************************

subroutine DefinePairDomain(irc,iPairDomain,iClass,Rmin,iDomain,RThr,Coord,nAtom,nOcc,nRThr)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nAtom, nOcc, iDomain(0:nAtom,nOcc), nRThr
integer(kind=iwp), intent(_OUT_) :: iPairDomain(0:nAtom,*), iClass(*)
real(kind=wp), intent(_OUT_) :: Rmin(*)
real(kind=wp), intent(in) :: RThr(*), Coord(3,nAtom)
integer(kind=iwp) :: i, iA, iAtom, iCount, ij, isThere, j, jA, jAtom, kAtom, lD, nnOcc
real(kind=wp) :: R
integer(kind=iwp), allocatable :: Union(:,:)

! Set return code.
! ----------------

irc = 0
if (nOcc < 2) return

! Set pair domains as union of individual domains: [ij]=[i]U[j] for i>=j.
! -----------------------------------------------------------------------

nnOcc = nOcc*(nOcc+1)/2
iPairDomain(:,1:nnOcc) = 0

call mma_allocate(Union,nAtom,nOcc,label='Union')

Union(:,:) = 0
do i=1,nOcc
  do iA=1,iDomain(0,i)
    iAtom = iDomain(iA,i)
    Union(iAtom,i) = 1
  end do
end do

ij = 0
do j=1,nOcc
  ij = ij+1
  lD = iDomain(0,j)
  iPairDomain(0:lD,ij) = iDomain(0:lD,j) ! case i=j
  do i=j+1,nOcc
    iCount = 0
    ij = ij+1
    do kAtom=1,nAtom
      isThere = Union(kAtom,j)+Union(kAtom,i)
      if (isThere > 0) then
        iCount = iCount+1
        iPairDomain(iCount,ij) = kAtom
      end if
    end do
    iPairDomain(0,ij) = iCount
  end do
end do

call mma_deallocate(Union)

! Set min. distance between any two atoms in pairs of domains.
! ------------------------------------------------------------

ij = 0
do j=1,nOcc
  ij = ij+1
  Rmin(ij) = Zero ! case i=j
  do i=j+1,nOcc
    ij = ij+1
    Rmin(ij) = huge(Rmin)
    do jA=1,iDomain(0,j)
      jAtom = iDomain(jA,j)
      do iA=1,iDomain(0,i)
        iAtom = iDomain(iA,i)
        R = sqrt((Coord(1,iAtom)-Coord(1,jAtom))**2+(Coord(2,iAtom)-Coord(2,jAtom))**2+(Coord(3,iAtom)-Coord(3,jAtom))**2)
        Rmin(ij) = min(Rmin(ij),R)
      end do
    end do
  end do
end do

! Set classification for each pair domain.
! Note that the caller must have ordered the thresholds according to
! RThr(1) < RThr(2) < RThr(3) < ... < RThr(nRThr)
! ------------------------------------------------------------------

if (nRThr > 0) then
  iClass(1:nnOcc) = nRThr
  do ij=1,nnOcc
    i = 0
    do while (i < nRThr)
      i = i+1
      if (Rmin(ij) <= RThr(i)) then
        iClass(ij) = i-1
        i = nRThr ! break while loop
      end if
    end do
  end do
end if

end subroutine DefinePairDomain
