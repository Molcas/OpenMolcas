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

subroutine RMVERT(SGS)
! Purpose: Remove vertices from a DRT table.

use gugx, only: SGStruct
use RasDef, only: nRas, nRsPrt, nRasEl
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
type(SGStruct), intent(inout) :: SGS
integer(kind=iwp) :: IC, ID, iRO, iSy, IV, L, Lev, N, NCHANGES, NLD, NV
logical(kind=iwp) :: Test
integer(kind=iwp), allocatable :: CONN(:), Lim(:)
integer(kind=iwp), parameter :: LTAB = 1, NTAB = 2

! Construct a restricted graph.
call mma_allocate(Lim,SGS%nLev,Label='Lim')
Lim(:) = 0
! Fill in the occupation limit table:
Lev = 0
do iRO=1,nRsPrt
  do iSy=1,SGS%nSym
    Lev = Lev+nRas(iSy,iRO)
  end do
  if (Lev > 0) Lim(Lev) = nRasEl(iRO)
end do

call mma_allocate(SGS%Ver,SGS%nVert0,Label='SGS%Ver')
call mma_allocate(CONN,SGS%nVert,Label='CONN')

! KILL VERTICES THAT DO NOT OBEY RESTRICTIONS.
do IV=1,SGS%nVert-1
  SGS%Ver(IV) = 1
  L = SGS%DRT0(IV,LTAB)
  N = SGS%DRT0(IV,NTAB)
  if (N < Lim(L)) SGS%Ver(IV) = 0
end do
SGS%Ver(SGS%nVert) = 1

NCHANGES = 1 ! Initiate first loop
do while (NCHANGES > 0)
  ! REMOVE ARCS HAVING A DEAD UPPER OR LOWER VERTEX.
  ! COUNT THE NUMBER OF ARCS REMOVED OR VERTICES KILLED.
  NCHANGES = 0
  do IV=1,SGS%nVert-1
    if (SGS%Ver(IV) == 0) then
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        if (ID > 0) then
          SGS%Down0(IV,IC) = 0
          NCHANGES = NCHANGES+1
        end if
      end do
    else
      NLD = 0
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        if (ID > 0) then
          if (SGS%Ver(ID) == 0) then
            SGS%Down0(IV,IC) = 0
            NCHANGES = NCHANGES+1
          else
            NLD = NLD+1
          end if
        end if
      end do
      if (NLD == 0) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES+1
      end if
    end if
  end do
  ! ALSO CHECK ON CONNECTIONS FROM ABOVE:
  CONN(:) = 0
  CONN(1) = SGS%Ver(1)
  do IV=1,SGS%nVert-1
    if (SGS%Ver(IV) == 1) then
      do IC=0,3
        ID = SGS%Down0(IV,IC)
        Test = ID > 0
        if (Test) Test = SGS%Ver(ID) == 1
        if (Test) CONN(ID) = 1
      end do
    end if
  end do
  do IV=1,SGS%nVert
    if ((SGS%Ver(IV) == 1) .and. (CONN(IV) == 0)) then
      SGS%Ver(IV) = 0
      NCHANGES = NCHANGES+1
    end if
  end do

end do

! IF NO CHANGES, THE REMAINING GRAPH IS VALID.
! EVERY VERTEX OBEYS THE RESTRICTIONS. EVERY VERTEX IS
! CONNECTED ABOVE AND BELOW (EXCEPTING THE TOP AND BOTTOM)
! TO OTHER CONFORMING VERTICES.
! THE PROCEDURE IS GUARANTEED TO FIND A STABLE SOLUTIONS,
! SINCE EACH ITERATION REMOVES ARCS AND/OR VERTICES FROM THE
! FINITE NUMBER WE STARTED WITH.

if (SGS%Ver(1) == 0) then
  write(u6,*) 'RASSI/RMVERT: Too severe restrictions.'
  write(u6,*) 'Not one single configuration is left.'
  call ABEND()
end if

NV = 0
do IV=1,SGS%nVert
  if (SGS%Ver(IV) == 1) then
    NV = NV+1
    SGS%Ver(IV) = NV
  end if
end do
SGS%nVert = NV

call mma_deallocate(CONN)
call mma_deallocate(Lim)

end subroutine RMVERT
