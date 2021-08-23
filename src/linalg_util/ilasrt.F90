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

subroutine ILASRT(ID,N,D,INFO)
! Variant of LAPACK's [SD]LASRT for sorting an integer array

use Definitions, only: iwp

implicit none
character, intent(in) :: ID
integer(kind=iwp), intent(in) :: N
integer(kind=iwp), intent(inout) :: D(N)
integer(kind=iwp), intent(out) :: INFO
!  =====================================================================
integer(kind=iwp) :: D1, D2, D3, DIR, DMNMX, ENDD, I, J, STACK(2,32), START, STKPNT, TMP
logical(kind=iwp), external :: LSAME
integer(kind=iwp), parameter :: slct = 20

! Test the input parameters.

INFO = 0
DIR = -1
if (LSAME(ID,'D')) then
  DIR = 0
else if (LSAME(ID,'I')) then
  DIR = 1
end if
if (DIR == -1) then
  INFO = -1
else if (N < 0) then
  INFO = -2
end if
if (INFO /= 0) then
  call XERBLA('ILASRT',-INFO)
  return
end if

! Quick return if possible

if (N <= 1) return

STKPNT = 1
STACK(1,1) = 1
STACK(2,1) = N
do
  START = STACK(1,STKPNT)
  ENDD = STACK(2,STKPNT)
  STKPNT = STKPNT-1
  if ((ENDD-START <= slct) .and. (ENDD-START > 0)) then

    ! Do Insertion sort on D( START:ENDD )

    if (DIR == 0) then

      ! Sort into decreasing order

      do I=START+1,ENDD
        do J=I,START+1,-1
          if (D(J) > D(J-1)) then
            DMNMX = D(J)
            D(J) = D(J-1)
            D(J-1) = DMNMX
          else
            exit
          end if
        end do
      end do

    else

      ! Sort into increasing order

      do I=START+1,ENDD
        do J=I,START+1,-1
          if (D(J) < D(J-1)) then
            DMNMX = D(J)
            D(J) = D(J-1)
            D(J-1) = DMNMX
          else
            exit
          end if
        end do
      end do

    end if

  else if (ENDD-START > slct) then

    ! Partition D( START:ENDD ) and stack parts, largest one first
    !
    ! Choose partition entry as median of 3

    D1 = D(START)
    D2 = D(ENDD)
    I = (START+ENDD)/2
    D3 = D(I)
    if (D1 < D2) then
      if (D3 < D1) then
        DMNMX = D1
      else if (D3 < D2) then
        DMNMX = D3
      else
        DMNMX = D2
      end if
    else
      if (D3 < D2) then
        DMNMX = D2
      else if (D3 < D1) then
        DMNMX = D3
      else
        DMNMX = D1
      end if
    end if

    if (DIR == 0) then

      ! Sort into decreasing order

      I = START-1
      J = ENDD+1
      do
        do
          J = J-1
          if (D(J) >= DMNMX) exit
        end do
        do
          I = I+1
          if (D(I) <= DMNMX) exit
        end do
        if (I < J) then
          TMP = D(I)
          D(I) = D(J)
          D(J) = TMP
        else
          exit
        end if
      end do
      if (J-START > ENDD-J-1) then
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = START
        STACK(2,STKPNT) = J
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = J+1
        STACK(2,STKPNT) = ENDD
      else
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = J+1
        STACK(2,STKPNT) = ENDD
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = START
        STACK(2,STKPNT) = J
      end if
    else

      ! Sort into increasing order

      I = START-1
      J = ENDD+1
      do
        do
          J = J-1
          if (D(J) <= DMNMX) exit
        end do
        do
          I = I+1
          if (D(I) >= DMNMX) exit
        end do
        if (I < J) then
          TMP = D(I)
          D(I) = D(J)
          D(J) = TMP
        else
          exit
        end if
      end do
      if (J-START > ENDD-J-1) then
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = START
        STACK(2,STKPNT) = J
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = J+1
        STACK(2,STKPNT) = ENDD
      else
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = J+1
        STACK(2,STKPNT) = ENDD
        STKPNT = STKPNT+1
        STACK(1,STKPNT) = START
        STACK(2,STKPNT) = J
      end if
    end if
  end if
  if (STKPNT <= 0) exit
end do

return

! End of ILASRT

end subroutine ILASRT
