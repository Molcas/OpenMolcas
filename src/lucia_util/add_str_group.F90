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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************

subroutine ADD_STR_GROUP(NSTADD,IOFADD,ISTADD,NSTB,NSTA,ISTRING,IELOF,NELADD,NELTOT)
! Part of assembling strings in individual types to
! super group of strings
!
! Copying strings belonging to a given type to supergroup of strings
!
! Jeppe Olsen, for once improving performance of LUCIA
!
! Input
! =====
! NSTADD : Number of strings to be added
! IOFADD : First string to be added
! ISTADD : Strings to be added
! NSTB   : Number of strings belonging to lower gasspaces
! NSTA   : Number of strings belonging to higher gasspaces
! ISTRING: Supergroup of strings under construction
! IELOF  : Place of first electron to be added
! NELADD : Number of electrons to be added
! NELTOT : Total number of electrons

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NSTADD, IOFADD, ISTADD(*), NSTB, NSTA, IELOF, NELADD, NELTOT
integer(kind=iwp), intent(inout) :: ISTRING(*)
integer(kind=iwp) :: IADD2, IISTR, IOFF0, IOFF1, IOFF2, IOFFX, IOFFY, ISTA, ISTB

!write(u6,*) '  Inside ADD ...'
if (NSTA > 1) then
  do IISTR=1,NSTADD
    ! Address of A(1,IISTR,1)
    ! A(I(after),Igas,I(before))
    IOFFY = (IOFADD-2+IISTR)*NELADD
    IOFF1 = (IISTR-1)*NSTA+1
    IADD2 = NSTADD*NSTA
    IOFF2 = IOFF1-IADD2
    do ISTB=1,NSTB
      ! Address of A(1,IISTR,ISTB)
      !IOFF2 = IOFF1+(ISTB-1)*NSTADD*NSTA
      IOFF2 = IOFF2+IADD2
      IOFFX = IELOF-1+(IOFF2-2)*NELTOT
      do ISTA=1,NSTA
        IOFFX = IOFFX+NELTOT
        ISTRING(IOFFX+1:IOFFX+NELADD) = ISTADD(IOFFY+1:IOFFY+NELADD)
        !ISTRING(IELOF:IELOF+NELADD-1,IOFF2-1+ISTA) = ISTADD(1:NELADD,IOFADD-1+IISTR)
      end do
    end do
  end do
else if (NSTA == 1) then
  ! Address of A(1,IISTR,1)
  ! A(I(after),Igas,I(before))
  do ISTB=1,NSTB
    IOFF0 = (ISTB-1)*NSTADD
    IOFFY = (IOFADD-2)*NELADD
    IOFFX = IELOF-1+(IOFF0-1)*NELTOT
    do IISTR=1,NSTADD
      ! Address of A(1,IISTR,ISTB)
      !IOFF2 = IISTR+IOFF0
      !IOFFX = IELOF-1+(IOFF2-1)*NELTOT
      IOFFX = IOFFX+NELTOT
      IOFFY = IOFFY+NELADD
      ISTRING(IOFFX+1:IOFFX+NELADD) = ISTADD(IOFFY+1:IOFFY+NELADD)
      !ISTRING((IOFF2-1)*NELTOT+IEOLOF:(IOFF2-1)*NELTOT+IEOLOF+NELADD-1) = &
      !  ISTADD((IOFADD-1+IISTR-1)*NELADD+1:(IOFADD+IISTR-1)*NELADD)
      !ISTRING(IELOF:IELOF+NELADD-1,IOFF2) = ISTADD(1:NELADD,IOFADD-1+IISTR)
      !write(u6,*) ' New string from ADD'
      !call IWRTMA(ISTRING(IOFFX+1),1,NELADD,1,NELADD)
      !write(u6,*) ' IOFFX, IOFFY, NELADD',IOFFX,IOFFY,NELADD
    end do
  end do
end if

end subroutine ADD_STR_GROUP
