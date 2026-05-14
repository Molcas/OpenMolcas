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

!#define _DEBUGPRINT_
subroutine DXTYP2_GAS(NDXTP,ITP,JTP,KTP,LTP,NOBTP,IL,IR,IPHGAS)
! Obtain types of I,J,K,l so
! <L!a+I a+K a L a J!R> is nonvanishing
! only combinations with type(I) >= type(K) and type(J) >= type(L)
! are included
!
! Intermediate occupations less than zero allowed for particle spaces
! (IPHGAS=2)

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: NDXTP
integer(kind=iwp), intent(_OUT_) :: ITP(*), JTP(*), KTP(*), LTP(*)
integer(kind=iwp), intent(in) :: NOBTP, IL(NOBTP), IR(NOBTP), IPHGAS(NOBTP)
integer(kind=iwp) :: IANNI1, IANNI2, ICREA1, ICREA2, IDIA, IJTP, IOBTP, KLTP, NANNI, NCREA, NDIF, NDIFT
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IDX

write(u6,*) ' DXTYP2_GAS in action'
write(u6,*) ' ===================='
write(u6,*) ' Occupation of left string'
call IWRTMA(IL,1,NOBTP,1,NOBTP)
write(u6,*) ' Occupation of right string'
call IWRTMA(IR,1,NOBTP,1,NOBTP)
#endif

! Number of differing occupations
NANNI = 0
NCREA = 0
NDIFT = 0

ICREA1 = 0
ICREA2 = 0
IANNI1 = 0
IANNI2 = 0
do IOBTP=1,NOBTP
  NDIFT = NDIFT+abs(IL(IOBTP)-IR(IOBTP))
  NDIF = IL(IOBTP)-IR(IOBTP)
  if (NDIF == 2) then
    ! two electrons of type IOBTP must be created
    ICREA1 = IOBTP
    ICREA2 = IOBTP
    NCREA = NCREA+2
  else if (NDIF == -2) then
    ! two electrons of type IOBTP must be annihilated
    IANNI1 = IOBTP
    IANNI2 = IOBTP
    NANNI = NANNI+2
  else if (NDIF == 1) then
    ! one electron of type IOBTP must be created
    if (NCREA == 0) then
      ICREA1 = IOBTP
    else
      ICREA2 = IOBTP
    end if
    NCREA = NCREA+1
  else if (NDIF == -1) then
    ! one electron of type IOBTP must be annihilated
    if (NANNI == 0) then
      IANNI1 = IOBTP
    else
      IANNI2 = IOBTP
    end if
    NANNI = NANNI+1
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' NCREA, NANNI ',NCREA,NANNI
write(u6,*) ' ICREA1, IANNI1 ',ICREA1,IANNI1
write(u6,*) ' ICREA2, IANNI2 ',ICREA2,IANNI2
#endif

NDXTP = 0
if (NDIFT <= 4) then
  if ((NCREA == 0) .and. (NANNI == 0)) then
    ! strings identical, include diagonal excitions  itp = jtp, ktp=ltp
    do IJTP=1,NOBTP
      if ((IR(IJTP) >= 1) .or. (IPHGAS(IJTP) == 2)) then
        do KLTP=1,IJTP
          if (((IJTP /= KLTP) .and. ((IR(KLTP) >= 1) .or. (IPHGAS(KLTP) == 2))) .or. &
              ((IJTP == KLTP) .and. ((IR(KLTP) >= 2) .or. (IPHGAS(KLTP) == 2)))) then
            NDXTP = NDXTP+1
            ITP(NDXTP) = IJTP
            JTP(NDXTP) = IJTP
            KTP(NDXTP) = KLTP
            LTP(NDXTP) = KLTP
          end if
        end do
      end if
    end do
  else if ((NCREA == 1) .and. (NANNI == 1)) then
    ! strings differ by single excitation
    ! diagonal excitation plus creation in ICREA1 and annihilation in IANNI1
    do IDIA=1,NOBTP
      if (((IDIA /= IANNI1) .and. ((IR(IDIA) >= 1) .or. (IPHGAS(IDIA) == 2))) .or. &
          ((IDIA == IANNI1) .and. ((IR(IDIA) >= 2) .or. (IPHGAS(IDIA) == 2)))) then
        NDXTP = NDXTP+1
        ITP(NDXTP) = max(ICREA1,IDIA)
        KTP(NDXTP) = min(ICREA1,IDIA)
        JTP(NDXTP) = max(IANNI1,IDIA)
        LTP(NDXTP) = min(IANNI1,IDIA)
      end if
    end do
  else if ((NCREA == 2) .and. (NANNI == 2)) then
    ! strings differ by double excitation
    NDXTP = 1
    ITP(1) = ICREA2
    KTP(1) = ICREA1
    JTP(1) = IANNI2
    LTP(1) = IANNI1
  end if
end if

#ifdef _DEBUGPRINT_
write(u6,'(A,I4)') ' Number of connecting double excitations ',NDXTP
if (NDXTP /= 0) then
  write(u6,*) '  ITYP KTYP LTYP JTYP'
  write(u6,*) '  ==================='
  do IDX=1,NDXTP
    write(u6,'(1X,4I5)') ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
  end do
end if
#endif

end subroutine DXTYP2_GAS
