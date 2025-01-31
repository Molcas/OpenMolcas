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

subroutine DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NOBTP,IL,IR)
! Obtain types of I,J,K,l so
! <L!a+I a+K a L a J!R> is nonvanishing
! only combinations with type(I) >= type(K) and type(J) >= type(L)
! are included

integer IL(NOBTP), IR(NOBTP)
integer ITP(*), JTP(*), KTP(*), LTP(*)

NTEST = 0
if (NTEST >= 100) then
  write(6,*) ' DXTYP_GAS in action'
  write(6,*) ' ==================='
  write(6,*) ' Occupation of left string'
  call IWRTMA(IL,1,NOBTP,1,NOBTP)
  write(6,*) ' Occupation of right string'
  call IWRTMA(IR,1,NOBTP,1,NOBTP)
end if

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

if (NTEST >= 1000) then
  write(6,*) ' NCREA, NANNI ',NCREA,NANNI
  write(6,*) ' ICREA1, IANNI1 ',ICREA1,IANNI1
  write(6,*) ' ICREA2, IANNI2 ',ICREA2,IANNI2
end if

NDXTP = 0
if (NDIFT > 4) then
  NDXTP = 0
else
  if ((NCREA == 0) .and. (NANNI == 0)) then
    ! strings identical, include diagonal excitions  itp = jtp, ktp=ltp
    do IJTP=1,NOBTP
      if (IR(IJTP) >= 1) then
        do KLTP=1,IJTP
          if (((IJTP /= KLTP) .and. (IR(KLTP) >= 1)) .or. ((IJTP == KLTP) .and. (IR(KLTP) >= 2))) then
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
      if (((IDIA /= IANNI1) .and. (IR(IDIA) >= 1)) .or. ((IDIA == IANNI1) .and. (IR(IDIA) >= 2))) then
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

if (NTEST /= 0) then
  write(6,'(A,I4)') ' Number of connecting double excitations ',NDXTP
  if (NDXTP /= 0) then
    write(6,*) '  ITYP KTYP LTYP JTYP'
    write(6,*) '  ==================='
    do IDX=1,NDXTP
      write(6,'(1X,4I5)') ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
    end do
  end if
end if

end subroutine DXTYP_GAS
