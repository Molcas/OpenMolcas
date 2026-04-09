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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

function NCASTR_MCLR(IAC,NSTTPI,NTPSTI,ICLSI,NOBATP,NOBTP,IELPTP)
! Number of allowed annihilation/creations from a given group of strings
!
! Jeppe Olsen, June 1994

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: NCASTR_MCLR
integer(kind=iwp), intent(in) :: IAC, NTPSTI, NSTTPI(NTPSTI), ICLSI, NOBTP, NOBATP(NOBTP), IELPTP(NOBTP,*)
integer(kind=iwp) :: ICLSO, IOBTP, ISTTP, ITPO, LCA, NENTRY

LCA = 0
!write(u6,*) ' NCASTR: IELPTP array'
!write(u6,*) ' ===================='
!write(u6,*) ' Number of occupation types ',NTPSTI
!call iwrtma(IELPTP,NOBTP,NTPSTI,NOBTP,NTPSTI)
!write(u6,*) ' NOBTP NTPSTI ',NOBTP,NTPSTI
!write(u6,*) ' Number of strings per type:'
!write(u6,*) (NSTTPI(I),I=1,NTPSTI)
do IOBTP=1,NOBTP
  do ISTTP=1,NTPSTI
    ! Type of resulting string
    call NEWTYP_MCLR(ICLSI,ISTTP,[IAC],[IOBTP],1,ICLSO,ITPO)
    !write(u6,*) ' IOBTP ISTTP => ITPO, ICLSO'
    !write(u6,*) IOBTP,ISTTP,ITPO,ICLSO
    !write(u6,*) ' IELPTP = ',IELPTP(IOBTP,ISTTP)
    !write(u6,*) ' NOBATP = ',NOBATP(IOBTP)
    if (IAC == 1) then
      NENTRY = IELPTP(IOBTP,ISTTP)
    else
      NENTRY = NOBATP(IOBTP)-IELPTP(IOBTP,ISTTP)
    end if
    !write(u6,*) ' NENTRY = ',NENTRY
    if (ITPO > 0) LCA = LCA+NENTRY*NSTTPI(ISTTP)

  end do
end do

!write(u6,*) ' Number of generated strings ',LCA
NCASTR_MCLR = LCA

return

end function NCASTR_MCLR
