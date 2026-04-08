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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************
! Some GAS routines
!
! Nomenclature
!
! Group of strings : set of strings with a given number of orbitals
!                   in a given GASspace
!
! Supergroup of strings : product of NGAS groups, i.e. a string with a
!                         given numb er of electrons in each GAS space
!
! Type of string : Type is defined by the total number of electrons
!                  in the string. A type will therefore in general
!                  consists of several supergroups

!#define _DEBUGPRINT_
subroutine SPGP_AC(INSPGRP,NINSPGRP,IOUTSPGRP,NOUTSPGRP,NGAS,MXPNGAS,IAC,ISPGRP_AC,IBASEIN,IBASEOUT)
! Annihilation/Creation mapping of strings
!
! Jeppe Olsen, April 1, 1997

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MXPNGAS, NINSPGRP, IOUTSPGRP(MXPNGAS,*), NOUTSPGRP, NGAS, IAC, IBASEIN, IBASEOUT
integer(kind=iwp), intent(inout) :: INSPGRP(MXPNGAS,*), ISPGRP_AC(NGAS,*)
integer(kind=iwp) :: IAMOKAY, IGAS, ISPGRP, ITO, JGAS, JSPGRP, NELIN, NELOUT

! Check first that supergroups + IAC information is consistent
NELIN = sum(INSPGRP(1:NGAS,IBASEIN))
NELOUT = sum(IOUTSPGRP(1:NGAS,IBASEOUT))
if (.not. (((IAC == 1) .and. (NELIN == NELOUT+1)) .or. ((IAC == 2) .and. (NELIN == NELOUT-1)))) then
  write(u6,*) ' Inconsistent data provided to SPGP_AC'
  write(u6,*) ' NELIN NELOUT IAC=',NELIN,NELOUT,IAC
  !stop ' Inconsistent data provided to SPGRP_AC'
  call SYSABENDMSG('lucia_util/spgp_ac','Internal error','')
end if

do ISPGRP=IBASEIN,IBASEIN+NINSPGRP-1
  do IGAS=1,NGAS
    if (IAC == 1) then
      INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP)-1
    else if (IAC == 2) then
      INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP)+1
    end if
    ! Find corresponding supergroup
    ITO = 0
    do JSPGRP=IBASEOUT,IBASEOUT+NOUTSPGRP-1
      IAMOKAY = 1
      do JGAS=1,NGAS
        if (INSPGRP(JGAS,ISPGRP) /= IOUTSPGRP(JGAS,JSPGRP)) IAMOKAY = 0
      end do
      if (IAMOKAY == 1) ITO = JSPGRP
    end do
    ISPGRP_AC(IGAS,ISPGRP) = ITO
    ! And clean up
    if (IAC == 1) then
      INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP)+1
    else if (IAC == 2) then
      INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP)-1
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Input supergroups'
call IWRTMA(INSPGRP(1,IBASEIN),NGAS,NINSPGRP,MXPNGAS,NINSPGRP)
write(u6,*) ' Output supergroups'
call IWRTMA(IOUTSPGRP(1,IBASEOUT),NGAS,NOUTSPGRP,MXPNGAS,NOUTSPGRP)

write(u6,*) ' Output from SPGP_AC'
write(u6,*) ' ==================='
write(u6,*)
write(u6,*) ' IAC = ',IAC
write(u6,*) ' Mapping :'
call IWRTMA(ISPGRP_AC(1,IBASEIN),NGAS,NINSPGRP,NGAS,NINSPGRP)
#endif

end subroutine SPGP_AC
