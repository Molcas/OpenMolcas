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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************

!#define_DEBUGPRINT_
subroutine ICISPC(MNRS10,MXRS30)
! Obtain internal CI spaces relevant for MRSDCI
!       /STRINP/+/LUCINP/ = > /CICISP/
! Jeppe Olsen, Dec 1990
!
! Input
! =====
! Information in STRINP
!
! Output
! ======
! Common block CICISP
!
! Internal CI spaces
!*******************************************************************
!   *  Basic space  * Allowed internal excit * Delta NA * Delta Nb *
!*******************************************************************
! 1 *  Zero order   *           0            *    0     *    0     *
!*******************************************************************
! ====================
!  Input Module
! ====================
!./Str_Info
!./MCLR_Data
! ====================
!  Output Module
! ====================
!module MCLR_Data
! NICISP : Number of internal CI spaces constructed
! IASTFI : Alpha string type for internal CI space
! IBSTFI : Beta string type for internal CI space
! IACTI  : Given internal space is active
! MXR3IC : Max number of elecs in RAS 3 space for internal CI space
! MNR1IC : Min number of elecs in RAS 1 space for internal CI space
! IZCI   : Internal zero order space
! NELCI : Number of electrons per CI space
! obtained by (IEX-1) fold internal excitation, with a NAEL + DELTAA
! alpha electrons and  NBEL + DELTAB beta electrons

use Str_Info, only: IAZTP, IBZTP, NELEC
use MCLR_Data, only: IACTI, IASTFI, IBSTFI, MNR1IC, MNR3IC, MXR1IC, MXR3IC, NAELCI, NBELCI, NELCI, NICISP, NORB1, NORB2
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MNRS10, MXRS30
integer(kind=iwp) :: ICI

ICI = 1
MNR1IC(ICI) = MNRS10
MXR3IC(ICI) = MXRS30
IASTFI(ICI) = IAZTP
IBSTFI(ICI) = IBZTP
NAELCI(ICI) = NELEC(IAZTP)
NBELCI(ICI) = NELEC(IBZTP)
NELCI(ICI) = NAELCI(ICI)+NBELCI(ICI)
IACTI(1) = 1
NICISP = ICI
! Number and distribution of electrons in each space
!do IEX=1,3
!  do IDA=-4,2
!    do IDB=-4,2
!      if (IRCI(IEX,IDA+5,IDB+5) /= 0) then
!        ICI = IRCI(IEX,IDA+5,IDB+5)
!        NAELCI(ICI) = NELEC(IASTFI(ICI))
!        NBELCI(ICI) = NELEC(IBSTFI(ICI))
!        NELCI(ICI) = NAELCI(ICI)+NBELCI(ICI)
!      end if
!    end do
!  end do
!end do

! Default max in RAS1 and min in RAS3
MXR1IC(1:NICISP) = min(2*NORB1,NELCI(1:NICISP))
MNR3IC(1:NICISP) = max(0,NELCI(1:NICISP)-2*(NORB1+NORB2))

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of internal CI spaces ',NICISP
write(u6,*) ' Space a-type b-type nael nbel mnrs1 mxrs1 mnrs3 mxrs3'
write(u6,*) ' ====================================================='
do ICI=1,NICISP
  if (IACTI(ICI) == 1) &
    write(u6,'(I5,2I7,2I5,4I6)') ICI,IASTFI(ICI),IBSTFI(ICI),NAELCI(ICI),NBELCI(ICI),MNR1IC(ICI),MXR1IC(ICI),MNR3IC(ICI),MXR3IC(ICI)
end do
#endif

return

end subroutine ICISPC
