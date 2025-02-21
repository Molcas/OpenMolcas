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
! Copyright (C) 1999, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CON_BLOCKS(IATP,IBTP,JATP,JBTP,IASM,IBSM,JASM,JBSM,ICONSPA,ICONSPB,NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)
! Does CI blocks IATP IBTP interact with blocks JATP JBTP
!
! Input
! =====
! IATP IBTP JATP JBTP : Supergroups, relative numbers
! IOCTPA, IOBTPB : Offset for type
! ICONSPA, ICONSPB : Connection matrices giving exciation
!                    level between two string types
! MXEXC : Largest excitation level
! IH_OCC_CONS : = 1 => Use only occupation conserving part of
!                      Hamiltonian
!
! Output
! ======
! INTERACT : =1 => The two blocks does interact
!
! Jeppe Olsen, April 99

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IATP, IBTP, JATP, JBTP, IASM, IBSM, JASM, JBSM, NOCTPA, ICONSPA(NOCTPA,NOCTPA), NOCTPB, &
                                 ICONSPB(NOCTPB,NOCTPB), MXEXC, IH_OCC_CONS
integer(kind=iwp), intent(inout) :: INTERACT
integer(kind=iwp) :: IA_EXC, IB_EXC

IA_EXC = ICONSPA(IATP,JATP)
IB_EXC = ICONSPB(IBTP,JBTP)
if (IH_OCC_CONS == 0) then
  ! Usual one- or two- electron operator
  if (MXEXC == 1) then
    if (((IA_EXC <= 1) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) .or. &
        ((IB_EXC <= 1) .and. (IATP == JATP) .and. (IASM == JASM))) INTERACT = 1
  else if (MXEXC == 2) then
    if (((IA_EXC <= 1) .and. (IB_EXC <= 1)) .or. ((IA_EXC == 2) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) .or. &
        ((IB_EXC == 2) .and. (IATP == JATP) .and. (IASM == JASM))) INTERACT = 1
  end if
!else
!  ! Orbital conserving part of  Hamiltonian
!  if ((IA_EXC == IB_EXC) .and. (IB_EXC <= 1)) then
!    IATP_ABS = IATP+IOCTPA-1
!    IBTP_ABS = IBTP+IOCTPB-1
!    JATP_ABS = JATP+IOCTPA-1
!    JBTP_ABS = JBTP+IOCTPB-1
!    ! Find Orb space where alpha strings differ
!    IPGAS = 0
!    IMGAS = 0
!    do IGAS=1,NGAS
!      IAEL = NELFSPGP(IGAS,IATP_ABS)
!      JAEL = NELFSPGP(IGAS,JATP_ABS)
!      if (IAEL-JAEL == 1) IPGAS = IGAS
!      if (IAEL-JAEL == -1)IMGAS = IGAS
!    end do
!    if (IPGAS /= 0) then
!      IPDIF = NELFSPGP(IPGAS,IBTP_ABS)-NELFSPGP(IPGAS,JBTP_ABS)
!    else
!      IPDIF = 0
!    end if
!    ! corresponding differences in beta
!    if (IMGAS /= 0) then
!      IMDIF = NELFSPGP(IMGAS,IBTP_ABS)-NELFSPGP(IMGAS,JBTP_ABS)
!    else
!      IMDIF = 0
!    end if
!    if ((IPGAS == 0) .and. (IMGAS == 0)) INTERACT = 1
!    if ((IPGAS /= 0) .and. (IMGAS /= 0)) then
!      if ((IPDIF == -1) .and. (IMDIF == 1)) INTERACT = 1
!    end if
!  end if
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from CONBLOCKS'
write(u6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
write(u6,*) ' IH_OCC_CONS, INTERACT = ',IH_OCC_CONS,INTERACT
#endif

end subroutine CON_BLOCKS
