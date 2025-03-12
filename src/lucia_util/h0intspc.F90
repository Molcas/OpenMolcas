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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine H0INTSPC(IH0SPC,NPTSPC,NOCTPA,NOCTPB,IOCA,IOCB,NGAS,MXPNGAS,INTH0SPC,NELFTP)
! Set up INTH0SPC : Division of CI space, so only determinants
!                   belonging to the same space  have nonvanishing
!                   matrix elements of H0
!
! =====
! Input
! =====
!
! IH0SPC : /= 0 : Interacting subspaces have been defined
!          == 0 : Interacting subspaces not defined, let
!                 everything interact
! NPTSPC : Number of subspaces defined
! NOCTPA :  Number of alpha occupation types
! NOCTPB : Number of beta occupation types
! IOCA : Occupation  of alpha string
! IOCB : Occupation  of beta string
!
! Jeppe Olsen, January 1996

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IH0SPC, NPTSPC, MXPNGAS, NOCTPA, NOCTPB, IOCA(MXPNGAS,NOCTPA), IOCB(MXPNGAS,NOCTPB), NGAS, &
                                 NELFTP(*)
integer(kind=iwp), intent(out) :: INTH0SPC(NOCTPA,NOCTPB)
integer(kind=iwp) :: IAMOKAY, IATP, IBTP, IEL, IGAS, ISPC

if (IH0SPC == 0) then
  ! All interactions allowed
  INTH0SPC(:,:) = 1
else
  ! Explicit construction of matrix giving partitionning of subspaces
  INTH0SPC(:,:) = 0

  do ISPC=1,NPTSPC
    do IATP=1,NOCTPA
      do IBTP=1,NOCTPB
        IAMOKAY = 1
        IEL = 0
        !write(u6,*) ' ISPC IATP IBTP ',ISPC,IATP,IBTP
        do IGAS=1,NGAS
          IEL = IEL+NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
          !write(u6,*) ' IGAS IEL ',IGAS,IEL
          !write(u6,*) ' Limits :',IOCPTSPC(1,IGAS,ISPC),IOCPTSPC(2,IGAS,ISPC)
          !if ((IEL < IOCPTSPC(1,IGAS,ISPC)) .or. (IEL > IOCPTSPC(2,IGAS,ISPC))) IAMOKAY = 0
          if (IEL /= 0) IAMOKAY = 0
        end do
        !write(u6,*) ' IAMOKAY = ',IAMOKAY
        ! Allowed
        if ((IAMOKAY == 1) .and. (INTH0SPC(IATP,IBTP) == 0)) INTH0SPC(IATP,IBTP) = ISPC
      end do
    end do
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ====================='
write(u6,*) ' Output from  H0INTSPC'
write(u6,*) ' ====================='
write(u6,*)
call IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
#endif

end subroutine H0INTSPC
