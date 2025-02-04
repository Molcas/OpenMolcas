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

subroutine H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,IOCA,IOCB,NGAS,MXPNGAS,INTH0SPC,NELFTP)
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
! IOCPTSPC : Allowed occumulated occupation of each subspace
! NOCTPA :  Number of alpha occupation types
! NOCTPB : Number of beta occupation types
! IOCA : Occupation  of alpha string
! IOCB : Occupation  of beta string
!
! Jeppe Olsen, January 1996

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
dimension IOCPTSPC(2,MXPNGAS,*)
dimension IOCA(MXPNGAS,*), IOCB(MXPNGAS,*)
dimension NELFTP(*)
! Output
dimension INTH0SPC(NOCTPA,NOCTPB)

if (IH0SPC == 0) then
  ! All interactions allowed
  IONE = 1
  call ISETVC(INTH0SPC,IONE,NOCTPA*NOCTPB)
else
  ! Explicit construction of matrix giving partitionning of subspaces
  IZERO = 0
  call ISETVC(INTH0SPC,IZERO,NOCTPA*NOCTPB)

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
          if ((IEL < IOCPTSPC(1,IGAS,ISPC)) .or. (IEL > IOCPTSPC(2,IGAS,ISPC))) IAMOKAY = 0
        end do
        !write(u6,*) ' IAMOKAY = ',IAMOKAY
        ! Allowed
        if ((IAMOKAY == 1) .and. (INTH0SPC(IATP,IBTP) == 0)) INTH0SPC(IATP,IBTP) = ISPC
      end do
    end do
  end do
end if

NTEST = 0
if (NTEST >= 10) then
  write(u6,*)
  write(u6,*) ' ====================='
  write(u6,*) ' Output from  H0INTSPC'
  write(u6,*) ' ====================='
  write(u6,*)
  call IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end if

end subroutine H0INTSPC
