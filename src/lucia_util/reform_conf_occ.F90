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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

subroutine REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
! Reform between two ways of writing occupations
!
! IOCC_EXP : Occupation in expanded form, i.e. the orbital for each
!            electron is given
!
! IOCC_PCK  : Occupation is given in packed form, i.e. each occupied
!             orbitals is given once, and a negative index indicates
!             a double occupattion
!
! IWAY = 1 Expanded to Packed form
! IWAY = 2 Packed to expanded form
!
! Jeppe Olsen, Nov. 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NEL, IOCC_EXP(NEL), NOCOB, IOCC_PCK(NOCOB), IWAY
integer(kind=iwp) :: IEL, IOCC, IORB, JORB, NTEST

if (IWAY == 1) then

  ! Expanded => Packed form

  ! Loop over electrons
  IEL = 1
  IOCC = 0
  do
    !IEL = IEL+1
    if (IEL < NEL) then
      if (IOCC_EXP(IEL) == IOCC_EXP(IEL+1)) then
        IOCC = IOCC+1
        IOCC_PCK(IOCC) = -IOCC_EXP(IEL)
        IEL = IEL+2
      else
        IOCC = IOCC+1
        IOCC_PCK(IOCC) = IOCC_EXP(IEL)
        IEL = IEL+1
      end if
    else
      ! Last occupation was not identical to previous, so single occupied
      IOCC = IOCC+1
      IOCC_PCK(IOCC) = IOCC_EXP(IEL)
      IEL = IEL+1
    end if
    if (IEL > NEL) exit
  end do

else if (IWAY == 2) then

  ! Packed to expanded form

  IEL = 0
  do IORB=1,NOCOB
    if (IOCC_PCK(IORB) < 0) then
      JORB = -IOCC_PCK(IORB)
      IEL = IEL+1
      IOCC_EXP(IEL) = JORB
      IEL = IEL+1
      IOCC_EXP(IEL) = JORB
    end if
  end do
else
  write(u6,*) ' REFORM_CONF... in error, IWAY = ',IWAY
  !stop ' REFORM_CONF... in error, IWAY'
  call SYSABENDMSG('lucia_util/reform_conv','Internal error','')
end if
! End of IWAY switch

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Reforming form of configuration'
  if (IWAY == 1) then
    write(u6,*) ' Expanded to packed form'
  else
    write(u6,*) ' Packed to expanded form'
  end if
  write(u6,*) ' IOCC_EXP :'
  call IWRTMA(IOCC_EXP,1,NEL,1,NEL)
  write(u6,*) ' IOCC_PCK :'
  call IWRTMA(IOCC_PCK,1,NOCOB,1,NOCOB)
end if

end subroutine REFORM_CONF_OCC
