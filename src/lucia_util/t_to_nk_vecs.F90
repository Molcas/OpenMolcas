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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

subroutine T_TO_NK_VECS(T,KORB,C,LUCIN,LUCOUT,NSSOA,NSSOB,NBLOCK,IBLOCK,NAEL,NBEL,IASTR,IBSTR,IBLTP,NSMST,ICISTR,NORB,IKAOCC,IKBOCC)
! Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
!
! Both files are initially rewound
!
! Jeppe Olsen, Feb. 1998

use lucia_data, only: IDISK
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: T
integer(kind=iwp), intent(in) :: KORB, LUCIN, LUCOUT, NSMST, NSSOA(NSMST,*), NSSOB(NSMST,*), NBLOCK, IBLOCK(8,NBLOCK), NAEL, NBEL, &
                                 IBLTP(*), ICISTR, NORB
real(kind=wp), intent(_OUT_) :: C(*)
integer(kind=iwp), intent(inout) :: IASTR(NAEL,*), IBSTR(NBEL,*)
integer(kind=iwp), intent(_OUT_) :: IKAOCC(*), IKBOCC(*)
integer(kind=iwp) :: IA, IAMPACK, IASM, IATP, IB, IBSM, IBTP, IDET, IDUM(1), IMZERO, IRESTR, JAEL, JBEL, JBLOCK, JSTR, KABOCC, &
                     KOCC, LDET, MINIA, NASTR1, NBSTR1, NIA, NIB
real(kind=wp) :: T2

IDISK(LUCIN) = 0
IDISK(LUCOUT) = 0

T2 = T**2
do JBLOCK=1,NBLOCK
  IATP = IBLOCK(1,JBLOCK)
  IBTP = IBLOCK(2,JBLOCK)
  IASM = IBLOCK(3,JBLOCK)
  IBSM = IBLOCK(4,JBLOCK)
  !write(u6,*) ' IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
  ! Obtain alpha strings of sym IASM and type IATP
  IDUM(1) = 0
  call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM,IDUM)
  ! Occupation of orb KORB
  do JSTR=1,NASTR1
    KOCC = 0
    do JAEL=1,NAEL
      if (IASTR(JAEL,JSTR) == KORB) KOCC = 1
    end do
    IKAOCC(JSTR) = KOCC
  end do
  !write(u6,*) ' IKAOCC array'
  !call IWRTMA(IKAOCC,1,NASTR1,1,NASTR1)

  ! Obtain Beta  strings of sym IBSM and type IBTP
  IDUM(1) = 0
  call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM,IDUM)
  !write(u6,*) ' After GETSTR, NBSTR1=',NBSTR1
  ! Occupation of orb KORB
  do JSTR=1,NBSTR1
    !write(u6,*) ' JSTR = ',JSTR
    KOCC = 0
    do JBEL=1,NBEL
      !write(u6,*) JBEL,IBSTR(JBEL,JSTR)
      if (IBSTR(JBEL,JSTR) == KORB) KOCC = 1
    end do
    IKBOCC(JSTR) = KOCC
  end do
  !write(u6,*) ' IKBOCC array'
  !call IWRTMA(IKBOCC,1,NBSTR1,1,NBSTR1)

  if (IBLTP(IASM) == 2) then
    IRESTR = 1
  else
    IRESTR = 0
  end if
  !write(u6,*) ' IBLTP ',IBLTP(IASM)

  NIA = NSSOA(IASM,IATP)
  NIB = NSSOB(IBSM,IBTP)
  !write(u6,*) ' NIA NIB ',NIA,NIB

  IMZERO = 0
  if (ICISTR >= 2) then
    ! Read in a Type-Type-symmetry block
    call IDAFILE(LUCIN,2,IDUM,1,IDISK(LUCIN))
    LDET = IDUM(1)
    call IDAFILE(LUCIN,2,IDUM,1,IDISK(LUCIN))
    call FRMDSC(C,LDET,-1,LUCIN,IMZERO,IAMPACK)
  end if
  if (IMZERO /= 1) then

    IDET = 0
    do IB=1,NIB
      if ((IRESTR == 1) .and. (IATP == IBTP)) then
        MINIA = IB
      else
        MINIA = 1
      end if
      do IA=MINIA,NIA

        IDET = IDET+1
        !write(u6,*) ' IA IB IDET',IA,IB,IDET
        KABOCC = IKAOCC(IA)+IKBOCC(IB)
        if (KABOCC == 1) then
          C(IDET) = T*C(IDET)
        else if (KABOCC == 2) then
          C(IDET) = T2*C(IDET)
        end if
      end do
      ! End of loop over alpha strings
    end do
    ! End of loop over beta strings

  end if
  ! End of if statement for nonvanishing blocks
  ! Save result on LUCOUT
  IDUM(1) = LDET
  call ITODS(IDUM,1,-1,LUCOUT)
  call TODSC(C,LDET,-1,LUCOUT)
end do
! End of loop over blocks
! This is the end, the end of every file my friend, the end
IDUM(1) = -1
call ITODS(IDUM,1,-1,LUCOUT)

end subroutine T_TO_NK_VECS
