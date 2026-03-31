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

subroutine PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,COEFF,SGM,PSI)
! Purpose: Add to the wave function SGM the result of applying
! an operator to PSI. The operator is a single creator (IMODE=1)
! or annihilator (IMODE=-1) multiplied by COEFF. Only the
! sector of SGM described by IFSBTAB1 will be updated.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IMODE, ISORB, IORBTAB(*), ISSTAB(*), IFSBTAB1(*), IFSBTAB2(*)
real(kind=wp) :: COEFF, SGM(*), PSI(*)
integer(kind=iwp) :: I, IBLKPOS1, IBLKPOS2, IFSB1, IFSB2, IPOS1, IPOS2, ISBS1, ISBS2, ISP, ISPART, ISST, ISST1, ISST2, &
                     ISSTARR(50), ISUM, J, KHSH2, KOINFO, KPOS, KSBS1, KSBS2, KSBSAN, KSBSCR, KSBSOP, KSORB, KSSTAN, KSSTCR, &
                     KSSTOP, KSSTTB, KSTARR1, KSTARR2, MORSBITS, NASPRT, NDI, NDJ, NFSB1, NHSH2, NPOP1, NSBS, NSBS1, NSBS2, NSSTP
real(kind=wp) :: CFFPHS, SCL
integer(kind=iwp), allocatable :: SBSET(:)

if (COEFF == Zero) return
! The orbital table:
NASPRT = IORBTAB(9)
KOINFO = 19
ISPART = IORBTAB(KOINFO+6+8*(ISORB-1))
KSORB = IORBTAB(KOINFO+7+8*(ISORB-1))
! The substring table:
MORSBITS = ISSTAB(6)
NSSTP = ISSTAB(7)
KSSTTB = 15
KSSTAN = ISSTAB(9)
KSSTCR = ISSTAB(10)
KSBSAN = ISSTAB(13)
KSBSCR = ISSTAB(14)
if (IMODE == 1) then
  KSSTOP = KSSTAN
  KSBSOP = KSBSAN
else
  KSSTOP = KSSTCR
  KSBSOP = KSBSCR
end if
! The FS blocks of the SGM wave function:
NFSB1 = IFSBTAB1(3)
KSTARR1 = 8
! The FS blocks of the PSI wave function:
NHSH2 = IFSBTAB2(6)
KHSH2 = IFSBTAB2(7)
KSTARR2 = 8
! Make an array with nr of earlier substrings for each
! substring type:
call mma_allocate(SBSET,NSSTP,Label='SBSET')
ISUM = 0
do ISST=1,NSSTP
  SBSET(ISST) = ISUM
  NSBS = ISSTAB(KSSTTB+5*(ISST-1))
  ISUM = ISUM+NSBS
end do

! Loop over FS blocks of the SGM wave function
do IFSB1=1,NFSB1
  KPOS = KSTARR1+(NASPRT+2)*(IFSB1-1)
  do ISP=1,NASPRT
    ISSTARR(ISP) = IFSBTAB1(KPOS-1+ISP)
  end do
  IBLKPOS1 = IFSBTAB1(KPOS+NASPRT+1)
  ! Initial values for lower and higher dimensions.
  ! Also, extra phase factor due to spin orbitals in higher substrings.
  NDI = 1
  do ISP=1,ISPART-1
    ISST1 = ISSTARR(ISP)
    NSBS1 = ISSTAB(KSSTTB+5*(ISST1-1))
    NDI = NDI*NSBS1
  end do
  NDJ = 1
  CFFPHS = COEFF
  do ISP=ISPART+1,NASPRT
    ISST1 = ISSTARR(ISP)
    NSBS1 = ISSTAB(KSSTTB+0+5*(ISST1-1))
    NPOP1 = ISSTAB(KSSTTB+1+5*(ISST1-1))
    NDJ = NDJ*NSBS1
    if (NPOP1 /= 2*(NPOP1/2)) CFFPHS = -CFFPHS
  end do
  ISST1 = ISSTARR(ISPART)
  NSBS1 = ISSTAB(KSSTTB+5*(ISST1-1))

  ! Modify the bra substring type by annih or creating ISORB
  ISST2 = ISSTAB(KSSTOP-1+KSORB+MORSBITS*(ISST1-1))
  if (ISST2 == 0) goto 200

  ! Determine dimensions for multiple daxpy:
  ! Dimension for earlier subpartitions is NDI
  ! Dimension for later   subpartitions is NDJ
  ! Dimensions for present subpartition are NSBS1,NSBS2

  NSBS2 = ISSTAB(KSSTTB+5*(ISST2-1))
  ISSTARR(ISPART) = ISST2
  ! Get the corresponding FS block number
  call HSHGET(ISSTARR,NASPRT,NASPRT+2,IFSBTAB2(KSTARR2),NHSH2,IFSBTAB2(KHSH2),IFSB2)
  ISSTARR(ISPART) = ISST1
  if (IFSB2 == 0) goto 200
  KPOS = KSTARR2+(NASPRT+2)*(IFSB2-1)
  IBLKPOS2 = IFSBTAB2(KPOS+NASPRT+1)
  ! Now loop over ket substrings in this subpartition
  do KSBS1=1,NSBS1
    ISBS1 = KSBS1+SBSET(ISST1)
    ISBS2 = ISSTAB(KSBSOP-1+KSORB+MORSBITS*(ISBS1-1))
    if (ISBS2 == 0) goto 100
    if (ISBS2 > 0) then
      SCL = CFFPHS
      ISBS2 = ISBS2
    else
      SCL = -CFFPHS
      ISBS2 = -ISBS2
    end if
    KSBS2 = ISBS2-SBSET(ISST2)

    ! CALL some multiple daxpy...
    if (NDI == 1) then
      if (NDJ == 1) then
        IPOS1 = IBLKPOS1+(KSBS1-1)
        IPOS2 = IBLKPOS2+(KSBS2-1)
        SGM(IPOS1) = SGM(IPOS1)+SCL*PSI(IPOS2)
      else
        do J=0,NDJ-1
          IPOS1 = IBLKPOS1+(KSBS1-1+NSBS1*J)
          IPOS2 = IBLKPOS2+(KSBS2-1+NSBS2*J)
          SGM(IPOS1) = SGM(IPOS1)+SCL*PSI(IPOS2)
        end do
      end if
    else
      if (NDJ == 1) then
        do I=0,NDI-1
          IPOS1 = IBLKPOS1+I+NDI*(KSBS1-1)
          IPOS2 = IBLKPOS2+I+NDI*(KSBS2-1)
          SGM(IPOS1) = SGM(IPOS1)+SCL*PSI(IPOS2)
        end do
      else
        do J=0,NDJ-1
          do I=0,NDI-1
            IPOS1 = IBLKPOS1+I+NDI*(KSBS1-1+NSBS1*J)
            IPOS2 = IBLKPOS2+I+NDI*(KSBS2-1+NSBS2*J)
            SGM(IPOS1) = SGM(IPOS1)+SCL*PSI(IPOS2)
          end do
        end do
      end if
    end if

100 continue
  end do

! End of loop over FS blocks
200 continue
end do
call mma_deallocate(SBSET)

end subroutine PRIMSGM
