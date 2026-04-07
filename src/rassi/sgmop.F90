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

subroutine SGMOP(IMODE,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,COEFF,SGM,PSI)
! Purpose: Add to the wave function SGM the result of applying
! an operator to PSI. The operator is a sum of creators (IMODE=1)
! or annihilators (IMODE=-1) multiplied by coefficients. Only the
! sector of SGM described by IFSBTAB1 will be updated.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IMODE, IORBTAB(*), ISSTAB(*), IFSBTAB1(*), IFSBTAB2(*)
real(kind=wp) :: COEFF(*), SGM(*), PSI(*)
integer(kind=iwp) :: I, IBLKPOS1, IBLKPOS2, IERR, IFSB1, IFSB2, IPH, IPHARR(50), IPOS1, IPOS2, ISBS1, ISBS2, ISORB, ISPART, ISST, &
                     ISST1, ISST2, ISSTARR(50), ISUM, J, KHSH2, KOINFO, KPOS, KSBS1, KSBS2, KSBSAN, KSBSCR, KSBSOP, KSORB, KSSTAN, &
                     KSSTCR, KSSTOP, KSSTTB, KSTARR1, KSTARR2, MORSBITS, NASORB, NASPRT, NDETS1, NDETS2, NDI, NDIARR(50), NDJ, &
                     NDJARR(50), NFSB1, NHSH2, NPOP1, NSBS1, NSBS2, NSSTP
real(kind=wp) :: CFFPHS, SCL
integer(kind=iwp), allocatable :: SBSET(:)

! The orbital table:
NASPRT = IORBTAB(9)
NASORB = IORBTAB(4)
KOINFO = 19
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
NDETS1 = IFSBTAB1(5)
KSTARR1 = 8
! The FS blocks of the PSI wave function:
NDETS2 = IFSBTAB2(5)
NHSH2 = IFSBTAB2(6)
KHSH2 = IFSBTAB2(7)
KSTARR2 = 8
! Make an array with nr of earlier substrings for each
! substring type:
call mma_allocate(SBSET,NSSTP,Label='SBSET')
ISUM = 0
do ISST=1,NSSTP
  SBSET(ISST) = ISUM
  ISUM = ISUM+ISSTAB(KSSTTB+5*(ISST-1))
end do

! Loop over FS blocks of the SGM wave function
do IFSB1=1,NFSB1
  KPOS = KSTARR1+(NASPRT+2)*(IFSB1-1)
  ISSTARR(1:NASPRT) = IFSBTAB1(KPOS:KPOS+NASPRT-1)
  IBLKPOS1 = IFSBTAB1(KPOS+NASPRT+1)
  !TEST write(u6,'(1x,a,8I8)') 'SGM FSB1,IBLKPOS1=',IBLKPOS1
  !TEST write(u6,'(1x,a,8I8)') 'SGM FSB1=',(ISSTARR(I),I=1,NASPRT)
  ! Initial values for lower and higher dimensions.
  ! Also, extra phase factor due to spin orbitals in higher substrings.
  NDI = 1
  do ISPART=1,NASPRT
    NDIARR(ISPART) = NDI
    ISST1 = ISSTARR(ISPART)
    NDI = NDI*ISSTAB(KSSTTB+5*(ISST1-1))
  end do
  NDJ = 1
  IPH = 1
  do ISPART=NASPRT,1,-1
    NDJARR(ISPART) = NDJ
    IPHARR(ISPART) = IPH
    ISST1 = ISSTARR(ISPART)
    NSBS1 = ISSTAB(KSSTTB+0+5*(ISST1-1))
    NPOP1 = ISSTAB(KSSTTB+1+5*(ISST1-1))
    NDJ = NDJ*NSBS1
    if (NPOP1 /= 2*(NPOP1/2)) IPH = -IPH
  end do
  ! Loop over active orbitals:
  do ISORB=1,NASORB
    if (COEFF(ISORB) == Zero) cycle
    ISPART = IORBTAB(KOINFO+6+8*(ISORB-1))
    CFFPHS = real(IPHARR(ISPART),kind=wp)*COEFF(ISORB)
    KSORB = IORBTAB(KOINFO+7+8*(ISORB-1))
    !TEST write(u6,'(1x,a,8I8)') 'ISORB,ISPART,KSORB:',ISORB,ISPART,KSORB
    ISST1 = ISSTARR(ISPART)
    NSBS1 = ISSTAB(KSSTTB+5*(ISST1-1))

    ! Modify the bra substring type by annih or creating ISORB
    ISST2 = ISSTAB(KSSTOP-1+KSORB+MORSBITS*(ISST1-1))
    if (ISST2 == 0) cycle

    !TEST write(u6,'(1x,a,8I8)') 'ISST1,ISST2:',ISST1,ISST2
    ! Determine dimensions for multiple daxpy:
    ! Dimension for earlier subpartitions is NDI
    ! Dimension for later   subpartitions is NDJ
    ! Dimensions for present subpartition are NSBS1,NSBS2
    NDI = NDIARR(ISPART)
    NDJ = NDJARR(ISPART)

    NSBS2 = ISSTAB(KSSTTB+5*(ISST2-1))
    ISSTARR(ISPART) = ISST2
    ! Get the corresponding FS block number
    !TEST write(u6,'(1x,a,8I8)') 'PSI FSB2=',(ISSTARR(I),I=1,NASPRT)
    call HSHGET(ISSTARR,NASPRT,NASPRT+2,IFSBTAB2(KSTARR2),NHSH2,IFSBTAB2(KHSH2),IFSB2)
    !TEST write(u6,'(1x,a,8I8)') 'IFSB1,IFSB2:',IFSB1,IFSB2
    ISSTARR(ISPART) = ISST1
    if (IFSB2 == 0) cycle
    KPOS = KSTARR2+(NASPRT+2)*(IFSB2-1)
    IBLKPOS2 = IFSBTAB2(KPOS+NASPRT+1)
    ! Now loop over ket substrings in this subpartition
    do KSBS1=1,NSBS1
      ISBS1 = KSBS1+SBSET(ISST1)
      ISBS2 = ISSTAB(KSBSOP-1+KSORB+MORSBITS*(ISBS1-1))
      if (ISBS2 == 0) then
        cycle
      else if (ISBS2 > 0) then
        SCL = CFFPHS
        ISBS2 = ISBS2
      else
        SCL = -CFFPHS
        ISBS2 = -ISBS2
      end if
      KSBS2 = ISBS2-SBSET(ISST2)

      ! CALL some multiple daxpy...
      do I=0,NDI-1
        do J=0,NDJ-1
          IPOS1 = IBLKPOS1+I+NDI*(KSBS1-1+NSBS1*J)
          IPOS2 = IBLKPOS2+I+NDI*(KSBS2-1+NSBS2*J)
          SGM(IPOS1) = SGM(IPOS1)+SCL*PSI(IPOS2)
          IERR = 0
          if ((IPOS1 < 1) .or. (IPOS1 > NDETS1)) IERR = 1
          if ((IPOS2 < 1) .or. (IPOS2 > NDETS2)) IERR = 1
          if (IERR /= 0) then
            write(u6,*) ' SGMOP addressing error.'
            write(u6,'(1x,a,8I8)') ' SGM dimension NDETS1=',NDETS1
            write(u6,'(1x,a,8I8)') ' Position IPOS1=',IPOS1
            write(u6,'(1x,a,8I8)') ' PSI dimension NDETS2=',NDETS2
            write(u6,'(1x,a,8I8)') ' Position IPOS2=',IPOS2
            call ABEND()
          end if
          ! Temporary test print:
          !TEST if (PSI(IPOS2) /= Zero) then
          !TEST   write(u6,'(1x,f16.8,2i8)') SCL,IPOS1,IPOS2
          !TEST   write(u6,'(1x,a,8I8)') 'IFSB1,IFSB2:',IFSB1,IFSB2
          !TEST   write(u6,'(1x,a,8I8)') 'IBLKPOS1,NDI,NSBS1:',IBLKPOS1,NDI,NSBS1
          !TEST   write(u6,'(1x,a,8I8)') 'IBLKPOS2,NDI,NSBS2:',IBLKPOS2,NDI,NSBS2
          !TEST   write(u6,'(1x,a,8I8)') 'I,KSBS1,J:',I,KSBS1,J
          !TEST   write(u6,'(1x,a,8I8)') 'I,KSBS2,J:',I,KSBS2,J
          !TEST end if
          ! End of test prints
        end do
      end do

    end do

    ! End of loop over orbitals
  end do
  ! End of loop over FS blocks
end do
call mma_deallocate(SBSET)

end subroutine SGMOP
