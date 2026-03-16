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

subroutine NEWSCTAB(MINOP,MAXOP,MLTPL,MS2,ICASE)

use definitions, only: iwp, u6
use stdalloc, only: mma_allocate
use rassi_global_arrays, only: TRANS1, TRANS2, TRANS
use rassi_global_arrays, only: SPNTAB1, SPNTAB2, SPNTAB

implicit none
integer(kind=iwp), intent(in) :: MINOP, MAXOP, MLTPL, MS2, ICASE
integer(kind=iwp), parameter :: ASPIN = 1, BSPIN = 0
integer(kind=iwp), parameter :: UPCPL = 1, DWNCPL = 0, NULLPTR = -1
integer(kind=iwp) IBLK, IOPEN, KSPCPL, KSPDET, LTRANS, LTRANS0, NA, NBLK, NCP, ND, NSPCPL, NSPDET, NTAB, NTRANS, NB
integer(kind=iwp), external :: NGENE, NOVERM

if ((MLTPL-MS2 < 1) .or. (MLTPL+MS2 < 1)) then
  write(u6,*) 'NewSCTab: Contradictory values of MLTPL vs. MS2.'
  write(u6,*) 'The function was invoked with the following arguments:'
  write(u6,'(1X,A,I9)') ' MINOP:',MINOP
  write(u6,'(1X,A,I9)') ' MAXOP:',MAXOP
  write(u6,'(1X,A,I9)') ' MLTPL:',MLTPL
  write(u6,'(1X,A,I9)') ' MS2  :',MS2
  call ABEND()
end if

! Run through construction loop twice. First get size of table
! and total nr of transformation coefficients:
NSPCPL = 0
NSPDET = 0
NTRANS = 0
IBLK = 0
do IOPEN=MINOP,MAXOP
  IBLK = IBLK+1
  NCP = NGENE(IOPEN,MLTPL)
  if (NCP == 0) cycle
  NA = (IOPEN+MS2)/2
  ND = NOVERM(IOPEN,NA)
  NSPCPL = NSPCPL+IOPEN*NCP
  NSPDET = NSPDET+IOPEN*ND
  NTRANS = NTRANS+ND*NCP
end do
NBLK = IBLK
NTAB = 8+6*NBLK+NSPCPL+NSPDET
! The spincoupling table will be NTAB integers long, and there
! will be NTRANS spin-coupling coefficients (real*8).
! The table consists of 8 header words, then an array (6,NBLK)
! with pointers and sizes to the spin coupling, spin determinant
! and spin-coupling coefficient arrays, and finally the
! spin coupling and spin determinant arrays themselves.
! The transformation coefficients are real*8 data and stored
! in a separate array.
select case (ICASE)
  case (1)
    call mma_allocate(SPNTAB1,NTAB,Label='SPNTAB1')
    SPNTAB => SPNTAB1(:)
    call mma_allocate(TRANS1,NTRANS,Label='TRANS1')
    TRANS => TRANS1(:)
  case (2)
    call mma_allocate(SPNTAB2,NTAB,Label='SPNTAB2')
    SPNTAB => SPNTAB2(:)
    call mma_allocate(TRANS2,NTRANS,Label='TRANS2')
    TRANS => TRANS2(:)
  case DEFAULT
end select
KSPCPL = 9+6*NBLK
KSPDET = KSPCPL+NSPCPL
! Table size
SPNTAB(1) = NTAB
! Table type identifier
SPNTAB(2) = 47
! Spin multiplicity
SPNTAB(3) = MLTPL
! Spin projection
SPNTAB(4) = MS2
! Min and max nr of open shells
SPNTAB(5) = MINOP
SPNTAB(6) = MAXOP
! Associated workspace array for Re*8 data (transf matrices)
SPNTAB(7) = -1 ! not used
SPNTAB(8) = NTRANS
! Individual information for each separate nr of open shells:
NTAB = 6
NTRANS = 0
IBLK = 0
LTRANS0 = 1
LTRANS = 1
do IOPEN=MINOP,MAXOP
  IBLK = IBLK+1
  NCP = NGENE(IOPEN,MLTPL)
  NA = (IOPEN+MS2)/2
  NB = IOPEN-NA
  if (NCP > 0) then
    ND = NOVERM(IOPEN,NA)
    SPNTAB(9+(IBLK-1)*6) = IOPEN
    SPNTAB(10+(IBLK-1)*6) = NCP
    SPNTAB(11+(IBLK-1)*6) = ND
    ! Compute spin couplings:
    call PROTOCSF(IOPEN,MLTPL,NCP,SPNTAB(KSPCPL:))
    SPNTAB(12+(IBLK-1)*6) = KSPCPL
    ! Compute spin determinants:
    call PROTOSD(NA,NB,ND,SPNTAB(KSPDET:))
    SPNTAB(13+(IBLK-1)*6) = KSPDET
    ! Compute spin coupling coefficients:
    call PROTOT(IOPEN,ND,SPNTAB(KSPDET:),NCP,SPNTAB(KSPCPL:),TRANS(LTRANS:))
    SPNTAB(14+(IBLK-1)*6) = LTRANS-LTRANS0+1
    KSPCPL = KSPCPL+IOPEN*NCP
    KSPDET = KSPDET+IOPEN*ND
    LTRANS = LTRANS+ND*NCP
  else
    SPNTAB(9+(IBLK-1)*6) = IOPEN
    SPNTAB(10+(IBLK-1)*6) = 0
    SPNTAB(11+(IBLK-1)*6) = 0
    SPNTAB(12+(IBLK-1)*6) = NULLPTR
    SPNTAB(13+(IBLK-1)*6) = NULLPTR
    SPNTAB(14+(IBLK-1)*6) = NULLPTR
  end if
end do

nullify(TRANS,SPNTAB)

end subroutine NEWSCTAB
