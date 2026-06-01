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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine PRWF1_CP2(NOCSF,IOCSF,NOW,IOW,ISYCI,CI,mCI,THR,nMidV)

use Symmetry_Info, only: Mul
use sguga, only: CIS, SGS
use caspt2_module, only: ISPIN, NSYM, PRSD
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nMidV, NOCSF(NSYM,NMIDV,NSYM), IOCSF(NSYM,NMIDV,NSYM), NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV), &
                                 ISYCI, mCI
real(kind=wp), intent(in) :: CI(mCI), THR
integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICONF, ICUP, ICUPOS, IDW0, IDWN, IDWNSV, IMS, ISY, ISYDWN, ISYUP, IUP, IUW0, K, LENCSF, &
                     LEV, MV, NCI, NDWN, nIpWlk, nLev, NNN, NUP
real(kind=wp) :: COEF
character(len=256) :: LINE
integer(kind=iwp), allocatable :: ICS(:), LEX(:)
character, parameter :: CODE(0:3) = ['0','u','d','2']

nLev = SGS%nLev
nIpWlk = CIS%nIpWlk

! NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
! ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
! THE   U P P E R   PART OF THE WALK.

call mma_allocate(ICS,NLEV,Label='ICS')

! SVC: set up a CSF string length as LENCSF
LINE = ' '
LENCSF = 0
ISY = 0
do LEV=1,NLEV
  if (ISY /= SGS%ISM(LEV)) then
    ISY = SGS%ISM(LEV)
    LENCSF = LENCSF+1
  end if
  LENCSF = LENCSF+1
end do
LENCSF = min(LENCSF,256)
LENCSF = max(LENCSF,10)

! Size of occup/spin coupling part of line:
write(u6,*) ' Occupation of active orbitals, and spin coupling'
write(u6,*) ' of open shells. (u,d: Spin up or down).'
write(u6,*) ' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
LINE(1:10) = 'Occupation'
write(u6,'(2X,A10,2X,A16,2X,A,2(2X,A13))') 'Conf','SGUGA info      ',LINE(1:LENCSF),'Coefficient','Weight'

! SVC2010:
! allocate scratch memory for determinant expansion
if (PRSD) call mma_allocate(LEX,NLEV,LABEL='LEX')

LINE = ''

! THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
! WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
do MV=1,NMIDV
  do ISYUP=1,NSYM
    NCI = NOCSF(ISYUP,MV,ISYCI)
    if (NCI == 0) cycle
    NUP = NOW(1,ISYUP,MV)
    ISYDWN = Mul(ISYUP,ISYCI)
    NDWN = NOW(2,ISYDWN,MV)
    ICONF = IOCSF(ISYUP,MV,ISYCI)
    IUW0 = 1-NIPWLK+IOW(1,ISYUP,MV)
    IDW0 = 1-NIPWLK+IOW(2,ISYDWN,MV)
    IDWNSV = 0
    do IDWN=1,NDWN
      do IUP=1,NUP
        ICONF = ICONF+1
        COEF = CI(ICONF)
        ! SKIP OR PRINT IT OUT?
        if (abs(COEF) < THR) cycle
        if (IDWNSV /= IDWN) then
          ICDPOS = IDW0+IDWN*NIPWLK
          ICDWN = CIS%ICASE(ICDPOS)
          ! UNPACK LOWER WALK.
          NNN = 0
          do LEV=1,SGS%MIDLEV
            NNN = NNN+1
            if (NNN == 16) then
              NNN = 1
              ICDPOS = ICDPOS+1
              ICDWN = CIS%ICASE(ICDPOS)
            end if
            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1
          end do
          IDWNSV = IDWN
        end if
        ICUPOS = IUW0+NIPWLK*IUP
        ICUP = CIS%ICASE(ICUPOS)
        ! UNPACK UPPER WALK:
        NNN = 0
        do LEV=SGS%MIDLEV+1,NLEV
          NNN = NNN+1
          if (NNN == 16) then
            NNN = 1
            ICUPOS = ICUPOS+1
            ICUP = CIS%ICASE(ICUPOS)
          end if
          IC1 = ICUP/4
          ICS(LEV) = ICUP-4*IC1
          ICUP = IC1
        end do
        ! PRINT IT!
        K = 0
        ISY = 0
        do LEV=1,NLEV
          if (ISY /= SGS%ISM(LEV)) then
            ISY = SGS%ISM(LEV)
            K = K+1
            LINE(K:K) = ' '
          end if
          K = K+1
          LINE(K:K) = CODE(ICS(LEV))
        end do
        write(u6,'(2X,I10,2X,"(",I2,":",I1,":",I4,"/",I4,")",2X,A,2(2X,F13.6))') ICONF,MV,ISYUP,IUP,IDWN,LINE(1:LENCSF),COEF,COEF**2
        ! SVC2010 experimental: add determinant expansion
        if (PRSD) then
          ! Specify projected spin in half integer units
          ! Default: use maximum spin projection
          IMS = ISPIN-1
          write(u6,*)
          call EXPCSF(ICS,NLEV,IMS,LEX,coef,0)
          write(u6,*)
        end if
      end do
    end do
  end do
end do

call mma_deallocate(ICS)

! SVC2010: free scratch for determinant expansion
if (PRSD) call mma_deallocate(LEX)
write(u6,*)

end subroutine PRWF1_CP2
