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

subroutine PRWF1(SGS,CIS,NLEV,NMIDV,ISM,ICS,NOCSF,IOCSF,NOW,IOW,ISYCI,CI,CITHR)

use definitions, only: iwp, wp, u6
use gugx, only: SGStruct, CIStruct
use Symmetry_Info, only: nSym => nIrrep, MUL

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in) :: NLEV, NMIDV
integer(kind=iwp), intent(in) :: ISM(NLEV)
integer(kind=iwp), intent(out) :: ICS(NLEV)
integer(kind=iwp), intent(in) :: NOCSF(NSYM,NMIDV,NSYM), IOCSF(NSYM,NMIDV,NSYM)
integer(kind=iwp), intent(in) :: NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
integer(kind=iwp), intent(in) :: ISYCI
real(kind=wp), intent(in) :: CI(*), CITHR
logical(kind=iwp), parameter :: SGINFO = .true.
character(len=80) LINE
character(len=1) :: CODE(0:3) = ['0','u','d','2']
integer(kind=iwp) IC1, ICDPOS, ICDWN, ICONF, ICUP, ICUPOS, IDW0, IDWN, IDWNSV, ISY, ISYDWN, ISYUP, IUP, IUW0, K, KNXT, KOCLAB, &
                  KOCSZ, KPAD1, KPAD2, LEV, MIDLEV, MV, NCI, NDWN, NIPWLK, NNN, NUP
real(kind=wp) COEF

! -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
! -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
! -- THE   U P P E R   PART OF THE WALK.
! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.

MIDLEV = SGS%MidLev
NIPWLK = CIS%nIpWlk

! Size of occup/spin coupling part of line:
write(u6,*) ' Occupation of active orbitals, and spin coupling'
write(u6,*) ' of open shells. (u,d: Spin up or down).'
LINE = ''
K = 0
ISY = 0
do LEV=1,NLEV
  if (ISY /= ISM(LEV)) then
    ISY = ISM(LEV)
    K = K+1
  end if
  K = K+1
end do
KOCLAB = 10
KOCSZ = max(K,KOCLAB)
KPAD1 = (KOCSZ-KOCLAB)/2
KPAD2 = (KOCSZ-K)/2
if (SGINFO) write(u6,*) ' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
LINE(1:7) = '  Conf '
K = 7
if (SGINFO) then
  LINE(8:22) = '  SGUGA info   '
  K = 22
end if
LINE(K+KPAD1:K+KPAD1+9) = 'Occupation'
K = K+KOCSZ
LINE(K:K+23) = '       Coef       Weight'
write(u6,*) LINE
LINE = ''

! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
do MV=1,NMIDV
  do ISYUP=1,NSYM
    NCI = NOCSF(ISYUP,MV,ISYCI)
    if (NCI == 0) cycle
    NUP = NOW(1,ISYUP,MV)
    ISYDWN = MUL(ISYUP,ISYCI)
    NDWN = NOW(2,ISYDWN,MV)
    ICONF = IOCSF(ISYUP,MV,ISYCI)
    IUW0 = 1-NIPWLK+IOW(1,ISYUP,MV)
    IDW0 = 1-NIPWLK+IOW(2,ISYDWN,MV)
    IDWNSV = 0
    do IDWN=1,NDWN
      do IUP=1,NUP
        ICONF = ICONF+1
        COEF = CI(ICONF)
        ! -- SKIP OR PRINT IT OUT?
        if (abs(COEF) < CITHR) cycle
        if (IDWNSV /= IDWN) then
          ICDPOS = IDW0+IDWN*NIPWLK
          ICDWN = CIS%ICase(ICDPOS)
          ! -- UNPACK LOWER WALK.
          NNN = 0
          do LEV=1,MIDLEV
            NNN = NNN+1
            if (NNN == 16) then
              NNN = 1
              ICDPOS = ICDPOS+1
              ICDWN = CIS%ICase(ICDPOS)
            end if
            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1
          end do
          IDWNSV = IDWN
        end if
        ICUPOS = IUW0+NIPWLK*IUP
        ICUP = CIS%ICase(ICUPOS)
        ! -- UNPACK UPPER WALK:
        NNN = 0
        do LEV=MIDLEV+1,NLEV
          NNN = NNN+1
          if (NNN == 16) then
            NNN = 1
            ICUPOS = ICUPOS+1
            ICUP = CIS%ICase(ICUPOS)
          end if
          IC1 = ICUP/4
          ICS(LEV) = ICUP-4*IC1
          ICUP = IC1
        end do
        ! -- PRINT IT!
        write(LINE(1:7),'(I6,1X)') ICONF
        K = 7
        if (SGINFO) then
          LINE(K+1:K+1) = '('
          write(LINE(K+2:K+3),'(I2)') MV
          LINE(K+4:K+4) = ':'
          write(LINE(K+5:K+5),'(I1)') ISYUP
          LINE(K+6:K+6) = ':'
          write(LINE(K+7:K+9),'(I3)') IUP
          LINE(K+10:K+10) = '/'
          write(LINE(K+11:K+13),'(I3)') IDWN
          LINE(K+14:K+14) = ')'
          K = K+14
        end if
        KNXT = K+KOCSZ
        K = K+KPAD2
        ISY = 0
        do LEV=1,NLEV
          if (ISY /= ISM(LEV)) then
            ISY = ISM(LEV)
            K = K+1
            LINE(K:K) = ' '
          end if
          K = K+1
          LINE(K:K) = CODE(ICS(LEV))
        end do
        K = KNXT
        K = K+1
        LINE(K:K+4) = '     '
        K = K+5
        write(LINE(K:K+7),'(F8.5)') COEF
        K = K+8
        LINE(K:K+4) = '     '
        K = K+5
        write(LINE(K:K+7),'(F8.5)') COEF**2
        write(u6,*) LINE(1:K+7)
      end do
    end do
  end do
end do
write(u6,*)
write(u6,*) repeat('*',80)

end subroutine PRWF1
