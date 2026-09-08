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

subroutine SG_PRWF(iState,ISYCI,CITHR,iSpin,CI,lCI,KeyPRSD,LuVecDet)

use sguga, only: CIS, nPack, SGS
use Symmetry_Info, only: MUL, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iState, ISYCI, iSpin, lCI, LuVecDet
logical(kind=iwp), intent(in) :: KeyPRSD
real(kind=wp), intent(in) :: CITHR, CI(lCI)
integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICONF, ICUP, ICUPOS, IDW0, IDWN, IDWNSV, IMS, ISY, ISYDWN, ISYUP, IUP, IUW0, K, KNXT, &
                     KOCLAB, KOCSZ, KPAD1, KPAD2, LEV, MV, NCI, NDWN, NNN, NUP
real(kind=wp) :: COEF
logical(kind=iwp) :: PrSym
character(len=200) :: LINE
integer(kind=iwp), allocatable :: ICS(:), Lex(:)
logical(kind=iwp), parameter :: SGINFO = .false.
character, parameter :: CODE(0:3) = ['0','u','d','2']

! scratch for determinant expansion
if (KeyPRSD) call mma_allocate(LEX,SGS(iState)%nLev,Label='LEX')

! -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
! -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
! -- THE   U P P E R   PART OF THE WALK.
! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.

call mma_allocate(ICS,SGS(iState)%nLev,Label='ICS')

! Size of occup/spin coupling part of line:
write(u6,*)
write(u6,100) 'Occupation of active orbitals, and spin coupling of open shells. (u,d: Spin up or down).'
write(u6,*)
LINE = ''
ISY = 0
PrSym = .false.
K = SGS(iState)%nLev
do LEV=1,SGS(iState)%nLev
  if (ISY /= SGS(iState)%ISM(LEV)) then
    ISY = SGS(iState)%ISM(LEV)
    if (ISY > 1) PrSym = .true.
    K = K+1
  end if
end do
KOCLAB = len(' Occupation')
KOCSZ = max(K,KOCLAB)
KPAD1 = KOCSZ-KOCLAB-(KOCSZ-KOCLAB)/2
KPAD2 = KOCSZ-K-(KOCSZ-K)/2
if (SGINFO) write(u6,100) 'SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
K = 0
LINE(K+1:K+7) = '   Conf'
K = K+7
if (SGINFO) then
  LINE(K+1:K+15) = '   SGUGA info  '
  K = K+15
end if
LINE(K+KPAD1+1:K+KPAD1+11) = ' Occupation'
K = K+KOCSZ
LINE(K+1:K+17) = '    Coeff  Weight'
write(u6,100) trim(LINE)
if (PrSym) then
  LINE = ''
  K = 0
  LINE(K+1:K+7) = '    Sym'
  K = K+7
  if (SGINFO) K = K+15
  K = K+KPAD2
  ISY = 0
  do LEV=1,SGS(iState)%nLev
    if (ISY /= SGS(iState)%ISM(LEV)) then
      ISY = SGS(iState)%ISM(LEV)
      K = K+1
    end if
    K = K+1
    write(LINE(K:K),'(I1)') ISY
  end do
  write(u6,100) trim(LINE)
end if

! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
!    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
do MV=1,CIS(iState)%nMidV
  do ISYUP=1,nIrrep
    NCI = CIS(iState)%NOCSF(ISYUP,MV,ISYCI)
    if (NCI == 0) cycle
    NUP = CIS(iState)%NOW(1,ISYUP,MV)
    ISYDWN = MUL(ISYUP,ISYCI)
    NDWN = CIS(iState)%NOW(2,ISYDWN,MV)
    ICONF = CIS(iState)%IOCSF(ISYUP,MV,ISYCI)
    IUW0 = 1-CIS(iState)%nIpWlk+CIS(iState)%IOW(1,ISYUP,MV)
    IDW0 = 1-CIS(iState)%nIpWlk+CIS(iState)%IOW(2,ISYDWN,MV)
    IDWNSV = 0
    do IDWN=1,NDWN
      do IUP=1,NUP
        ICONF = ICONF+1
        COEF = CI(ICONF)
        ! -- SKIP OR PRINT IT OUT?
        if (abs(COEF) < CITHR) cycle
        if (IDWNSV /= IDWN) then
          ICDPOS = IDW0+IDWN*CIS(iState)%nIpWlk
          ICDWN = CIS(iState)%ICase(ICDPOS)
          ! -- UNPACK LOWER WALK.
          NNN = 0
          do LEV=1,SGS(iState)%MidLev
            NNN = NNN+1
            if (NNN == nPack+1) then
              NNN = 1
              ICDPOS = ICDPOS+1
              ICDWN = CIS(iState)%ICase(ICDPOS)
            end if
            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1
          end do
          IDWNSV = IDWN
        end if
        ICUPOS = IUW0+CIS(iState)%nIpWlk*IUP
        ICUP = CIS(iState)%ICase(ICUPOS)
        ! -- UNPACK UPPER WALK:
        NNN = 0
        do LEV=SGS(iState)%MidLev+1,SGS(iState)%nLev
          NNN = NNN+1
          if (NNN == nPack+1) then
            NNN = 1
            ICUPOS = ICUPOS+1
            ICUP = CIS(iState)%ICase(ICUPOS)
          end if
          IC1 = ICUP/4
          ICS(LEV) = ICUP-4*IC1
          ICUP = IC1
        end do
        ! -- PRINT IT!
        LINE = ''
        K = 0
        write(LINE(K+1:K+7),'(I7)') ICONF
        K = K+7
        if (SGINFO) then
          write(LINE(K+1:K+15),'(1x,"(",I2,":",I1,":",I3,"/",I3,")")') MV,ISYUP,IUP,IDWN
          K = K+15
        end if
        KNXT = K+KOCSZ
        K = K+KPAD2
        ISY = 0
        do LEV=1,SGS(iState)%nLev
          if (ISY /= SGS(iState)%ISM(LEV)) then
            ISY = SGS(iState)%ISM(LEV)
            K = K+1
          end if
          K = K+1
          LINE(K:K) = CODE(ICS(LEV))
        end do
        K = KNXT
        write(LINE(K+1:K+17),'(1x,F8.5,F8.5)') COEF,COEF**2
        K = K+17
        write(u6,100) trim(LINE)
        if (KeyPRSD) then
          ! use maximum spin projection value
          IMS = ISPIN-1
          write(u6,*)
          call EXPCSF(ICS,SGS(iState)%nLev,IMS,LEX,coef,LuVecDet)
          write(u6,*)
        end if
      end do
    end do
  end do
end do
write(u6,*)
write(u6,100) repeat('*',120)

call mma_deallocate(ICS)

! free memory for determinant expansion
if (KeyPRSD) call mma_deallocate(LEX)

100 format(6x,a)

end subroutine SG_PRWF
