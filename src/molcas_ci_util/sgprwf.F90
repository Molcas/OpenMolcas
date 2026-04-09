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

subroutine SGPRWF(SGS,CIS,LSYM,PRWTHR,iSpin,CI,lCI,KeyPRSD,LUVECDET)
! PURPOSE: PRINT THE WAVEFUNCTION (SPIN COUPLING AND OCCUPATIONS)
!
! NOTE:    THIS ROUTINE USES THE SPLIT GRAPH GUGA CONVENTION, I.E.,
!          CI BLOCKS ARE MATRICES CI(I,J), WHERE THE  FIRST INDEX
!          REFERS TO THE UPPER PART OF THE WALK.

use gugx, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp), intent(in) :: LSYM, iSpin, lCI, LUVECDET
real(kind=wp), intent(in) :: PRWTHR, CI(lCI)
logical(kind=iwp), intent(in) :: KeyPrSD
integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICONF, ICS(50), ICUP, ICUPOS, IDW0, IDWN, IDWNSV, IMS, iOff, ISYDWN, ISYM, ISYUP, IUP, &
                     IUW0, LEV, MV, NCI, NDWN, NNN, NUP
real(kind=wp) :: COEF
character(len=400) :: Line
integer(kind=iwp), allocatable :: Lex(:)

! RECONSTRUCT THE CASE LIST

if (.not. allocated(CIS%iCase)) call MKCLIST(SGS,CIS)

! scratch for determinant expansion
if (KeyPRSD) call mma_allocate(LEX,SGS%nLev,Label='LEX')

Line(1:16) = '      conf/sym  '
iOff = 16
iSym = SGS%ISm(1)
do Lev=1,SGS%nLev
  if (SGS%ISm(Lev) /= iSym) iOff = iOff+1
  write(Line(iOff+Lev:),'(I1)') SGS%ISm(Lev)
  if (SGS%ISm(Lev) /= iSym) iSym = SGS%ISm(Lev)
end do
iOff = iOff+SGS%nLev+3
Line(iOff:iOff+15) = '   Coeff  Weight'
write(u6,'(A)') Line(1:iOff+15)
Line = ' '

! ENTER THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
! WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.

do MV=1,CIS%nMidV
  do ISYUP=1,SGS%nSym
    NCI = CIS%NOCSF(ISYUP,MV,LSYM)
    if (NCI == 0) cycle
    NUP = CIS%NOW(1,ISYUP,MV)
    ISYDWN = 1+ieor(ISYUP-1,LSYM-1)
    NDWN = CIS%NOW(2,ISYDWN,MV)
    ICONF = CIS%IOCSF(ISYUP,MV,LSYM)
    IUW0 = 1-CIS%nIpWlk+CIS%IOW(1,ISYUP,MV)
    IDW0 = 1-CIS%nIpWlk+CIS%IOW(2,ISYDWN,MV)
    IDWNSV = 0
    do IDWN=1,NDWN
      do IUP=1,NUP
        ICONF = ICONF+1
        COEF = CI(ICONF)
        ! SKIP OR PRINT IT OUT?
        if (abs(COEF) < PRWTHR) cycle
        if (IDWNSV /= IDWN) then
          ICDPOS = IDW0+IDWN*CIS%nIpWlk
          ICDWN = CIS%iCase(ICDPOS)
          ! UNPACK LOWER WALK.
          NNN = 0
          do LEV=1,SGS%MidLev
            NNN = NNN+1
            if (NNN == 16) then
              NNN = 1
              ICDPOS = ICDPOS+1
              ICDWN = CIS%iCase(ICDPOS)
            end if
            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1
          end do
          IDWNSV = IDWN
        end if
        ICUPOS = IUW0+CIS%nIpWlk*IUP
        ICUP = CIS%iCase(ICUPOS)
        ! UNPACK UPPER WALK:
        NNN = 0
        do LEV=SGS%MidLev+1,SGS%nLev
          NNN = NNN+1
          if (NNN == 16) then
            NNN = 1
            ICUPOS = ICUPOS+1
            ICUP = CIS%iCase(ICUPOS)
          end if
          IC1 = ICUP/4
          ICS(LEV) = ICUP-4*IC1
          ICUP = IC1
        end do
        ! PRINT IT!
        write(Line(1:),'(I8)') iConf
        iOff = 10
        iSym = SGS%ISm(1)
        do Lev=1,SGS%nLev
          if (SGS%ISm(Lev) /= iSym) iOff = iOff+1

          select case (ICS(Lev))
            case (3)
              write(Line(iOff+Lev:),'(A1)') '2'
            case (2)
              write(Line(iOff+Lev:),'(A1)') 'd'
            case (1)
              write(Line(iOff+Lev:),'(A1)') 'u'
            case (0)
              write(Line(iOff+Lev:),'(A1)') '0'
            case default
              call Abend()
          end select

          if (SGS%ISm(Lev) /= iSym) iSym = SGS%ISm(Lev)

        end do
        iOff = iOff+SGS%nLev+3
        write(Line(iOff:),'(2F8.5)') COEF,COEF**2
        write(u6,'(6X,A)') Line(1:iOff+15)
        if (KeyPRSD) then
          ! use maximum spin projection value
          IMS = ISPIN-1
          write(u6,*)
          call EXPCSF(ICS,SGS%nLev,IMS,LEX,coef,LuVecDet)
          write(u6,*)
        end if
        Line = ' '
      end do
    end do
  end do
end do

! free memory for determinant expansion
if (KeyPRSD) call mma_deallocate(LEX)

end subroutine SGPRWF
