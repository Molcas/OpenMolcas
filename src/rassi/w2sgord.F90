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

subroutine W2SGORD(SGS,CIS,MWS2W,NLIST,KWALK,ICNUM)
! Purpose: Given a list of bit-packed total walks,
! translate this into a list of elements of a CI array.
! MXCPI= Max nr of case numbers packed in one integer.
!
! For a wave function in Split GUGA storage structure,
! given KWALK(J,I) with J=1..MIPWLK that contains the
! complete Guga walk, as a packed array of case numbers, construct
! ICNUM(I), which is the index of this configuration in the
! Split-Guga scheme.
! Use MAWS to WLK table, MWS2W

use gugx, only: CIStruct, SGStruct
use Symmetry_Info, only: MUL
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in) :: MWS2W(*), NLIST, KWALK(*)
integer(kind=iwp), intent(out) :: ICNUM(NLIST)
integer(kind=iwp) :: IC, ICONF, IDV, IDW, IOFF, ISYCI, ISYDWN, ISYUP, IUV, IUW, LDIM, LEV, MAWSD, MAWSU, MIPWLK, MV
integer(kind=iwp), allocatable :: ICS(:)
integer(kind=iwp), parameter :: MXCPI = 15

! Nr of integers used to store each total walk:
MIPWLK = 1+(SGS%nLev-1)/MXCPI
! Allocate scratch space for case numbers:
call mma_allocate(ICS,SGS%nLev,Label='ICS')

do ICONF=1,NLIST
  ! Unpack total walk to ICS()
  call UPKWLK(SGS%nLev,MIPWLK,1,KWALK((ICONF-1)*MIPWLK+1),ICS)

  ! Follow upper walk down to SGS%MidLev:
  MAWSU = 0
  IUV = 1
  ISYUP = 1
  IDV = -1000000000
  do LEV=SGS%nLev,SGS%MidLev+1,-1
    IC = ICS(LEV)
    if (((IC+1)/2) == 1) ISYUP = MUL(SGS%ISM(LEV),ISYUP)
    IDV = SGS%DOWN(IUV,IC)
    MAWSU = MAWSU+SGS%MAW(IUV,IC)
    IUV = IDV
  end do
  ! We have found the midvertex number:
  MV = IDV+1-SGS%MVSta
  ! Follow lower walk down to SGS%MidLev:
  ISYDWN = 1
  MAWSD = 0
  do LEV=SGS%MidLev,1,-1
    IC = ICS(LEV)
    if (((IC+1)/2) == 1) ISYDWN = MUL(SGS%ISM(LEV),ISYDWN)
    IDV = SGS%DOWN(IUV,IC)
    MAWSD = MAWSD+SGS%MAW(IUV,IC)
    IUV = IDV
  end do
  IUW = MWS2W(MAWSU)
  IDW = MWS2W(MAWSD)
  ! Subtract the offsets:
  IUW = IUW-CIS%IOW(1,ISYUP,MV)/CIS%nIpWlk
  IDW = IDW-CIS%IOW(2,ISYDWN,MV)/CIS%nIpWlk
  ! Split-Guga storage scheme: an element in a set of matrices.
  ! Offset to the matrix we want is CIS%IOCSF(ISYUP,MV,ISYCI).
  ! Leading dimension=nr of upwalks in this block.
  ISYCI = MUL(ISYUP,ISYDWN)
  IOFF = CIS%IOCSF(ISYUP,MV,ISYCI)
  LDIM = CIS%NOW(1,ISYUP,MV)
  ICNUM(ICONF) = IOFF+IUW+LDIM*(IDW-1)
end do

call mma_deallocate(ICS)

end subroutine W2SGORD
