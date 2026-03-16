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

subroutine W2SGORD1(NLEV,NVERT,NMIDV,NIPWLK,ISM,MIDLEV,MVSTA,IOCSF,NOW,IOW,IDOWN,MAW,ICS,MWS2W,MIPWLK,NLIST,KWALK,ICNUM)
! Purpose: For a wave function in Split GUGA storage structure,
! given KWALK(J,I) with J=1..MIPWLK that contains the
! complete Guga walk, as a packed array of case numbers, construct
! ICNUM(I), which is the index of this configuration in the
! Split-Guga scheme.
! Use MAWS to WLK table, MWS2W

use definitions, only: iwp
use Symmetry_Info, only: nSym => nIrrep, MUL

implicit none
integer(kind=iwp), intent(in) :: NLEV, NVERT, NMIDV, NIPWLK
integer(kind=iwp), intent(in) :: ISM(NLEV)
integer(kind=iwp), intent(in) :: MIDLEV, MVSTA
integer(kind=iwp), intent(in) :: IOCSF(NSYM,NMIDV,NSYM)
integer(kind=iwp), intent(in) :: NOW(2,NSYM,NMIDV), IOW(2,NSYM,NMIDV)
integer(kind=iwp), intent(in) :: IDOWN(NVERT,0:3), MAW(NVERT,0:3)
integer(kind=iwp), intent(out) :: ICS(NLEV)
integer(kind=iwp), intent(in) :: MWS2W(*)
integer(kind=iwp), intent(in) :: MIPWLK, NLIST
integer(kind=iwp), intent(in) :: KWALK(MIPWLK,NLIST)
integer(kind=iwp), intent(out) :: ICNUM(NLIST)
integer(kind=iwp) IC, ICONF, IDV, IDW, IOFF, ISYCI, ISYDWN, ISYUP, IUV, IUW, LDIM, LEV, MAWSD, MAWSU, MV

do ICONF=1,NLIST
  ! Unpack total walk to ICS()
  call UPKWLK(NLEV,MIPWLK,1,KWALK(1,ICONF),ICS)

  ! Follow upper walk down to MIDLEV:
  MAWSU = 0
  IUV = 1
  ISYUP = 1
  IDV = -1000000000
  do LEV=NLEV,MIDLEV+1,-1
    IC = ICS(LEV)
    if (((IC+1)/2) == 1) ISYUP = MUL(ISM(LEV),ISYUP)
    IDV = IDOWN(IUV,IC)
    MAWSU = MAWSU+MAW(IUV,IC)
    IUV = IDV
  end do
  ! We have found the midvertex number:
  MV = IDV+1-MVSTA
  ! Follow lower walk down to MIDLEV:
  ISYDWN = 1
  MAWSD = 0
  do LEV=MIDLEV,1,-1
    IC = ICS(LEV)
    if (((IC+1)/2) == 1) ISYDWN = MUL(ISM(LEV),ISYDWN)
    IDV = IDOWN(IUV,IC)
    MAWSD = MAWSD+MAW(IUV,IC)
    IUV = IDV
  end do
  IUW = MWS2W(MAWSU)
  IDW = MWS2W(MAWSD)
  ! Subtract the offsets:
  IUW = IUW-IOW(1,ISYUP,MV)/NIPWLK
  IDW = IDW-IOW(2,ISYDWN,MV)/NIPWLK
  ! Split-Guga storage scheme: an element in a set of matrices.
  ! Offset to the matrix we want is IOCSF(ISYUP,MV,ISYCI).
  ! Leading dimension=nr of upwalks in this block.
  ISYCI = MUL(ISYUP,ISYDWN)
  IOFF = IOCSF(ISYUP,MV,ISYCI)
  LDIM = NOW(1,ISYUP,MV)
  ICNUM(ICONF) = IOFF+IUW+LDIM*(IDW-1)
end do

end subroutine W2SGORD1
