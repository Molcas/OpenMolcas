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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKSD(DREF,NDREF,PREF,NPREF)
! Set up the matrix SD(tuP,xyQ),P and Q are 1 or 2,
! Formulae used:
!    SD(tu1,xy1)=2*(Gutxy + dxt Duy)
!    SD(tu2,xy1)= -(Gutxy + dxt Duy)
!    SD(tu2,xy2)= -Gxtuy +2*dxt Duy

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: MTU
use EQSOLV, only: IDSMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: NASHT, NINDEP, NSYM, NTU, NTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
integer(kind=iwp) :: ID, IDISK, IP, IS11, IS12, IS21, IS22, ISYM, ITABS, ITU, ITU2, ITUABS, IUABS, IUTP, IUYP, IXABS, IXTP, IXY, &
                     IXY2, IXYABS, IXYP, IYABS, NAS, NIN, NSD
real(kind=wp) :: DUY, GUTXY, GXTUY, S11, S22
real(kind=wp), allocatable :: SD(:)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,5)
  if (NIN == 0) cycle
  NAS = NTU(ISYM)
  NSD = nTri_Elem(2*NAS)
  if (NSD > 0) call mma_allocate(SD,NSD,LABEL='SD')
  do ITU=1,NAS
    ITU2 = ITU+NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    do IXY=1,ITU
      IXY2 = IXY+NAS
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      IS11 = iTri(ITU,IXY)
      IS21 = iTri(ITU2,IXY)
      IS12 = iTri(IXY2,ITU)
      IS22 = iTri(ITU2,IXY2)
      IUTP = IUABS+NASHT*(ITABS-1)
      IXYP = IXABS+NASHT*(IYABS-1)
      IP = iTri(IUTP,IXYP)
      GUTXY = Two*PREF(IP)
      IXTP = IXABS+NASHT*(ITABS-1)
      IUYP = IUABS+NASHT*(IYABS-1)
      IP = iTri(IXTP,IUYP)
      GXTUY = Two*PREF(IP)
      S11 = Two*GUTXY
      S22 = -GXTUY
      if (IXABS == ITABS) then
        ID = iTri(IUABS,IYABS)
        DUY = DREF(ID)
        S11 = S11+Two*DUY
        S22 = S22+Two*DUY
      end if
      ! SD(tu1,xy1)=2*(Gutxy + dtx Duy)
      SD(IS11) = S11
      ! SD(tu2,xy1)= -(Gutxy + dtx Duy)
      SD(IS21) = -Half*S11
      SD(IS12) = -Half*S11
      ! SD(tu2,xy2)= -Gxtuy +2*dtx Duy
      SD(IS22) = S22
    end do
  end do

  ! Write to disk
  if (NSD > 0) then
    if (NINDEP(ISYM,5) > 0) then
      IDISK = IDSMAT(ISYM,5)
      call DDAFILE(LUSBT,1,SD,NSD,IDISK)
    end if
    call mma_deallocate(SD)
  end if
end do

end subroutine MKSD
