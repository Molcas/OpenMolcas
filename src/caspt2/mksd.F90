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

use definitions, only: iwp, wp
use constants, only: Half, Two
use SUPERINDEX, only: MTU
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTU, NTUES, NASHT

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), allocatable :: SD(:)
integer(kind=iwp) ISYM, NIN, NAS, NSD, ITU, ITU2, ITUABS, ITABS, IUABS, IXY, IXY2, IXYABS, IXABS, IYABS, IS11, IS21, IS12, IS22, &
                  IUTP, IXYP, IP1, IP2, IP, IXTP, IUYP, ID, ID1, ID2, IDISK
real(kind=wp) GUTXY, GXTUY, S11, S22, DUY

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,5)
  if (NIN == 0) cycle
  NAS = NTU(ISYM)
  NSD = (2*NAS*(2*NAS+1))/2
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
      IS11 = (ITU*(ITU-1))/2+IXY
      IS21 = (ITU2*(ITU2-1))/2+IXY
      IS12 = (IXY2*(IXY2-1))/2+ITU
      IS22 = (ITU2*(ITU2-1))/2+IXY2
      IUTP = IUABS+NASHT*(ITABS-1)
      IXYP = IXABS+NASHT*(IYABS-1)
      IP1 = max(IUTP,IXYP)
      IP2 = min(IUTP,IXYP)
      IP = (IP1*(IP1-1))/2+IP2
      GUTXY = Two*PREF(IP)
      IXTP = IXABS+NASHT*(ITABS-1)
      IUYP = IUABS+NASHT*(IYABS-1)
      IP1 = max(IXTP,IUYP)
      IP2 = min(IXTP,IUYP)
      IP = (IP1*(IP1-1))/2+IP2
      GXTUY = Two*PREF(IP)
      S11 = Two*GUTXY
      S22 = -GXTUY
      if (IXABS == ITABS) then
        ID1 = max(IUABS,IYABS)
        ID2 = min(IUABS,IYABS)
        ID = (ID1*(ID1-1))/2+ID2
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
