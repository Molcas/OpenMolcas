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

subroutine MKBD(DREF,NDREF,PREF,NPREF,FD,FP)

use definitions, only: iwp, wp
use constants, only: Half, Two, Four
use SUPERINDEX, only: MTU
use caspt2_global, only: ipea_shift
use caspt2_global, only: LUSBT
use EQSOLV, only: IDSMAT, IDBMAT
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTU, EASUM, NASHT, NTUES, EPSA

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(in) :: FD(NDREF), FP(NPREF)
real(kind=wp), allocatable :: BD(:), S(:), SD(:)
integer(kind=iwp) ISYM, NIN, NAS, NBD, NS2, NAS2, IDS, IDIAG, I, ITU, ITU2, ITUABS, ITABS, IUABS, IXY, IXY2, IXYABS, IXABS, IYABS, &
                  IB11, IB21, IB12, IB22, IUTP, IXYP, IP1, IP2, IP, ID, ID1, ID2, IDISK, IDT, IDU, IUYP, IXTP
real(kind=wp) ET, EX, B11, B22, ETX, FUY, DUY

! Set up the matrix BD(tuP,xyQ),P and Q are 1 or 2,
! Formulae used:
!    BD(tu1,xy1)=
!      2*Futxy + 2*(Ex+Et-A)*Gutxy + 2*dxt (Fuy + (Et-A)*Duy)
!    BD(tu2,xy1)= -BD(tu1,xy1)/2
!    BD(tu2,xy2)=
!       -Fxtuy - (Ex+Et-A)*Gxtuy + 2*dxt (Fuy + (Ex-A)*Duy)
! where A=EASUM=Sum(w) of (Ew*Dww)

! Loop over superindex symmetry.
do ISYM=1,NSYM
  NIN = NINDEP(ISYM,5)
  if (NIN == 0) cycle
  NAS = NTU(ISYM)
  NBD = (2*NAS*(2*NAS+1))/2
  if (NBD > 0) then
    call mma_allocate(BD,NBD,Label='BD')
    !GG.Nov03  Load in SD the diagonal elements of SD matrix:
    NS2 = (2*NAS*(2*NAS+1))/2
    NAS2 = 2*NAS
    call mma_allocate(S,NS2,Label='S')
    call mma_allocate(SD,NAS2,Label='SD')
    IDS = IDSMAT(ISYM,5)
    call DDAFILE(LUSBT,2,S,NS2,IDS)
    IDIAG = 0
    do I=1,NAS2
      IDIAG = IDIAG+I
      SD(I) = S(IDIAG)
    end do
    call mma_deallocate(S)
    !GG End
  end if
  do ITU=1,NAS
    ITU2 = ITU+NAS
    ITUABS = ITU+NTUES(ISYM)
    ITABS = MTU(1,ITUABS)
    IUABS = MTU(2,ITUABS)
    ET = EPSA(ITABS)
    do IXY=1,ITU
      IXY2 = IXY+NAS
      IXYABS = IXY+NTUES(ISYM)
      IXABS = MTU(1,IXYABS)
      IYABS = MTU(2,IXYABS)
      EX = EPSA(IXABS)
      IB11 = (ITU*(ITU-1))/2+IXY
      IB21 = (ITU2*(ITU2-1))/2+IXY
      IB12 = (IXY2*(IXY2-1))/2+ITU
      IB22 = (ITU2*(ITU2-1))/2+IXY2
      IUTP = IUABS+NASHT*(ITABS-1)
      IXYP = IXABS+NASHT*(IYABS-1)
      IP1 = max(IUTP,IXYP)
      IP2 = min(IUTP,IXYP)
      IP = (IP1*(IP1-1))/2+IP2
      ETX = ET+EX
      B11 = Four*(FP(IP)+(ETX-EASUM)*PREF(IP))
      IXTP = IXABS+NASHT*(ITABS-1)
      IUYP = IUABS+NASHT*(IYABS-1)
      IP1 = max(IXTP,IUYP)
      IP2 = min(IXTP,IUYP)
      IP = (IP1*(IP1-1))/2+IP2
      B22 = -Two*(FP(IP)+(ETX-EASUM)*PREF(IP))
      if (IXABS == ITABS) then
        ID1 = max(IUABS,IYABS)
        ID2 = min(IUABS,IYABS)
        ID = (ID1*(ID1-1))/2+ID2
        FUY = FD(ID)
        DUY = DREF(ID)
        B11 = B11+Two*(FUY+(ET-EASUM)*DUY)
        B22 = B22+Two*(FUY+(EX-EASUM)*DUY)
      end if
      BD(IB11) = B11
      BD(IB21) = -Half*B11
      BD(IB12) = -Half*B11
      BD(IB22) = B22
      !GG.Nov03
      if (ITU == IXY) then
        IDT = (ITABS*(ITABS+1))/2
        IDU = (IUABS*(IUABS+1))/2
        BD(IB11) = BD(IB11)+ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU)
        BD(IB22) = BD(IB22)+ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
      end if
      !GG End
    end do
  end do

  ! Write to disk
  if ((NBD > 0) .and. (NINDEP(ISYM,5) > 0)) then
    IDISK = IDBMAT(ISYM,5)
    call DDAFILE(LUSBT,1,BD,NBD,IDISK)
    call mma_deallocate(BD)
    !GG.Nov03 DisAlloc SD
    call mma_deallocate(SD)
    !GG End
  end if
end do

end subroutine MKBD
