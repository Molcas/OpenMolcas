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

subroutine MKBE(DREF,NDREF,FD)

use EQSOLV, only: IDBMAT, IDSMAT
use caspt2_global, only: ipea_shift, LUSBT
use caspt2_module, only: EASUM, EPSA, NAES, NASH, NINDEP, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF
real(kind=wp), intent(in) :: DREF(NDREF), FD(NDREF)
integer(kind=iwp) :: I, IBE, ID, IDIAG, IDISK, IDS, IDT, ISYM, IT, ITABS, IX, IXABS, NAS, NBE, NINM, NINP, NS
real(kind=wp) :: ET, EX, Val
real(kind=wp), allocatable :: BE(:), S(:), SD(:)

! Set up the matrix BE(t,x)
! Formula used:
!    BE(t,x)=-Ftx + (EASUM-Ex-Et)*Dtx + 2dtx Ex

do ISYM=1,NSYM
  NINP = NINDEP(ISYM,6)
  if (NINP == 0) cycle
  NINM = NINDEP(ISYM,7)
  NAS = NASH(ISYM)
  NBE = (NAS*(NAS+1))/2
  if (NBE > 0) then
    call mma_allocate(BE,NBE,LABEL='BE')
    !GG.Nov03  Load in SD the diagonal elements of SE matrix:
    NS = (NAS*(NAS+1))/2
    call mma_allocate(S,NS,Label='S')
    call mma_allocate(SD,NAS,Label='SD')
    IDS = IDSMAT(ISYM,6)
    call DDAFILE(LUSBT,2,S,NS,IDS)
    IDIAG = 0
    do I=1,NAS
      IDIAG = IDIAG+I
      SD(I) = S(IDIAG)
    end do
    call mma_deallocate(S)
    !GG End
  end if
  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    ET = EPSA(ITABS)
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      EX = EPSA(IXABS)
      IBE = (IT*(IT-1))/2+IX
      ID = (ITABS*(ITABS-1))/2+IXABS
      Val = -FD(ID)+(EASUM-EX-ET)*DREF(ID)
      if (ITABS == IXABS) Val = Val+Two*EX
      !GG.Nov03
      if (IT == IX) then
        IDT = (ITABS*(ITABS+1))/2
        Val = Val+ipea_shift*Half*DREF(IDT)*SD(IT)
      end if
      !GG End
      BE(IBE) = Val
    end do
  end do

  ! Write to disk
  if ((NBE > 0) .and. (NINDEP(ISYM,6) > 0)) then
    IDISK = IDBMAT(ISYM,6)
    call DDAFILE(LUSBT,1,BE,NBE,IDISK)
    if ((NINM > 0) .and. (NINDEP(ISYM,7) > 0)) then
      IDISK = IDBMAT(ISYM,7)
      call DDAFILE(LUSBT,1,BE,NBE,IDISK)
    end if
    call mma_deallocate(BE)
    !GG.Nov03 DisAlloc SD
    call mma_deallocate(SD)
    !GG End
  end if
end do

end subroutine MKBE
