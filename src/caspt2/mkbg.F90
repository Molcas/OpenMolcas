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

subroutine MKBG(DREF,NDREF,FD)

use Index_Functions, only: iTri, nTri_Elem
use EQSOLV, only: IDBMAT, IDSMAT
use caspt2_global, only: ipea_shift, LUSBT
use general_data, only: NASH
use caspt2_module, only: EASUM, NAES, NINDEP, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF
real(kind=wp), intent(in) :: DREF(NDREF), FD(NDREF)
integer(kind=iwp) :: I, IBG, ID, IDIAG, IDISK, IDS, IDT, ISYM, IT, ITABS, IX, IXABS, NAS, NBG, NINM, NINP
real(kind=wp) :: Val
real(kind=wp), allocatable :: BG(:), S(:), SD(:)

! Set up the matrix BG(t,x)
! Formula used:
! BG(t,x)= Ftx -EASUM*Dtx

do ISYM=1,NSYM
  NINP = NINDEP(ISYM,10)
  if (NINP == 0) cycle
  NINM = NINDEP(ISYM,11)
  NAS = NASH(ISYM)
  NBG = nTri_Elem(NAS)
  if (NBG > 0) then
    call mma_Allocate(BG,NBG,LABEL='BG')
    !GG.Nov03  Load in SD the diagonal elements of SG matrix:
    call mma_allocate(S,NBG,Label='S')
    call mma_allocate(SD,NAS,Label='SD')
    IDS = IDSMAT(ISYM,10)
    call DDAFILE(LUSBT,2,S,NBG,IDS)
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
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      IBG = iTri(IT,IX)
      ID = iTri(ITABS,IXABS)
      !GG.Nov03
      Val = FD(ID)-EASUM*DREF(ID)
      if (IT == IX) then
        IDT = nTri_Elem(ITABS)
        Val = Val+ipea_shift*half*(two-DREF(IDT))*SD(IT)
      end if
      BG(IBG) = Val
      !BG(BG) = FD(ID)-EASUM*DREF(ID)
      !GG End
    end do
  end do

  ! Write to disk, and save size and address.
  if (NBG > 0) then
    if (NINDEP(ISYM,10) > 0) then
      IDISK = IDBMAT(ISYM,10)
      call DDAFILE(LUSBT,1,BG,NBG,IDISK)
    end if
    if ((NINM > 0) .and. (NINDEP(ISYM,11) > 0)) then
      IDISK = IDBMAT(ISYM,11)
      call DDAFILE(LUSBT,1,BG,NBG,IDISK)
    end if
    !GG.Nov03 DisAlloc SD
    call mma_deallocate(SD)
    !GG End
    call mma_deallocate(BG)
  end if
end do

end subroutine MKBG
