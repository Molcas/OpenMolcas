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

subroutine MKSG(DREF,NDREF)
! Set up the matrix SG(t,x)
! Formula used:
!    SG(t,x)= Dtx

use Index_Functions, only: iTri, nTri_Elem
use EQSOLV, only: IDSMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: NAES, NASH, NINDEP, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF
real(kind=wp), intent(in) :: DREF(NDREF)
integer(kind=iwp) :: ID, IDISK, ISG, ISYM, IT, ITABS, IX, IXABS, NAS, NINM, NINP, NSG
real(kind=wp), allocatable :: SG(:)

do ISYM=1,NSYM
  NINP = NINDEP(ISYM,10)
  if (NINP == 0) cycle
  NINM = NINDEP(ISYM,11)
  NAS = NASH(ISYM)
  NSG = nTri_Elem(NAS)
  if (NSG > 0) call mma_allocate(SG,NSG,Label='SG')
  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    do IX=1,IT
      IXABS = IX+NAES(ISYM)
      ISG = iTri(IT,IX)
      ID = iTri(ITABS,IXABS)
      SG(ISG) = DREF(ID)
    end do
  end do

  ! Write to disk
  if ((NSG > 0) .and. (NINDEP(ISYM,10) > 0)) then
    IDISK = IDSMAT(ISYM,10)
    call DDAFILE(LUSBT,1,SG,NSG,IDISK)
    if ((NINM > 0) .and. (NINDEP(ISYM,11) > 0)) then
      IDISK = IDSMAT(ISYM,11)
      call DDAFILE(LUSBT,1,SG,NSG,IDISK)
    end if
    call mma_deallocate(SG)
  end if
end do

end subroutine MKSG
