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
! Copyright (C) 1992, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************

#include "compiler_features.h"

!#define _DEBUGPRINT_
subroutine vOO2OV_inner(v1,n1,v2,n2,iD)
!***********************************************************************
!                                                                      *
!     purpose: converts vector of dim nOO (e.g. gradient) to vector    *
!              of lower dimension nOV (occ-virt block only)            *
!              if n1==nOO && n2==nOV  -> compress                      *
!              else if if n1==nOV && n2==nOO -> decompress             *
!                                                                      *
!     output:                                                          *
!       v2      : compressed or decompressed vector, dep. on input     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. Schuetz                                                       *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!     Roland Lindh                                                     *
!     Uppsala University, Sweden, 2017                                 *
!     Fortran pointer structure                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

#ifndef POINTER_REMAP
use, intrinsic :: iso_c_binding
#endif
use InfSCF, only: nOO, nSym, nFro, nOcc, nOrb, kOV

implicit none
! declaration subroutine parameters
integer n1, n2, iD
real*8, target :: v1(n1), v2(n2)
real*8, dimension(:,:), pointer :: pv1, pv2
! declaration local variables
integer iSym, ii, ia, ioffs, ivoffs
integer nia1, nia2, nii1, nii2, nv1, nv2
#ifdef _DEBUGPRINT_
integer iOff, nO, nV
#endif

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
write(6,*) 'n1,n2,nOO,kOV(:)=',n1,n2,nOO,kOV(:)
iOff = 1
do iSym=1,nSym
  if (n1 == nOO) then
    nO = nOrb(iSym)
    call RecPrt('vOO2OV: v1',' ',v1(ioff),nO,nO)
    iOff = iOff+nO**2
  else if (n1 == kOV(iD)) then
    nO = nOcc(iSym,iD)-nFro(iSym)
    nV = nOrb(iSym)-nOcc(iSym,iD)-nFro(iSym)
    call RecPrt('vOO2OV: v1',' ',v1(ioff),nO,nV)
    iOff = iOff+nO*nV
  end if
end do
#endif
ioffs = 0
ivoffs = 1
do iSym=1,nSym
  ! range for occupied non-frozen orbitals
  nii1 = nFro(iSym)+1
  nii2 = nOcc(iSym,iD)
  ! range for virtual orbitals
  nia1 = nOcc(iSym,iD)+1
  nia2 = nOrb(iSym)
  ! size of the full block
  nv1 = nia2**2
  ! size of the virtual-occupied block
  nv2 = (nia2-nia1+1)*(nii2-nii1+1)

  if ((n1 == nOO) .and. (n2 == kOV(iD))) then

    ! compress

#   ifdef POINTER_REMAP
    pv1(1:nia2,1:nia2) => v1(ioffs+1:ioffs+nv1)
    pv2(nia1:nia2,nii1:nii2) => v2(ivoffs:ivoffs+nv2-1)

    do ii=nii1,nii2
      do ia=nia1,nia2
        !if (pv1(ia,ii) /= -pv1(ii,ia)) then
        !  write(6,*) 'inconsistency in gradient'
        !  call Abend()
        !end if
        pv2(ia,ii) = pv1(ia,ii)
      end do
    end do
#   else
    call c_f_pointer(c_loc(v1(ioffs+1)),pv1,[nia2,nia2])
    call c_f_pointer(c_loc(v2(ivoffs)),pv2,[nia2-nia1+1,nii2-nii1+1])

    do ii=nii1,nii2
      do ia=nia1,nia2
        !If (pv1(ia,ii) /= -pv1(ii,ia)) then
        !  write(6,*) 'inconsistency in gradient'
        !  call Abend()
        !end if
        pv2(ia-nia1+1,ii-nii1+1) = pv1(ia,ii)
      end do
    end do
#   endif

  else if ((n1 == kOV(iD)) .and. (n2 == nOO)) then

    ! decompress

#   ifdef POINTER_REMAP
    pv1(nia1:nia2,nii1:nii2) => v1(ivoffs:ivoffs+nv2-1)
    pv2(1:nia2,1:nia2) => v2(ioffs+1:ioffs+nv1)

    do ii=nii1,nii2
      do ia=nia1,nia2
        pv2(ia,ii) = pv1(ia,ii)
        pv2(ii,ia) = -pv1(ia,ii)
      end do
    end do
#   else
    call c_f_pointer(c_loc(v1(ivoffs)),pv1,[nia2-nia1+1,nii2-nii1+1])
    call c_f_pointer(c_loc(v2(ioffs+1)),pv2,[nia2,nia2])

    do ii=nii1,nii2
      do ia=nia1,nia2
        pv2(ia,ii) = pv1(ia-nia1+1,ii-nii1+1)
        pv2(ii,ia) = -pv1(ia-nia1+1,ii-nii1+1)
      end do
    end do
#   endif
  end if

  nullify(pv1,pv2)
  ivoffs = ivoffs+nv2
  ioffs = ioffs+(nOrb(iSym)*nOrb(iSym))
end do

#ifdef _DEBUGPRINT_
write(6,*) 'n1,n2,nOO,kOV=',n1,n2,nOO,kOV(:)
iOff = 1
do iSym=1,nSym
  if (n2 == nOO) then
    nO = nOrb(iSym)
    call RecPrt('vOO2OV: v2',' ',v2(ioff),nO,nO)
    iOff = iOff+nO**2
  else if (n2 == kOV(iD)) then
    nO = nOcc(iSym,iD)-nFro(iSym)
    nV = nOrb(iSym)-nOcc(iSym,iD)-nFro(iSym)
    call RecPrt('vOO2OV: v2',' ',v2(ioff),nO,nV)
    iOff = iOff+nO*nV
  end if
end do
#endif

return

end subroutine vOO2OV_inner
