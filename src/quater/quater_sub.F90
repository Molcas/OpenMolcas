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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine quater_sub(nAtoms,G1,G2,ireturn)

use Quater_globals, only: debug, rotate, translate, ngeoms, list
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: G1(3*nAtoms)
real(kind=wp), intent(inout) :: G2(3*nAtoms)
integer(kind=iwp), intent(out) :: ireturn
real(kind=wp) :: U1(3), U2(3), V1(3), V2(3), V1best(3), V2best(3), Q(0:3), Vtrans(3)
integer(kind=iwp), parameter :: XYZ(3) = [1, 2, 3]

#ifdef _DEBUGPRINT_
debug = .true.
#else
debug = .false.
#endif
rotate = .true.
translate = .true.
ngeoms = 1
list(1)%title = 'ORIGINAL'
list(2)%title = 'GEOHYPER'
list(3)%title = 'bALIGNED'

call quaterinit()

list(1)%nat = nAtoms
list(2)%nat = nAtoms
list(3)%nat = nAtoms

call mma_allocate(list(1)%geo,3,nAtoms,label='GEO0')
call mma_allocate(list(1)%geolbl,nAtoms,label='GEO0')
call mma_allocate(list(2)%geo,3,nAtoms,label='GEO1')
call mma_allocate(list(2)%geolbl,nAtoms,label='GEO1')
call mma_allocate(list(3)%geo,3,nAtoms,label='GEO2')
call mma_allocate(list(3)%geolbl,nAtoms,label='GEO2')

call dcopy_(3*nAtoms,G1,1,list(1)%geo,1)
call dcopy_(3*nAtoms,G2,1,list(2)%geo,1)

call SetVect(nAtoms,list(1)%geo,XYZ,V1,V2)
call SetVect(nAtoms,list(2)%geo,XYZ,U1,U2)

if (debug) then
  write(u6,*) 'Reference axis'
  call RecPrt('U1',' ',U1,3,1)
  call RecPrt('U2',' ',U2,3,1)
  write(u6,*) 'New axis'
  call RecPrt('V1',' ',V1,3,1)
  call RecPrt('V2',' ',V2,3,1)
end if

call QuaterSolve(U1,U2,V1,V2,Q)

if (debug) then
  write(u6,*) 'Normalized Reference axis'
  call RecPrt('U1',' ',U1,3,1)
  call RecPrt('U2',' ',U2,3,1)
  write(u6,*) 'Normalized New axis'
  call RecPrt('V1',' ',V1,3,1)
  call RecPrt('V2',' ',V2,3,1)
  call QuaterRotation(Q,U1,V1best)
  call QuaterRotation(Q,U2,V2best)
  call RecPrt('Best V1',' ',V1best,3,1)
  call RecPrt('Best V2',' ',V2best,3,1)
end if

list(3)%geo(:,:) = list(2)%geo(:,:)

call RotateGeoms(Q)
call SetVectTrans(list(1)%nat,list(1)%geo,XYZ,list(3)%nat,list(3)%geo,XYZ,Vtrans)
call TranslateGeoms(Vtrans)

if (debug) then
  call PrintGeom(-1_iwp,list(3)%nat,list(3)%title,list(3)%geo,list(3)%geolbl)
end if

call dcopy_(3*nAtoms,list(3)%geo,1,G2,1)

call mma_deallocate(list(1)%geo)
call mma_deallocate(list(1)%geolbl)
call mma_deallocate(list(2)%geo)
call mma_deallocate(list(2)%geolbl)
call mma_deallocate(list(3)%geo)
call mma_deallocate(list(3)%geolbl)

ireturn = 0

return

end subroutine quater_sub
