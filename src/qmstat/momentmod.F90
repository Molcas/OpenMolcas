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

subroutine MomentMod(Re,NRe,Cmo,nBRe,nBNRe,LindMOs,iS1,iS2,First,DiffMax)

use qmstat_global, only: iPrint
use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBRe, nBNRe, iS1, iS2
real(kind=wp), intent(in) :: Re(nTri_Elem(nBRe)), NRe(nBNRe,nBNRe), Cmo(nBNRe,nBNRe)
logical(kind=iwp), intent(in) :: LindMOs(nBNRe)
logical(kind=iwp), intent(inout) :: First
real(kind=wp), intent(out) :: DiffMax
integer(kind=iwp) :: i, icomp, iopt, irc, iSmLbl, j, kaunt1, nSize1, nSize2
character(len=8) :: Label
real(kind=wp) :: Diffx, Diffy, Diffz, DipRe(3), DipNRe(3)
real(kind=wp), allocatable :: D(:), Dsq(:,:), DxM(:,:), DxRe(:), DyM(:,:), DyRe(:), DzM(:,:), DzRe(:), TEMP(:,:)
real(kind=wp), external :: Ddot_

if (First .and. (iPrint >= 5)) then
  write(u6,*)
  write(u6,*) '     Modifications of dipoles by renormalization and basis reduction.'
  write(u6,*)
  write(u6,*) '     State pair    |  Difference '
  write(u6,*) '     --------------|---------------------'
  First = .false.
end if

nSize1 = nTri_Elem(nBNRe)
nSize2 = nTri_Elem(nBRe)
call mma_allocate(D,nSize1,label='Dip')
call mma_allocate(DxRe,nSize2,label='DipXre')
call mma_allocate(DyRe,nSize2,label='DipYre')
call mma_allocate(DzRe,nSize2,label='DipZre')
call mma_allocate(Dsq,nBNRe,nBNre,label='Dipsq')
call mma_allocate(DxM,nBNRe,nBNRe,label='DipXm')
call mma_allocate(DyM,nBNRe,nBNRe,label='DipYm')
call mma_allocate(DzM,nBNRe,nBNRe,label='DipZm')
call mma_allocate(TEMP,nBNRe,nBNRe,label='TEMP')
irc = -1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
iSmLbl = 0
Label = 'Mltpl  1'
! X
icomp = 1
call RdOne(irc,iopt,Label,icomp,D,iSmLbl)
call Square(D,Dsq,1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Cmo,nBNRe,Dsq,nBNRe,Zero,TEMP,nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,TEMP,nBNRe,Cmo,nBNRe,Zero,DxM,nBNRe)
! Y
icomp = 2
call RdOne(irc,iopt,Label,icomp,D,iSmLbl)
call Square(D,Dsq,1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Cmo,nBNRe,Dsq,nBNRe,Zero,TEMP,nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,TEMP,nBNRe,Cmo,nBNRe,Zero,DyM,nBNRe)
! Z
icomp = 3
call RdOne(irc,iopt,Label,icomp,D,iSmLbl)
call Square(D,Dsq,1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Cmo,nBNRe,Dsq,nBNRe,Zero,TEMP,nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,TEMP,nBNRe,Cmo,nBNRe,Zero,DzM,nBNRe)
! Triangularize and reduce.
kaunt1 = 0
do i=1,nBNRe
  if (.not. LindMOs(i)) cycle
  do j=1,i
    if (.not. LindMOs(j)) cycle
    kaunt1 = kaunt1+1
    DxRe(kaunt1) = DxM(j,i)
    DyRe(kaunt1) = DyM(j,i)
    DzRe(kaunt1) = DzM(j,i)
  end do
end do
! Density
DipNRe(1) = Ddot_(nBNRe**2,DxM,1,NRe,1)
DipNRe(2) = Ddot_(nBNRe**2,DyM,1,NRe,1)
DipNRe(3) = Ddot_(nBNRe**2,DzM,1,NRe,1)
DipRe(1) = Ddot_(nSize2,DxRe,1,Re,1)
DipRe(2) = Ddot_(nSize2,DyRe,1,Re,1)
DipRe(3) = Ddot_(nSize2,DzRe,1,Re,1)
Diffx = abs(DipRe(1)-DipNRe(1))
Diffy = abs(DipRe(2)-DipNRe(2))
Diffz = abs(DipRe(3)-DipNRe(3))
if (iPrint >= 5) write(u6,99) iS1,iS2,'(',Diffx,',',Diffy,',',Diffz,')'
! Return number
DiffMax = Diffy
if (Diffx >= Diffy) DiffMax = Diffx
if ((Diffz >= Diffx) .and. (Diffz >= Diffy)) DiffMax = Diffz
! Deallocate en masse.
call mma_deallocate(D)
call mma_deallocate(DxRe)
call mma_deallocate(DyRe)
call mma_deallocate(DzRe)
call mma_deallocate(Dsq)
call mma_deallocate(DxM)
call mma_deallocate(DyM)
call mma_deallocate(DzM)
call mma_deallocate(TEMP)

return

99 format('     ',2I3,'          ',3(A,F10.7),A)

end subroutine MomentMod
