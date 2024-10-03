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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RmLDep(AMat,lDm,lth)
!***********************************************************************
!                                                                      *
!     purpose: Remove linear dependencies                              *
!                                                                      *
!     input:                                                           *
!       AMat    : matrix with linear dependencies (lDm,lDm)            *
!                                                                      *
!     output:                                                          *
!       AMat    : inverted matrix without linear dependencies (lDm,lDm)*
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: lDm, lth
real(kind=wp), intent(inout) :: AMat(lDm,lDm)
real(kind=wp) :: Dummy
integer(kind=iwp) :: i, iDum, iErr, ij, lthS, lthT, nFound
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: j
#endif
real(kind=wp), allocatable :: ATri(:), EVal(:), EVec(:), Scr(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

lthT = nTri_Elem(lth)
lthS = lth**2
call mma_allocate(ATri,lthT,Label='ATri')
call mma_allocate(EVec,lthS,Label='EVec')
call mma_allocate(EVal,lth,Label='EVal')

! Put a unit matrix into the eigenvectors work space
call unitmat(EVec,lth)

! Copy triangular part of AMat to work space
ij = 1
do i=1,lth
  ATri(ij:ij+i-1) = AMat(i,1:i)
  ij = ij+i
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' Squared A-matrix in RmLDep:'
do i=1,lth
  write(u6,'(5(1x,es13.6))') (AMat(i,j),j=1,lth)
end do
write(u6,*) ' Triangular A-matrix:'
write(u6,'(5(1x,es13.6))') (ATri(i),i=1,lthT)
#endif

! Diagonalize
call mma_allocate(Scr,lth**2,Label='Scr')
Dummy = Zero
iDum = 0
call Diag_Driver('V','A','L',lth,ATri,Scr,lth,Dummy,Dummy,iDum,iDum,EVal,EVec,lth,1,0,'J',nFound,iErr)
call mma_deallocate(Scr)
#ifdef _DEBUGPRINT_
write(u6,*) ' Eigenvalues of A-matrix in RLnDep:'
write(u6,'(5(1x,es13.6))') (EVal(i),i=1,lth)
#endif

! Form the inverse
AMat(:,1:lth) = Zero
do i=1,lth
  if (EVal(i) > 1.0e-12_wp) then
    AMat(i,i) = One/EVal(i)
  else
    !write(u6,*) ' Eigenvalue',i,' smaller then 1.0e-12'
    AMat(i,i) = Zero
  end if
end do
call mma_allocate(Scr,lthS,Label='Scr')
call DGEMM_('N','T',lth,lth,lth, &
            One,AMat,lDm, &
            EVec,lth, &
            Zero,Scr,lth)
call DGEMM_('N','N',lth,lth,lth, &
            One,EVec,lth, &
            Scr,lth, &
            Zero,AMat,lDm)
call mma_deallocate(Scr)

call mma_deallocate(EVal)
call mma_deallocate(EVec)
call mma_deallocate(ATri)

end subroutine RmLDep
