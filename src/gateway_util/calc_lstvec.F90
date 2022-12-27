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
! Copyright (C) 2009, Roland Lindh                                     *
!               2010, Mickael G. Delcey                                *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine calc_LSTvec(mynRP,Reac,Prod,TanVec,Invar)

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mynRP
real(kind=wp), intent(in) :: Reac(mynRP), Prod(mynRP)
real(kind=wp), intent(out) :: TanVec(mynRP)
logical(kind=iwp), intent(in) :: Invar
integer(kind=iwp) :: i, iAt, iCnt, iProdA, iReacA, j, jTmp, mAt, nAt, nData, nsc
real(kind=wp) :: norm, RMax, RMSD
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: iStab(:)
real(kind=wp), allocatable :: W(:), XYZ(:,:)
real(kind=wp), external :: DDot_

!***********************************************************************
!                                                                      *
!      Object : compute the LST vector, that is Reac-Prod              *
!      Of course it is not that simple, because of metric problems     *
!                                                                      *
!      M.G. Delcey     June 2010                                       *
!      Lund University                                                 *
!                                                                      *
!***********************************************************************
!
! Get the symmetry stabilizers for each center

nAt = mynRP/3
call mma_Allocate(iStab,nAt,label='iStab')
iAt = 1
nsc = 0
do i=1,nCnttp
  do iCnt=1,dbsc(i)%nCntr
    nsc = nsc+1
    if (.not. (dbsc(i)%pChrg .or. dbsc(i)%Frag .or. dbsc(i)%Aux)) then
      jTmp = 0
      do j=1,dc(nsc)%nStab-1
        jTmp = ior(jTmp,dc(nsc)%iStab(j))
      end do
      iStab(iAt) = jTmp
      iAt = iAt+1
    end if
  end do
end do

! Align the structures and compute the vector (un-weighted)

call mma_allocate(XYZ,3*nAt*8,2)
iReacA = 1
iProdA = 2
call Expand_Coor(Reac,nAt,XYZ(1,iReacA),mAt)
call Expand_Coor(Prod,nAt,XYZ(1,iProdA),mAt)
call Qpg_dArray('Weights',Found,nData)
if (Found .and. (nData >= mAt)) then
  call mma_allocate(W,nData,label='W')
  call Get_dArray('Weights',W,nData)
else
  call SysAbendMsg('calc_LSTvec','No or wrong weights were found in the RUNFILE.','')
end if
if (Invar) then
  call Superpose_w(XYZ(1,iReacA),XYZ(1,iProdA),W,mAt,RMSD,RMax)
  call Fix_Symmetry(XYZ(1,iReacA),nAt,iStab)
end if
TanVec(:) = XYZ(1:mynRP,iReacA)-XYZ(1:mynRP,iProdA)
call mma_deallocate(XYZ)
call mma_deallocate(iStab)
call mma_deallocate(W)

! And normalize it

norm = DDot_(mynRP,TanVec,1,TanVec,1)
TanVec(:) = TanVec/sqrt(norm)
!call RecPrt('TanVec',' ',TanVec,3,nAt)

end subroutine calc_LSTvec
