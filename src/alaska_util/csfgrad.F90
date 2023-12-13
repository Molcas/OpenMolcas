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

subroutine CSFGRad(Grad,nGrad)
!***********************************************************************
!                                                                      *
! Object: to compute the CSF component of the non-adiabatic derivative *
!         coupling                                                     *
!                                                                      *
! This routine assumes C1 symmetry and needs the 'D1ao-' matrix stored *
! in the runfile                                                       *
!                                                                      *
!***********************************************************************

!use Basis_Info, only: nBas
use Grd_interface, only: grd_kernel, grd_mem
use NAC, only: IsCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: nD, lOper(1)
real(kind=wp) :: CCoor(3)
logical(kind=iwp) :: Found
character(len=80) :: Label
real(kind=wp), allocatable :: aDAO(:)
procedure(grd_kernel) :: OvrGrd
procedure(grd_mem) :: OvrMmG

Grad(:) = Zero

!nB = nBas(0)
call Qpg_dArray('D1ao-',Found,nD)
call mma_allocate(aDAO,nD)
call Get_dArray('D1ao-',aDAO,nD)
!call TriPrt('DAO-','',aDAO,nB)

!IFG Compute the CSF contribution to the coupling vector.
!    Inner product of S[x] and D^A (antisymmetric component of transition density matrix)
!    This is the same as the product of S[x]^A and D
isCSF = .true.
CCoor(:) = Zero
lOper(1) = 1
Label = 'The CSF Contribution'
call OneEl_g(OvrGrd,OvrMmG,Grad,nGrad,.false.,CCoor,aDAO,nD,lOper,1,0,Label)
isCSF = .false.

call mma_deallocate(aDAO)

end subroutine CSFGRad
