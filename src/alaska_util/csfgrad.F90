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

implicit none
#include "stdalloc.fh"
#include "real.fh"
#include "nac.fh"
integer nGrad
real*8 Grad(nGrad)
integer nD, lOper(1)
real*8 CCoor(3)
real*8, dimension(:), allocatable :: aDAO
logical Found
character(len=80) Label
external OvrGrd, OvrMmG

call DCopy_(nGrad,[Zero],0,Grad,1)

!nB = nBas(0)
call Qpg_dArray('D1ao-',Found,nD)
call mma_allocate(aDAO,nD)
call Get_dArray('D1ao-',aDAO,nD)
!call TriPrt('DAO-','',aDAO,nB)

!IFG Compute the CSF contribution to the coupling vector.
!    Inner product of S[x] and D^A (antisymmetric component of transition density matrix)
!    This is the same as the product of S[x]^A and D
isCSF = .true.
call DCopy_(3,[Zero],0,CCoor,1)
lOper(1) = 1
Label = 'The CSF Contribution'
call OneEl_g(OvrGrd,OvrMmG,Grad,nGrad,.false.,CCoor,aDAO,nD,lOper,1,0,Label)
isCSF = .false.

call mma_deallocate(aDAO)

end subroutine CSFGRad
