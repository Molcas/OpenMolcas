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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1991, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MkSrt0(iSquar,nIrrep,nBas,iSkip)
!***********************************************************************
!                                                                      *
!     Purpose: Set up all information needed to compute 2el integral   *
!              addresses in the subroutine PLF2 and INDSFT2.           *
!                                                                      *
!     Called from: Sort0                                               *
!                                                                      *
!     Calls to : none                                                  *
!                                                                      *
!     Calling parameters:                                              *
!     iSquar  : ordering flag                                          *
!     nIrrep  : number of irreducible representations                  *
!     nBas    : number of basis functions per irred. rep.              *
!     iSkip   : flag to exlude symmetry combinations                   *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use sort_data, only: DimSyB, TriSyB, mxSyP, nBs, nSkip, nSyOp, Square
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSquar, nIrrep, nBas(nIrrep), iSkip(nIrrep)
#include "print.fh"
integer(kind=iwp) :: iBsi, iPrint, iRout, iSymi, jBsj, jSymj

iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 10) write(u6,*) ' >>> Enter MKSRT0 <<<'
!----------------------------------------------------------------------*
!     Gather information on desired ordering scheme                    *
!----------------------------------------------------------------------*

Square = .true.
if (iSquar == 0) Square = .false.

!----------------------------------------------------------------------*
!     Gather data on the number of symmetry operations                 *
!----------------------------------------------------------------------*

nSyOp = nIrrep
mxSyP = nSyOp*(nSyOp+1)/2

!----------------------------------------------------------------------*
!     Gather data on the number of basis functions                     *
!----------------------------------------------------------------------*

do iSymi=1,nSyOp
  nBs(iSymi) = nBas(iSymi)
end do

!----------------------------------------------------------------------*
!     Put flags to exclude symmetry combinations into common block     *
!----------------------------------------------------------------------*

do iSymi=1,nSyOp
  nSkip(iSymi) = iSkip(iSymi)
end do

!----------------------------------------------------------------------*
!     Precompute the dimension of the symmetry blocks                  *
!     and symmetry block numbers for pairs of symmtry indices          *
!----------------------------------------------------------------------*

do iSymi=1,nSyOp
  iBsi = nBs(iSymi)
  DimSyB(iSymi,iSymi) = iBsi*(iBsi+1)/2
  TriSyB(iSymi,iSymi) = iSymi*(iSymi+1)/2
  do jSymj=1,iSymi-1
    jBsj = nBs(jSymj)
    DimSyB(iSymi,jSymj) = iBsi*jBsj
    DimSyB(jSymj,iSymi) = jBsj*iBsi
    TriSyB(iSymi,jSymj) = jSymj+iSymi*(iSymi-1)/2
    TriSyB(jSymj,iSymi) = jSymj+iSymi*(iSymi-1)/2
  end do
end do

return

end subroutine MkSrt0
