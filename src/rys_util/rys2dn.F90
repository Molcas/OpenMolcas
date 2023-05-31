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
! Copyright (C) 1990,1991,2023, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
subroutine Rys2DN(xyz2D,xyz2DN,nArg,lRys,nabMax,ncdMax,P,Q,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the 2-dimensional integrals of the Rys            *
!         quadrature. The z components are assumed to be pre-          *
!         conditioned with the weights of the roots of the             *
!         Rys polynomial.                                              *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
! Modified loop structure for RISC 1991 R. Lindh Dept. of Theoretical  *
! Chemistry, University of Lund, Sweden.                               *
! VV: improve loop structure                                           *
! Modified for use to compute generalized 2D integrals to compute      *
! integrals for the Breit-Dirac and the Berit-Pauli Hamiltonian.       *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nArg, lRys, nabMax, ncdMax
real(kind=wp), intent(inout) :: xyz2D(nArg*lRys,3,0:nabMax+nOrdOp,0:ncdMax+nOrdOp)
real(kind=wp), intent(out) ::   xyz2DN(nArg*lRys,3,nOrdOp,0:nabMax,0:ncdMax)
real(kind=wp), intent(in) :: P(3), Q(3)
integer(kind=iwp) :: ie, if

#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iOrdOp
character(LEN=30) Label
if (nabMax > 0) call RecPrt('P',' ',P,1,3)
if (ncdMax > 0) call RecPrt('Q',' ',Q,1,3)
#endif

xyz2DN(:,:,:,:) = Zero

do ie=0,nabMax
  do if=0,ncdMax
    xyz2DN(:,1,1,ie,if) = xyz2D(:,1,ie+1,if) - xyz2D(:,1,ie,if+1) + (P(1)-Q(1))*xyz2D(:,1,ie,if)
    xyz2DN(:,2,1,ie,if) = xyz2D(:,2,ie+1,if) - xyz2D(:,2,ie,if+1) + (P(2)-Q(2))*xyz2D(:,2,ie,if)
    xyz2DN(:,3,1,ie,if) = xyz2D(:,3,ie+1,if) - xyz2D(:,3,ie,if+1) + (P(3)-Q(3))*xyz2D(:,3,ie,if)
  end do
end do

If (nOrdOp==2) Then
do ie=0,nabMax
  do if=0,ncdMax
    xyz2DN(:,1,2,ie,if) = xyz2DN(:,1,2,ie+1,if) - xyz2DN(:,1,2,ie,if+1) + (P(1)-Q(1))*xyz2DN(:,1,2,ie,if)
    xyz2DN(:,2,2,ie,if) = xyz2DN(:,2,2,ie+1,if) - xyz2DN(:,2,2,ie,if+1) + (P(2)-Q(2))*xyz2DN(:,2,2,ie,if)
    xyz2DN(:,3,2,ie,if) = xyz2DN(:,3,2,ie+1,if) - xyz2DN(:,3,2,ie,if+1) + (P(3)-Q(3))*xyz2DN(:,3,2,ie,if)
  end do
end do
End If

#ifdef _DEBUGPRINT_
write(u6,*) ' Generalized 2D integrals: 2D(:,ie,if,iOrdOp)'
do iOrdOp=1, nOrdOp
do ie=0,nabMax
  do if=0,ncdMax
    write(Label,'(A,I2,A,I2,A)') ' 2DN(',ie,',',if,',',iOrdOp')(x)'
    call RecPrt(Label,' ',xyz2DN(:,1,iOrdOp,ie,if),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') ' 2DN(',ie,',',if,',',iOrdOp')(y)'
    call RecPrt(Label,' ',xyz2DN(:,2,iOrdOp,ie,if),lRys,nArg)
    write(Label,'(A,I2,A,I2,A)') ' 2DN(',ie,',',if,',',iOrdOp')(z)'
    call RecPrt(Label,' ',xyz2DN(:,3,iOrdOp,ie,if),lRys,nArg)
  end do
end do
end do
#endif

return

end subroutine Rys2DN
