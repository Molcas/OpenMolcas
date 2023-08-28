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

subroutine Rys2DN(xyz2D,xyz2DN,nArg,lRys,nabMax,ncdMax,CoorAC)
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
real(kind=wp), intent(inout) :: xyz2D(lRys,nArg,3,0:nabMax+2,0:ncdMax+2)
real(kind=wp), intent(out) ::   xyz2DN(lRys,nArg,3,2,0:nabMax+1,0:ncdMax+1)
real(kind=wp), intent(in) :: CoorAC(3,2)
integer(kind=iwp) :: ie, if, iRys

#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iOrdOp, iab, icd
character(LEN=30) Label
call RecPrt('Rys2Dn: AC',' ',CoorAC,3,2)
do iab = 0, nabMax+2
   do icd = 0, ncdMax+2
      write(Label,'(A,I3,A,I3,A)') ' In Rys2Dn: xyz2D(x)(',iab,',',icd,')'
      call RECPRT(Label,' ',xyz2D(:,:,1,iab,icd),lRys,nArg)
      write(Label,'(A,I3,A,I3,A)') ' In Rys2Dn: xyz2D(y)(',iab,',',icd,')'
      call RECPRT(Label,' ',xyz2D(:,:,2,iab,icd),lRys,nArg)
      write(Label,'(A,I3,A,I3,A)') ' In Rys2Dn: xyz2D(z)(',iab,',',icd,')'
      call RECPRT(Label,' ',xyz2D(:,:,3,iab,icd),lRys,nArg)
  end do
end do
#endif

xyz2DN(:,:,:,:,:,:) = Zero

do ie=0,nabMax+1
  do if=0,ncdMax+1
    do iRys = 1, lRys
      xyz2DN(iRys,:,1,1,ie,if) = xyz2D(iRys,:,1,ie+1,if) - xyz2D(iRys,:,1,ie,if+1) + (CoorAC(1,1)-CoorAC(1,2))*xyz2D(iRys,:,1,ie,if)
      xyz2DN(iRys,:,2,1,ie,if) = xyz2D(iRys,:,2,ie+1,if) - xyz2D(iRys,:,2,ie,if+1) + (CoorAC(2,1)-CoorAC(2,2))*xyz2D(iRys,:,2,ie,if)
      xyz2DN(iRys,:,3,1,ie,if) = xyz2D(iRys,:,3,ie+1,if) - xyz2D(iRys,:,3,ie,if+1) + (CoorAC(3,1)-CoorAC(3,2))*xyz2D(iRys,:,3,ie,if)
    end do
  end do
end do
do ie=0,nabMax
  do if=0,ncdMax
    do iRys = 1, lRys
      xyz2DN(iRys,:,1,2,ie,if) = xyz2DN(iRys,:,1,1,ie+1,if) - xyz2DN(iRys,:,1,1,ie,if+1) &
                               + (CoorAC(1,1)-CoorAC(1,2))*xyz2DN(iRys,:,1,1,ie,if)
      xyz2DN(iRys,:,2,2,ie,if) = xyz2DN(iRys,:,2,1,ie+1,if) - xyz2DN(iRys,:,2,1,ie,if+1) &
                               + (CoorAC(2,1)-CoorAC(2,2))*xyz2DN(iRys,:,2,1,ie,if)
      xyz2DN(iRys,:,3,2,ie,if) = xyz2DN(iRys,:,3,1,ie+1,if) - xyz2DN(iRys,:,3,1,ie,if+1) &
                               + (CoorAC(3,1)-CoorAC(3,2))*xyz2DN(iRys,:,3,1,ie,if)
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Generalized 2D integrals: 2D(:,ie,if,iOrdOp)'
do iOrdOp=1, 2
do ie=0,nabMax
  do if=0,ncdMax
    write(Label,'(A,I2,A,I2,A,I2,A)') 'Rys2Dn: 2DN(',ie,',',if,',',iOrdOp,')(x)'
    call RecPrt(Label,' ',xyz2DN(:,:,1,iOrdOp,ie,if),lRys,nArg)
    write(Label,'(A,I2,A,I2,A,I2,A)') 'Rys2Dn: 2DN(',ie,',',if,',',iOrdOp,')(y)'
    call RecPrt(Label,' ',xyz2DN(:,:,2,iOrdOp,ie,if),lRys,nArg)
    write(Label,'(A,I2,A,I2,A,I2,A)') 'Rys2Dn: 2DN(',ie,',',if,',',iOrdOp,')(z)'
    call RecPrt(Label,' ',xyz2DN(:,:,3,iOrdOp,ie,if),lRys,nArg)
  end do
end do
end do
#endif

return

end subroutine Rys2DN
