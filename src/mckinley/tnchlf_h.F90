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
! Copyright (C) 1990,1994,1996, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine Tnchlf_h(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,lZeta,nVec,IncVec,A1,A2,A3,Indij)
!***********************************************************************
!                                                                      *
! Object: to do a half transformation. The loop over the two matrix-   *
!         matrix multiplications is segmented such that the end of the *
!         intermediate matrix will not push the start of the same out  *
!         from the cache.                                              *
!                                                                      *
! Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!                                                                      *
!             Modified to decontraction May 1996, by R. Lindh          *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCntr1, nPrm1, nCntr2, nPrm2, lZeta, nVec, IncVec, Indij(lZeta)
real(kind=wp), intent(in) :: Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2), A1(nCntr1,nCntr2,nVec)
real(kind=wp), intent(out) :: A2(nCntr2,IncVec*nPrm1), A3(nVec,lZeta)
integer(kind=iwp) :: iCntr1, iCntr2, iiVec, ijVec, iPrm1, iPrm2, iZeta, mVec
logical(kind=iwp) :: Seg1, Seg2

! Check if the basis set is segmented

Seg1 = .false.
loop1: do iPrm1=nPrm1,1,-1
  do iCntr1=nCntr1,1,-1
    if (Coeff1(iPrm1,iCntr1) == Zero) then
      Seg1 = .true.
      exit loop1
    end if
  end do
end do loop1

Seg2 = .false.
loop2: do iPrm2=nPrm2,1,-1
  do iCntr2=nCntr2,1,-1
    if (Coeff2(iPrm2,iCntr2) == Zero) then
      Seg2 = .true.
      exit loop2
    end if
  end do
end do loop2

! Set output matrix to zero

A3(:,:) = Zero

! Loop sectioning

do iiVec=1,nVec,IncVec
  mVec = min(IncVec,nVec-iiVec+1)
  ! Set intermediate matrix to zero
  A2(:,1:mVec*nPrm1) = Zero

  if (Seg1) then

    ! First quarter transformation

    do iPrm1=1,nPrm1
      ijVec = mVec*(iPrm1-1)+1
      do iCntr1=1,nCntr1
        ! Check for zero due to segmented basis
        if (abs(Coeff1(iPrm1,iCntr1)) > Zero) then
          do iCntr2=1,nCntr2
            A2(iCntr2,ijVec:ijVec+mVec-1) = A2(iCntr2,ijVec:ijVec+mVec-1)+Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iiVec:iiVec+mVec-1)
          end do
        end if
      end do
    end do

  else  ! Seg1

    ! First quarter transformation

    do iPrm1=1,nPrm1
      ijVec = mVec*(iPrm1-1)+1
      do iCntr1=1,nCntr1
        do iCntr2=1,nCntr2
          A2(iCntr2,ijVec:ijVec+mVec-1) = A2(iCntr2,ijVec:ijVec+mVec-1)+Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iiVec:iiVec+mVec-1)
        end do
      end do
    end do

  end if  ! Seg1

  if (Seg2) then

    ! Second quarter transformation

    do iCntr2=1,nCntr2
      do iZeta=1,lZeta
        iPrm2 = (Indij(iZeta)-1)/nPrm1+1
        iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
        ! Check for zero due to segmented basis
        if (abs(Coeff2(iPrm2,iCntr2)) > Zero) then
          ijVec = mVec*(iPrm1-1)+1
          A3(iiVec:iiVec+mVec-1,iZeta) = A3(iiVec:iiVec+mVec-1,iZeta)+Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec:ijVec+mVec-1)
        end if
      end do
    end do

  else

    ! Second quarter transformation

    do iCntr2=1,nCntr2
      do iZeta=1,lZeta
        iPrm2 = (Indij(iZeta)-1)/nPrm1+1
        iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
        ijVec = mVec*(iPrm1-1)+1
        A3(iiVec:iiVec+mVec-1,iZeta) = A3(iiVec:iiVec+mVec-1,iZeta)+Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec:ijVec+mVec-1)
      end do
    end do

  end if

  ! End of loop sectioning

end do

return

end subroutine Tnchlf_h
