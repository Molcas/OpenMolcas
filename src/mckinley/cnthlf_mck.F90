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
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine Cnthlf_mck(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,nZeta,lZeta,nVec,First,IncVec,A1,A2,A3,Indij)
!***********************************************************************
!                                                                      *
! Object: to do a half transformation. The loop over the two matrix-   *
!         matrix multiplications is segmented such that the end of the *
!         intermediate matrix will not push the start of the same out  *
!         from the cache.                                              *
!                                                                      *
! Author:     Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCntr1, nPrm1, nCntr2, nPrm2, nZeta, lZeta, nVec, IncVec, Indij(nZeta)
real(kind=wp), intent(in) :: Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2), A1(lZeta,nVec)
logical(kind=iwp), intent(in) :: First
real(kind=wp), intent(out) :: A2(nPrm2,IncVec*nCntr1)
real(kind=wp), intent(inout) :: A3(nVec,nCntr1,nCntr2)
integer(kind=iwp) :: iCntr1, iCntr2, iiVec, ijVec, iPrm1, iPrm2, iZeta, jZeta, mVec
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

if (First) A3(:,:,:) = Zero

! Loop sectioning

do iiVec=1,nVec,IncVec
  mVec = min(IncVec,nVec-iiVec+1)
  ! Set intermediate matrix to zero
  A2(:,1:nCntr1*mVec) = Zero

  if (Seg1) then

    ! First quarter transformation

    do iPrm1=1,nPrm1
      do iCntr1=1,nCntr1
        ijVec = mVec*(iCntr1-1)+1
        ! Check for zero due to segmented basis
        if (abs(Coeff1(iPrm1,iCntr1)) > Zero) then
          do iPrm2=1,nPrm2
            iZeta = (iPrm2-1)*nPrm1+iPrm1
            jZeta = Indij(iZeta)
            ! Skip due to screening
            if (jZeta > 0) then
              A2(iPrm2,ijVec:ijVec+mVec-1) = A2(iPrm2,ijVec:ijVec+mVec-1)+Coeff1(iPrm1,iCntr1)*A1(jZeta,iiVec:iiVec+mVec-1)
            end if
          end do
        end if
      end do
    end do

  else

    ! First quarter transformation

    do iPrm1=1,nPrm1
      do iCntr1=1,nCntr1
        ijVec = mVec*(iCntr1-1)+1
        do iPrm2=1,nPrm2
          iZeta = (iPrm2-1)*nPrm1+iPrm1
          jZeta = Indij(iZeta)
          ! Skip due to screening
          if (jZeta > 0) then
            A2(iPrm2,ijVec:ijVec+mVec-1) = A2(iPrm2,ijVec:ijVec+mVec-1)+Coeff1(iPrm1,iCntr1)*A1(jZeta,iiVec:iiVec+mVec-1)
          end if
        end do
      end do
    end do

  end if

  if (Seg2) then

    ! Second quarter transformation

    do iPrm2=1,nPrm2
      do iCntr2=1,nCntr2
        ! Check for zero due to segmented basis
        if (abs(Coeff2(iPrm2,iCntr2)) > Zero) then
          do iCntr1=1,nCntr1
            ijVec = mVec*(iCntr1-1)+1
            A3(iiVec:iiVec+mVec-1,iCntr1,iCntr2) = A3(iiVec:iiVec+mVec-1,iCntr1,iCntr2)+ &
                                                   Coeff2(iPrm2,iCntr2)*A2(iPrm2,ijVec:ijVec+mVec-1)
          end do
        end if
      end do
    end do

  else

    ! Second quarter transformation

    do iPrm2=1,nPrm2
      do iCntr2=1,nCntr2
        do iCntr1=1,nCntr1
          ijVec = mVec*(iCntr1-1)+1
          A3(iiVec:iiVec+mVec-1,iCntr1,iCntr2) = A3(iiVec:iiVec+mVec-1,iCntr1,iCntr2)+ &
                                                 Coeff2(iPrm2,iCntr2)*A2(iPrm2,ijVec:ijVec+mVec-1)
        end do
      end do
    end do

  end if

  ! End of loop sectioning

end do

return

end subroutine Cnthlf_mck
