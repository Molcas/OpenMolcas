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

subroutine Tnchlf_h(Coeff1,nCntr1,nPrm1,Coeff2,nCntr2,nPrm2,mZeta,lZeta,nVec,IncVec,A1,A2,A3,Indij)
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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 Coeff1(nPrm1,nCntr1), Coeff2(nPrm2,nCntr2), A1(nCntr1,nCntr2,nVec), A2(nCntr2,IncVec*nPrm1), A3(nVec,lZeta)
integer Indij(lZeta)
logical Seg1, Seg2

! Check if the basis set is segmented

Seg1 = .false.
do iPrm1=nPrm1,1,-1
  do iCntr1=nCntr1,1,-1
    if (Coeff1(iPrm1,iCntr1) == Zero) then
      Seg1 = .true.
      Go To 10
    end if
  end do
end do
10 continue

Seg2 = .false.
do iPrm2=nPrm2,1,-1
  do iCntr2=nCntr2,1,-1
    if (Coeff2(iPrm2,iCntr2) == Zero) then
      Seg2 = .true.
      Go To 20
    end if
  end do
end do
20 continue

! Set output matrix to zero

call FZero(A3,nVec*lZeta)

! Loop sectioning

do iiVec=1,nVec,IncVec
  mVec = min(IncVec,nVec-iiVec+1)
  ! Set intermediate matrix to zero
  call dcopy_(nCntr2*mVec*nPrm1,[Zero],0,A2,1)

  if (Seg1) then

    ! First quarter transformation

    do iPrm1=1,nPrm1
      do iCntr1=1,nCntr1
        ! Check for zero due to segmented basis
        if (abs(Coeff1(iPrm1,iCntr1)) > Zero) then
          do iCntr2=1,nCntr2
            do iVec=iiVec,iiVec+mVec-1
              ijVec = mVec*(iPrm1-1)+(iVec-iiVec+1)
              A2(iCntr2,ijVec) = A2(iCntr2,ijVec)+Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iVec)
            end do
          end do
        end if
      end do
    end do

  else  ! Seg1

    ! First quarter transformation

    do iPrm1=1,nPrm1
      do iCntr1=1,nCntr1
        do iCntr2=1,nCntr2
          do iVec=iiVec,iiVec+mVec-1
            ijVec = mVec*(iPrm1-1)+(iVec-iiVec+1)
            A2(iCntr2,ijVec) = A2(iCntr2,ijVec)+Coeff1(iPrm1,iCntr1)*A1(iCntr1,iCntr2,iVec)
          end do
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
          do iVec=iiVec,iiVec+mVec-1
            ijVec = mVec*(iPrm1-1)+(iVec-iiVec+1)
            A3(iVec,iZeta) = A3(iVec,iZeta)+Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec)
          end do
        end if
      end do
    end do

  else

    ! Second quarter transformation

    do iCntr2=1,nCntr2
      do iZeta=1,lZeta
        iPrm2 = (Indij(iZeta)-1)/nPrm1+1
        iPrm1 = Indij(iZeta)-(iPrm2-1)*nPrm1
        do iVec=iiVec,iiVec+mVec-1
          ijVec = mVec*(iPrm1-1)+(iVec-iiVec+1)
          A3(iVec,iZeta) = A3(iVec,iZeta)+Coeff2(iPrm2,iCntr2)*A2(iCntr2,ijVec)
        end do
      end do
    end do

  end if

  ! End of loop sectioning

end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(mZeta)

end subroutine Tnchlf_h
