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

subroutine RICD_Helper(Do_nacCD_Basis,nTest,iAngMin_,iAngMax_,jAngMin_,jAngMax_,nBS,iAng,jAng,list,nBS_Max)

implicit real*8(a-h,o-z)
logical Do_nacCD_Basis
integer iAngMin_(0:nBS_Max-1), iAngMax_(0:nBS_Max-1)
parameter(iTabMx=15)
integer jAngMin_(0:nBS_Max-1,0:nBS_Max-1), jAngMax_(0:nBS_Max-1,0:nBS_Max-1)
integer list(2,0:((nTest+1)*(nTest+2))/2,0:nTest*2)
integer list2(0:nTest**2)

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Do_nacCD_Basis) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nBS = (nTest+2)/2
  do iBS=0,nBS-1
    iAngMin_(iBS) = iBS
    iAngMax_(iBS) = nTest-iBS
    do iAng=0,iAngMax_(iBS)
      jAngMax_(iBS,iAng) = min(iAng,iAngMin_(iBS))
      if (iAng == iAngMax_(iBS)) jAngMax_(iBS,iAng) = iAngMax_(iBS)
      if (iAng < iAngMin_(iBS)) jAngMax_(iBS,iAng) = 0
      jAngMin_(iBS,iAng) = iAngMin_(iBS)
      if (iAng <= iAngMin_(iBS)) jAngMin_(iBS,iAng) = 0
      do jAng=jAngMin_(iBS,iAng),jAngMax_(iBS,iAng)
        list(1,0,iAng) = iAng
        list(2,0,iAng) = jAng
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  nBS = 1
  iPair = 0
  do iBS=0,nBS-1
    iAngMax_(iBS) = nTest*2
    do iAng=iAngMin_(iBS),iAngMax_(iBS)
      jAngMax_(iBS,iAng) = 0
      jAngMin_(iBS,iAng) = 0
      do jAng=jAngMin_(iBS,iAng),jAngMax_(iBS,iAng)
        list2(iAng) = 0
        do k=0,nTest
          do l=0,k
            do m=iAng,0,-2
              n = k-l
              if ((n == m) .and. (k+l >= iAng)) then
                iPair = list2(iAng)
                list(1,iPair,iAng) = l
                list(2,iPair,iAng) = k
                list2(iAng) = list2(iAng)+1
              end if
            end do ! m
          end do   ! l
        end do     ! k
      end do       ! jAng
    end do         ! iAng
  end do           ! iBS
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine RICD_Helper
