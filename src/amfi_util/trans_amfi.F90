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

subroutine trans_amfi(coeffs,idim1,idim2,ich,nolds1,nolds2,nolds3,nolds4,nnew1,nnew2,nnew3,nnew4,array1,array2)
!bs makes the transformation for the ich-th index
!coeffs                      : (nolds(ith),nnew(ith)) modified contraction coefficients
!idim1                       :  first dimension
!idim2                       :  second dimension
!ich                         : index to be changed
!nolds1,nolds2,nolds3,nolds4 : old dimensions
!nnew1,nnew2,nnew3,nnew4     : new dimensions
!array1                      : array of size (nolds1,nolds2,nolds3,nolds4)
!array2                      : array of size (nnew1,nnew2,nnew3,nnew4)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idim1, idim2, ich, nolds1, nolds2, nolds3, nolds4, nnew1, nnew2, nnew3, nnew4
real(kind=wp), intent(in) :: coeffs(idim1,idim2), array1(nolds1,nolds2,nolds3,nolds4)
real(kind=wp), intent(out) :: array2(nnew1,nnew2,nnew3,nnew4)
integer(kind=iwp) :: ind1, ind2, ind3, ind4, ind5

!write(u6,*) 'begin trans ',ich
!write(u6,'(8I5)') nolds1,nolds2,nolds3,nolds4,nnew1,nnew2,nnew3,nnew4
array2(:,:,:,:) = Zero
if (ich == 1) then
  do ind4=1,nnew4
    do ind3=1,nnew3
      do ind2=1,nnew2
        do ind5=1,nnew1
          do ind1=1,nolds1
            array2(ind5,ind2,ind3,ind4) = array2(ind5,ind2,ind3,ind4)+coeffs(ind1,ind5)*array1(ind1,ind2,ind3,ind4)
          end do
        end do
      end do
    end do
  end do
else if (ich == 2) then
  !write(6,*) 'transform second index '
  do ind4=1,nnew4
    do ind3=1,nnew3
      do ind5=1,nnew2
        do ind2=1,nolds2
          array2(:,ind5,ind3,ind4) = array2(:,ind5,ind3,ind4)+coeffs(ind2,ind5)*array1(1:nnew1,ind2,ind3,ind4)
        end do
      end do
    end do
  end do
  !write(u6,*) 'end  to transform second index '
else if (ich == 3) then
  do ind4=1,nnew4
    do ind5=1,nnew3
      do ind3=1,nolds3
        array2(:,:,ind5,ind4) = array2(:,:,ind5,ind4)+coeffs(ind3,ind5)*array1(1:nnew1,1:nnew2,ind3,ind4)
      end do
    end do
  end do
else if (ich == 4) then
  do ind5=1,nnew4
    do ind4=1,nolds4
      array2(:,:,:,ind5) = array2(:,:,:,ind5)+coeffs(ind4,ind5)*array1(1:nnew1,1:nnew2,1:nnew3,ind4)
    end do
  end do
end if
!write(u6,*) 'end  trans '

return

end subroutine trans_amfi
