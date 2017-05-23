************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine trans_amfi(
cbs   makes the transformation for the ich-th index
     *coeffs, !(nolds(ith),nnew(ith)) modified contraction coefficients
     *idim1,  !  first dimension
     *idim2,  !  second dimension
     *ich,    ! index to be changed
     *nolds1,nolds2,nolds3,nolds4,  ! old dimensions
     *nnew1,nnew2,nnew3,nnew4, ! new dimensions
     *array1, ! array of size (nolds1,nolds2,nolds3,nolds4)
     *array2  ! array of size (nnew1,nnew2,nnew3,nnew4)
     *)
      implicit real*8 (a-h,o-z)
      dimension coeffs(idim1,idim2),
     *array1(nolds1,nolds2,nolds3,nolds4),
     *array2(nnew1,nnew2,nnew3,nnew4)
c     write(6,*) 'begin trans ' ,ich
c     write(6,'(8I5)') nolds1,nolds2,nolds3,nolds4,
c    *nnew1,nnew2,nnew3,nnew4
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind4)=0d0
      enddo
      enddo
      enddo
      enddo
      if (ich.eq.1) then
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind5=1,nnew1
      do ind1=1,nolds1
      array2(ind5,ind2,ind3,ind4)=array2(ind5,ind2,ind3,ind4)+
     *coeffs(ind1,ind5)*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.2) then
c     write(6,*) 'transform second index '
      do ind4=1,nnew4
      do ind3=1,nnew3
      do ind5=1,nnew2
      do ind2=1,nolds2
      coeff=coeffs(ind2,ind5)
      do ind1=1,nnew1
      array2(ind1,ind5,ind3,ind4)=array2(ind1,ind5,ind3,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
c     write(6,*) 'end  to transform second index '
      elseif (ich.eq.3) then
      do ind4=1,nnew4
      do ind5=1,nnew3
      do ind3=1,nolds3
      coeff=coeffs(ind3,ind5)
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind5,ind4)=array2(ind1,ind2,ind5,ind4)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      elseif (ich.eq.4) then
      do ind5=1,nnew4
      do ind4=1,nolds4
      coeff=coeffs(ind4,ind5)
      do ind3=1,nnew3
      do ind2=1,nnew2
      do ind1=1,nnew1
      array2(ind1,ind2,ind3,ind5)=array2(ind1,ind2,ind3,ind5)+
     *coeff*array1(ind1,ind2,ind3,ind4)
      enddo
      enddo
      enddo
      enddo
      enddo
      endif
c     write(6,*) 'end  trans '
      return
      end
