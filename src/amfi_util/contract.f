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
      subroutine contract(
     *coeffs1, !(nprim(1),ncont(1)) modified contraction coefficients
     *coeffs2, !(nprim(2),ncont(2)) modified contraction coefficients
     *coeffs3, !(nprim(3),ncont(3)) modified contraction coefficients
     *coeffs4, !(nprim(4),ncont(4)) modified contraction coefficients
     *ncont,   ! i-th element is number of contracted functions i. index
     *nprim,   ! i-th element is number of primitive functions  i. index
cbs  array one contains at the beginning the uncontracted integrals
     *arr1,  ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *arr2   ! array of size (nprim(1)*nprim(2)*nprim(3)*nprim(4))
     *)
      implicit real*8 (a-h,o-z)
      dimension coeffs1(*),coeffs2(*),coeffs3(*),coeffs4(*),
     *arr1(*),arr2(*),ncont(4),nprim(4),nolds(4),nnew(4)
C
cbs   makes four indextransformations in a row....
cbs   try to find out, which indices should be transformed first...
c
      ratio1=DBLE(nprim(1))/DBLE(ncont(1))
      ratio2=DBLE(nprim(2))/DBLE(ncont(2))
      ratio3=DBLE(nprim(3))/DBLE(ncont(3))
      ratio4=DBLE(nprim(4))/DBLE(ncont(4))
      do IBM=1,4
      nolds(IBM)=nprim(IBM)
      nnew(IBM)=nprim(IBM)
      enddo
cbs   determine first, second,third and last index
************************************************************************
cbs   determine the first
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then
         ifirst=1
         ratio1=0d0
         nnew(ifirst)=ncont(ifirst)
         call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio2) then
         ifirst=2
         ratio2=0d0
         nnew(ifirst)=ncont(ifirst)
         call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio3) then
         ifirst=3
         ratio3=0d0
         nnew(ifirst)=ncont(ifirst)
         call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio4) then
         ifirst=4
         ratio4=0d0
         nnew(ifirst)=ncont(ifirst)
         call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else
         ifirst=0
         write (6,*) 'Contract: you should not be here!'
         call abend()
      endif
      nolds(ifirst)=nnew(ifirst)
************************************************************************
cbs   determine the second
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then
         isec=1
         ratio1=0d0
         nnew(isec)=ncont(isec)
         call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio2) then
         isec=2
         ratio2=0d0
         nnew(isec)=ncont(isec)
         call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio3) then
         isec=3
         ratio3=0d0
         nnew(isec)=ncont(isec)
         call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio4) then
         isec=4
         ratio4=0d0
         nnew(isec)=ncont(isec)
         call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else
         isec=0
         write (6,*) 'Contract: you should not be here!'
         call abend()
      endif
      nolds(isec)=nnew(isec)
************************************************************************
cbs   determine the third
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then
         ithird=1
         ratio1=0d0
         nnew(ithird)=ncont(ithird)
         call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio2) then
         ithird=2
         ratio2=0d0
         nnew(ithird)=ncont(ithird)
         call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio3) then
         ithird=3
         ratio3=0d0
         nnew(ithird)=ncont(ithird)
         call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else if (xmax.eq.ratio4) then
         ithird=4
         ratio4=0d0
         nnew(ithird)=ncont(ithird)
         call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr1,arr2)
      else
         ithird=0
         write (6,*) 'Contract: you should not be here!'
         call abend()
      endif
      nolds(ithird)=nnew(ithird)
************************************************************************
cbs   determine the last
      xmax=max(ratio1,ratio2,ratio3,ratio4)
      if (xmax.eq.ratio1) then
         ifourth=1
         ratio1=0d0
         nnew(ifourth)=ncont(ifourth)
         call trans_amfi(coeffs1,nprim(1),ncont(1),1,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio2) then
         ifourth=2
         ratio2=0d0
         nnew(ifourth)=ncont(ifourth)
         call trans_amfi(coeffs2,nprim(2),ncont(2),2,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio3) then
         ifourth=3
         ratio3=0d0
         nnew(ifourth)=ncont(ifourth)
         call trans_amfi(coeffs3,nprim(3),ncont(3),3,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else if (xmax.eq.ratio4) then
         ifourth=4
         ratio4=0d0
         nnew(ifourth)=ncont(ifourth)
         call trans_amfi(coeffs4,nprim(4),ncont(4),4,nolds(1),nolds(2),
     *      nolds(3),nolds(4),nnew(1),nnew(2),nnew(3),nnew(4),arr2,arr1)
      else
         ifourth=0
         write (6,*) 'Contract: you should not be here!'
         call abend()
      endif
cbs   contracted integrals are now on
cbs   arr1(ncont1,ncont2,ncont3,ncont4)
*
      return
      end
