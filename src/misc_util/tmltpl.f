************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      subroutine tmltpl(inp,lpole,maxlab,
     *                  labs,
     *                  ndim,prvec,t,temp)
c
c     Purpose: transformation of maxlab cartesian l-th
c              moments into the corresponding l-pole cartesian
c              moments; l=lpole
c
c
c     inp                    if inp.eq.0 the t matrix will be
c                            calculated. Otherwise the t matrix
c                            is transferred from the previous
c                            call
c     lpole                  l value for the l-pole moment
c     maxlab                 is the number of cartesian components
c                            of the l-pole moment
c     labs(1:maxlab)         labels for components of the l-pole
c                            moment
c     ndim                   defines the number of rows which will
c                            be transformed in the table
c                            prvec(1:ndim,1:maxlab)
c     t(1:maxlab,1:maxlab)   is the transformation matrix generated
c                            in this program for lpole=2,3,4
c     temp(1:maxlab)         is a temporary strorage area
c
************************************************************************
      implicit real*8 (a-h,o-z)
      character*1 l1,l2,l3,l4
      character*14 l14
      character*16 labs(1:maxlab)
      dimension prvec(1:ndim,1:maxlab)
      dimension t(1:maxlab,1:maxlab)
      dimension temp(1:maxlab)
      dimension irr(1:3,1:3),ilab(1:6,1:3),irrrr(1:6,1:3)
c
      k=0 ! dummy initialize
      l=0 ! dummy initialize
c
c     building of transformation matrices for the transformation
c     from cartesina l-th moments to cartesian l-pole moments
c     limited to l=2,3, and 4
c
      if (inp.eq.1) go to 98
      go to (100,200,300), lpole-1
c
c     quadrupole moments
c
  100 continue
c
      do 101 i=1,maxlab
        do 102 j=1,maxlab
          t(i,j)=0.0d+00
  102   continue
        t(i,i)=t(i,i)+1.5d+00
        read (labs(i),'(a14,2a1)') l14,l1,l2
        if (l1.eq.l2) then
          t(i,1)=t(i,1)-0.5d+00
          t(i,4)=t(i,4)-0.5d+00
          t(i,6)=t(i,6)-0.5d+00
        endif
  101 continue
      go to 99
c
c     octupole moments
c
  200 continue
c
      do 205 i=1,3
        do 206 j=1,3
          irr(i,j)=0
  206   continue
        irr(i,i)=2
  205 continue
c
      do 201 i=1,maxlab
        do 202 j=1,maxlab
          t(i,j)=0.0d+00
  202   continue
        t(i,i)=t(i,i)+2.5d+00
        read (labs(i),'(a13,3a1)') l14,l1,l2,l3
        if (l1.eq.l2) then
          do 203 i1=1,3
            do 204 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  204       continue
  203     continue
          k=0
          if (l3.eq.'X') k=1
          if (l3.eq.'Y') k=2
          if (l3.eq.'Z') k=3
          do 207 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ind=(3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.5d+00
  207     continue
        endif
        if (l2.eq.l3) then
          do 208 i1=1,3
            do 209 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  209       continue
  208     continue
          k=0
          if (l1.eq.'X') k=1
          if (l1.eq.'Y') k=2
          if (l1.eq.'Z') k=3
          do 210 j=1,3
            ilab(j,k)=irr(j,k)+1
            ind=(3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.5d+00
  210     continue
        endif
        if (l1.eq.l3) then
          do 211 i1=1,3
            do 212 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  212       continue
  211     continue
          k=0
          if (l2.eq.'X') k=1
          if (l2.eq.'Y') k=2
          if (l2.eq.'Z') k=3
          do 213 j=1,3
            ilab(j,k)=irr(j,k)+1
            ind=(3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.5d+00
  213     continue
        endif
  201 continue
      go to 99
c
c     hexadexapole moments
c
  300 continue
c
      do 305 i=1,3
        do 306 j=1,3
          irr(i,j)=0
  306   continue
        do 307 j=1,6
          irrrr(j,i)=0
  307   continue
        irr(i,i)=2
  305 continue
      irrrr(1,1)=4
      irrrr(2,1)=2
      irrrr(2,2)=2
      irrrr(3,1)=2
      irrrr(3,3)=2
      irrrr(4,2)=4
      irrrr(5,2)=2
      irrrr(5,3)=2
      irrrr(6,3)=4
c
      do 301 i=1,maxlab
        do 302 j=1,maxlab
          t(i,j)=0.0d+00
  302   continue
        t(i,i)=t(i,i)+4.375d+00
        read (labs(i),'(a12,4a1)') l14,l1,l2,l3,l4
        if (l1.eq.l2) then
          do 303 i1=1,3
            do 304 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  304       continue
  303     continue
          k=0
          l=0
          if (l3.eq.'X') k=1
          if (l3.eq.'Y') k=2
          if (l3.eq.'Z') k=3
          if (l4.eq.'X') l=1
          if (l4.eq.'Y') l=2
          if (l4.eq.'Z') l=3
          do 308 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  308     continue
        endif
        if (l1.eq.l3) then
          do 309 i1=1,3
            do 310 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  310       continue
  309     continue
          k=0
          l=0
          if (l2.eq.'X') k=1
          if (l2.eq.'Y') k=2
          if (l2.eq.'Z') k=3
          if (l4.eq.'X') l=1
          if (l4.eq.'Y') l=2
          if (l4.eq.'Z') l=3
          do 311 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  311     continue
        endif
        if (l1.eq.l4) then
          do 312 i1=1,3
            do 313 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  313       continue
  312     continue
          k=0
          l=0
          if (l2.eq.'X') k=1
          if (l2.eq.'Y') k=2
          if (l2.eq.'Z') k=3
          if (l3.eq.'X') l=1
          if (l3.eq.'Y') l=2
          if (l3.eq.'Z') l=3
          do 314 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  314     continue
        endif
        if (l2.eq.l3) then
          do 315 i1=1,3
            do 316 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  316       continue
  315     continue
          k=0
          l=0
          if (l1.eq.'X') k=1
          if (l1.eq.'Y') k=2
          if (l1.eq.'Z') k=3
          if (l4.eq.'X') l=1
          if (l4.eq.'Y') l=2
          if (l4.eq.'Z') l=3
          do 317 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  317     continue
        endif
        if (l2.eq.l4) then
          do 318 i1=1,3
            do 319 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  319       continue
  318     continue
          k=0
          l=0
          if (l1.eq.'X') k=1
          if (l1.eq.'Y') k=2
          if (l1.eq.'Z') k=3
          if (l3.eq.'X') l=1
          if (l3.eq.'Y') l=2
          if (l3.eq.'Z') l=3
          do 320 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  320     continue
        endif
        if (l3.eq.l4) then
          do 321 i1=1,3
            do 322 i2=1,3
              ilab(i1,i2)=irr(i1,i2)
  322       continue
  321     continue
          k=0
          l=0
          if (l1.eq.'X') k=1
          if (l1.eq.'Y') k=2
          if (l1.eq.'Z') k=3
          if (l2.eq.'X') l=1
          if (l2.eq.'Y') l=2
          if (l2.eq.'Z') l=3
          do 323 j=1,3
            ilab(j,k)=ilab(j,k)+1
            ilab(j,l)=ilab(j,l)+1
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)-0.625d+00
  323     continue
        endif
        do 324 i1=1,6
          do 325 i2=1,3
            ilab(i1,i2)=irrrr(i1,i2)
  325     continue
  324   continue
        if (l1.eq.l2.and.l3.eq.l4) then
          do 326 j=1,6
            f=0.125d+00
            if (j.eq.2.or.j.eq.3.or.j.eq.5) f=2.0d+00*f
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)+f
  326     continue
        endif
        if (l1.eq.l3.and.l2.eq.l4) then
          do 327 j=1,6
            f=0.125d+00
            if (j.eq.2.or.j.eq.3.or.j.eq.5) f=2.0d+00*f
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)+f
  327     continue
        endif
        if (l2.eq.l3.and.l1.eq.l4) then
          do 328 j=1,6
            f=0.125d+00
            if (j.eq.2.or.j.eq.3.or.j.eq.5) f=2.0d+00*f
            ind=(4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
            T(i,ind)=T(i,ind)+f
  328     continue
        endif
  301 continue
c
   99 continue
c
c     print the transformation matrix
c
c      write (*,'(//1x,a,i2/)') 'transformation matrix:  lpole=',lpole
c      do 401 i=1,maxlab
c        write (*,'(15f7.3)') (t(i,j),j=1,maxlab)
c  401 continue
c
   98 continue
c
c     transform cartesian moment to multipole moments
c
      do 1 icount=1,ndim
        do 2 k=1,maxlab
          temp(k)=prvec(icount,k)
    2   continue
        do 3 k=1,maxlab
          sum=+0.0d+00
          do 4 l=1,maxlab
            sum=sum+t(k,l)*temp(l)
    4     continue
          prvec(icount,k)=sum
    3   continue
    1 continue
c
      return
#ifdef _WARNING_WORKAROUND_
      if (.false.) call unused_character(l14)
#endif
      end
