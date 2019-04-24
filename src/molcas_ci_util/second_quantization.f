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
* Copyright (C) 2013, Steven Vancoillie                                *
************************************************************************
      module second_quantization
!     Implements determinants of k electrons in n (spin)orbitals as combinations
!     (k,n) represented as bitstrings, using 32-bit integers.  Combinations are
!     represented in lexicographic ordering, and their rank can be computed using
!     a binomial table. The 32st bit is used for to indicate sign. all-one-bits
!     are used to indicate annihilation (destroyed state). Excitation operators
!     can be applied to the determinants and return the resulting determinant.

!     Additionally, a wavefunction is represtented by a structure (derived type)
!     consisting of determinant coefficients resulting from the direct product of
!     separate alpha and beta determinants.

!     Steven Vancoillie, November 2013, Lund

      implicit none
      integer, save :: onebits(0:255), ranktbl(0:255,64)

      contains

      recursive integer function gcd(i,j) result (k)
      implicit none
      integer, intent(in) :: i, j
      if (j.eq.0) then
        k = i
      else
        k = gcd(j,mod(i,j))
      end if
      end function gcd

      integer function binom_coef(k,n)
      ! compute binomial coefficient
      ! returns #k-subsets out of n choices
      implicit none
      integer, intent(in) :: k, n
      integer :: i, frac(2), div
      if (k .gt. n) then
        binom_coef = 0
      else
        frac = 1
        do i=1,k
          frac(1) = frac(1) * (n-k+i)
          frac(2) = frac(2) * i
          div = gcd(frac(1),frac(2))
          if (div.gt.1) then
            frac = frac / div
          end if
        end do
        binom_coef = frac(1)/frac(2)
      end if
      end function binom_coef

      integer function lex_init(k,n) result(c)
      ! initialize the first combination
      implicit none
      integer, intent(in) :: k, n
      if (k .gt. n) then
        c = 0
      else
        c = 2**k - 1
      end if
      end function lex_init

      integer function lex_next(c) result(d)
      ! generate the next lexicographic bit combination
      implicit none
      integer, intent(in) :: c
      integer :: t, s
      t = ior(c,c-1)+1
      s = iand(not(t-1),t)-1
      d = ior(t,ishft(s,-(trailz_(c)+1)))
      end function lex_next

      integer function trailz_(c)
      implicit none
      integer, intent(in) :: c
#include "compiler_features.h"
#ifdef TRAILING_ZEROS
      trailz_ = trailz(c)
#else
      ! simple trailing zero computation, only used when the trailz
      ! intrinsic is not supported.
      integer :: t
      t = c
      trailz_ = 0
      if (iand(t,z'FFFF') == 0) then
        t = ishft(t,-16)
        trailz_ = trailz_ + 16
      end if
      if (iand(t,z'FF') == 0) then
        t = ishft(t,-8)
        trailz_ = trailz_ + 8
      end if
      if (iand(t,z'F') == 0) then
        t = ishft(t,-4)
        trailz_ = trailz_ + 4
      end if
      if (iand(t,z'3') == 0) then
        t = ishft(t,-2)
        trailz_ = trailz_ + 2
      end if
      if (iand(t,z'1') == 0) then
        t = ishft(t,-1)
        trailz_ = trailz_ + 1
      end if
      trailz_ = trailz_ + iand(not(t),1)
#endif
      end function trailz_

      subroutine rank_init
      ! initializes a rank table indexed as:
      !      1 2   3    4        8
      !      1 1-8 1-16 1-24 ... 1-56
      !  0
      !  1
      !  ...
      !  255
      ! the row index is a hash table using the
      ! bit representation of a byte section.
      ! the columns refer to the position of
      ! a byte section and the number of possible
      ! preceding one-bits
      integer :: irow, icol, ibyte, ioffset
      integer :: ibit, ipos
      integer :: irank

      do irow = 0, 255
        onebits(irow) = 0
        do ibit = 0, 7
          if (btest(irow,ibit)) then
            onebits(irow) = onebits(irow) + 1
          end if
        end do
      end do

      do irow = 0, 255
        ! compute rank for first byte
        irank = 0
        ipos = 0
        do ibit = 0, 7
          if (btest(irow,ibit)) then
            ipos = ipos + 1
            irank = irank + binom_coef(ipos,ibit)
          end if
        end do
        icol = 1
        ranktbl(irow,icol) = irank
        ! compute rank for subsequent bytes
        ! depending on the number of preceding
        ! one-bits in the previous bytes
        do ibyte = 2, 4
          do ioffset = 0, 8*(ibyte-1)
            icol = icol + 1
            ! compute rank
            irank = 0
            ipos = ioffset
            do ibit = 0, 7
              if (btest(irow,ibit)) then
                ipos = ipos + 1
                irank = irank + binom_coef(ipos,ibit+8*(ibyte-1))
              end if
            end do
            ranktbl(irow,icol) = irank
          end do
        end do
      end do
      end subroutine rank_init

      integer function lexrank(c)
      integer, intent(in) :: c
      integer :: byte(4), ones(4)
      lexrank = 0
      if (c.eq.-1) return
      byte(1) = ibits(c, 0,8)
      byte(2) = ibits(c, 8,8)
      byte(3) = ibits(c,16,8)
      byte(4) = ibits(c,24,6)
      ones(1) = onebits(byte(1))
      ones(2) = onebits(byte(2)) + ones(1)
      ones(3) = onebits(byte(3)) + ones(2)
      lexrank = 1 + ranktbl(byte(1), 1) +
     & ranktbl(byte(2), 2+ones(1)) +
     & ranktbl(byte(3),11+ones(2)) +
     & ranktbl(byte(4),28+ones(3))
      end function lexrank

      integer function fase(c)
      integer, intent(in) :: c
      if (btest(c,31)) then
        fase=-1
      else
        fase=1
      end if
      end function fase

      integer function ex1(p,q,c) result(d)
      implicit none
      integer, intent(in) :: p, q, c
      integer :: t
      d = c
      if (.not.btest(d,q-1)) then
        d = -1
        return
      end if
      d = ibclr(d,q-1)
      if (btest(d,p-1)) then
        d = -1
        return
      end if
      d = ibset(d,p-1)
      if(p.gt.q) then
        t=ibits(d,q,p-q-1)
      else if (p.lt.q) then
        t=ibits(d,p,q-p-1)
      else
        return
      end if
      do while (t.ne.0)
        d = ieor(ishft(iand(t,1),31),d)
        t = ishft(t,-1)
      end do
      end function ex1

      integer function ann(p,c) result(d)
      ! operates on determinant c with (a_p)
      ! and returns the resulting determinant with
      ! the sign in the highest bit (1 = negative)
      integer, intent(in) :: p, c
      integer :: t, b_1111, z_6996
      data b_1111 /b'1111'/
      data z_6996 /z'6996'/
      if (.not.btest(c,p-1)) then
        d = -1
        return
      end if
      d = ibclr(c,p-1)
      t = ibits(d,0,p-1)
      t = ieor(t,ishft(t,-16))
      t = ieor(t,ishft(t,-8))
      t = ieor(t,ishft(t,-4))
      t = iand(t,b_1111)
      t = iand(ishft(z_6996,-t),1)
      d = ieor(ishft(t,31),d)
      end function ann

      integer function cre(p,c) result(d)
      ! operates on determinant c with (a+_p)
      ! and returns the resulting determinant with
      ! the sign in the highest bit (1 = negative)
      integer, intent(in) :: p, c
      integer :: t, b_1111, z_6996
      data b_1111 /b'1111'/
      data z_6996 /z'6996'/
      if (btest(c,p-1)) then
        d = -1
        return
      end if
      d = ibset(c,p-1)
      t = ibits(d,0,p-1)
      t = ieor(t,ishft(t,-16))
      t = ieor(t,ishft(t,-8))
      t = ieor(t,ishft(t,-4))
      t = iand(t,b_1111)
      t = iand(ishft(z_6996,-t),1)
      d = ieor(ishft(t,31),d)
      end function cre

      integer function ann2(p,q,c) result(d)
      ! operates on determinant c with (a_p a_q)
      ! and returns the resulting determinant with
      ! the sign in the highest bit (1 = negative)
      implicit none
      integer, intent(in) :: p, q, c
      integer :: t, b_1111, z_6996
      data b_1111 /b'1111'/
      data z_6996 /z'6996'/
      if (.not.(btest(c,q-1).and.btest(c,p-1))) then
        d = -1
        return
      end if
      d = ibclr(c,q-1)
      d = ibclr(d,p-1)
      if(q.lt.p) then
        t=ibits(d,q,p-q-1)
      else if (p.lt.q) then
        t=ibits(d,p,q-p-1)
        d=ieor(ishft(1,31),d)
      else
        d = -1
        return
      end if
      t = ieor(t,ishft(t,-16))
      t = ieor(t,ishft(t,-8))
      t = ieor(t,ishft(t,-4))
      t = iand(t,b_1111)
      t = iand(ishft(z_6996,-t),1)
      d = ieor(ishft(t,31),d)
      end function ann2

      integer function cre2(p,q,c) result(d)
      ! operates on determinant c with (a+_p a+_q)
      ! and returns the resulting determinant with
      ! the sign in the highest bit (1 = negative)
      implicit none
      integer, intent(in) :: p, q, c
      integer :: t, b_1111, z_6996
      data b_1111 /b'1111'/
      data z_6996 /z'6996'/
      if ((btest(c,q-1).or.btest(c,p-1))) then
        d = -1
        return
      end if
      d = ibset(c,q-1)
      d = ibset(d,p-1)
      if(q.lt.p) then
        t=ibits(d,q,p-q-1)
        d=ieor(ishft(1,31),d)
      else if (p.lt.q) then
        t=ibits(d,p,q-p-1)
      else
        d = -1
        return
      end if
      t = ieor(t,ishft(t,-16))
      t = ieor(t,ishft(t,-8))
      t = ieor(t,ishft(t,-4))
      t = iand(t,b_1111)
      t = iand(ishft(z_6996,-t),1)
      d = ieor(ishft(t,31),d)
      end function cre2

#ifdef __PGI
      ! PGI Fortran compiler version of ibits does not work correctly when
      ! LEN == 0
      integer function ibits(c,pos,len) result(d)
      implicit none
      integer, intent(in) :: c, pos, len
      d = ishft(not(0),-(bit_size(d)-len))
      d = iand(ishft(c,-pos),d)
      return
      end function
#endif

      end module second_quantization
