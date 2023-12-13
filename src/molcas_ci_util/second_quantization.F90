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
! Copyright (C) 2013, Steven Vancoillie                                *
!***********************************************************************

module second_quantization
! Implements determinants of k electrons in n (spin)orbitals as combinations
! (k,n) represented as bitstrings, using 32-bit integers.  Combinations are
! represented in lexicographic ordering, and their rank can be computed using
! a binomial table. The 32nd bit is used for to indicate sign. all-one-bits
! are used to indicate annihilation (destroyed state). Excitation operators
! can be applied to the determinants and return the resulting determinant.
!
! Additionally, a wavefunction is represtented by a structure (derived type)
! consisting of determinant coefficients resulting from the direct product of
! separate alpha and beta determinants.
!
! Steven Vancoillie, November 2013, Lund

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: onebits(0:255), ranktbl(0:255,64)
integer(kind=iwp), parameter :: b_1111 = int(b'1111',kind=iwp), z_6996 = int(z'6996',kind=iwp)

public :: ann, ann2, binom_coef, cre, cre2, ex1, fase, lex_init, lex_next, lexrank, rank_init

contains

#include "compiler_features.h"

#ifndef TRAILING_ZEROS
#define trailz trailz_
function trailz_(c) result(d)
! simple trailing zero computation, only used when the trailz intrinsic is not supported.

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: c
  integer(kind=iwp) :: t

  d = 0
  do t=0,bit_size(c)-1
    if (btest(c,t)) exit
    d = d+1
  end do

end function trailz_
#endif

#ifndef IBITS_LEN_ZERO
#define ibits ibits_
function ibits_(c,pos,len) result(d)
! Workaround for ibits not working with zero length

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: c, pos, len

  d = ishft(not(0),-(bit_size(d)-len))
  d = iand(ishft(c,-pos),d)

end function ibits_
#endif

recursive function gcd(i,j) result(k)

  integer(kind=iwp) :: k
  integer(kind=iwp), intent(in) :: i, j

  if (j == 0) then
    k = i
  else
    k = gcd(j,mod(i,j))
  end if

end function gcd

function binom_coef(k,n)
! compute binomial coefficient
! returns #k-subsets out of n choices

  integer(kind=iwp) :: binom_coef
  integer(kind=iwp), intent(in) :: k, n
  integer(kind=iwp) :: i, frac(2), div

  if (k > n) then
    binom_coef = 0
  else
    frac = 1
    do i=1,k
      frac(1) = frac(1)*(n-k+i)
      frac(2) = frac(2)*i
      div = gcd(frac(1),frac(2))
      if (div > 1) then
        frac = frac/div
      end if
    end do
    binom_coef = frac(1)/frac(2)
  end if

end function binom_coef

function lex_init(k,n) result(c)
! initialize the first combination

  integer(kind=iwp) :: c
  integer(kind=iwp), intent(in) :: k, n

  if (k > n) then
    c = 0
  else
    c = 2**k-1
  end if

end function lex_init

function lex_next(c) result(d)
! generate the next lexicographic bit combination

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: c
  integer(kind=iwp) :: t, s

  t = ior(c,c-1)+1
  s = iand(not(t-1),t)-1
  d = ior(t,ishft(s,-min(trailz(c)+1,storage_size(s))))

end function lex_next

subroutine rank_init()
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

  integer(kind=iwp) :: irow, icol, ibyte, ioffset, ibit, ipos, irank

  do irow=0,255
    onebits(irow) = 0
    do ibit=0,7
      if (btest(irow,ibit)) then
        onebits(irow) = onebits(irow)+1
      end if
    end do
  end do

  do irow=0,255
    ! compute rank for first byte
    irank = 0
    ipos = 0
    do ibit=0,7
      if (btest(irow,ibit)) then
        ipos = ipos+1
        irank = irank+binom_coef(ipos,ibit)
      end if
    end do
    icol = 1
    ranktbl(irow,icol) = irank
    ! compute rank for subsequent bytes
    ! depending on the number of preceding
    ! one-bits in the previous bytes
    do ibyte=2,4
      do ioffset=0,8*(ibyte-1)
        icol = icol+1
        ! compute rank
        irank = 0
        ipos = ioffset
        do ibit=0,7
          if (btest(irow,ibit)) then
            ipos = ipos+1
            irank = irank+binom_coef(ipos,ibit+8*(ibyte-1))
          end if
        end do
        ranktbl(irow,icol) = irank
      end do
    end do
  end do

end subroutine rank_init

function lexrank(c)

  integer(kind=iwp) :: lexrank
  integer(kind=iwp), intent(in) :: c
  integer(kind=iwp) :: byte(4), ones(4)

  lexrank = 0
  if (c == -1) return
  byte(1) = ibits(c,0,8)
  byte(2) = ibits(c,8,8)
  byte(3) = ibits(c,16,8)
  byte(4) = ibits(c,24,6)
  ones(1) = onebits(byte(1))
  ones(2) = onebits(byte(2))+ones(1)
  ones(3) = onebits(byte(3))+ones(2)
  lexrank = 1+ranktbl(byte(1),1)+ranktbl(byte(2),2+ones(1))+ranktbl(byte(3),11+ones(2))+ranktbl(byte(4),28+ones(3))

end function lexrank

function fase(c)

  integer(kind=iwp) :: fase
  integer(kind=iwp), intent(in) :: c

  if (btest(c,31)) then
    fase = -1
  else
    fase = 1
  end if

end function fase

function ex1(p,q,c) result(d)

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: p, q, c
  integer(kind=iwp) :: t

  d = c
  if (.not. btest(d,q-1)) then
    d = -1
    return
  end if
  d = ibclr(d,q-1)
  if (btest(d,p-1)) then
    d = -1
    return
  end if
  d = ibset(d,p-1)
  if (p > q) then
    t = ibits(d,q,p-q-1)
  else if (p < q) then
    t = ibits(d,p,q-p-1)
  else
    return
  end if
  do while (t /= 0)
    d = ieor(ishft(iand(t,1),31),d)
    t = ishft(t,-1)
  end do

end function ex1

function ann(p,c) result(d)
! operates on determinant c with (a_p)
! and returns the resulting determinant with
! the sign in the highest bit (1 = negative)

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: p, c
  integer(kind=iwp) :: t

  if (.not. btest(c,p-1)) then
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

function cre(p,c) result(d)
! operates on determinant c with (a+_p)
! and returns the resulting determinant with
! the sign in the highest bit (1 = negative)

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: p, c
  integer(kind=iwp) :: t

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

function ann2(p,q,c) result(d)
! operates on determinant c with (a_p a_q)
! and returns the resulting determinant with
! the sign in the highest bit (1 = negative)

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: p, q, c
  integer(kind=iwp) :: t

  if (.not. (btest(c,q-1) .and. btest(c,p-1))) then
    d = -1
    return
  end if
  d = ibclr(c,q-1)
  d = ibclr(d,p-1)
  if (q < p) then
    t = ibits(d,q,p-q-1)
  else if (p < q) then
    t = ibits(d,p,q-p-1)
    d = ieor(ishft(1,31),d)
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

function cre2(p,q,c) result(d)
! operates on determinant c with (a+_p a+_q)
! and returns the resulting determinant with
! the sign in the highest bit (1 = negative)

  integer(kind=iwp) :: d
  integer(kind=iwp), intent(in) :: p, q, c
  integer(kind=iwp) :: t

  if ((btest(c,q-1) .or. btest(c,p-1))) then
    d = -1
    return
  end if
  d = ibset(c,q-1)
  d = ibset(d,p-1)
  if (q < p) then
    t = ibits(d,q,p-q-1)
    d = ieor(ishft(1,31),d)
  else if (p < q) then
    t = ibits(d,p,q-p-1)
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

end module second_quantization
