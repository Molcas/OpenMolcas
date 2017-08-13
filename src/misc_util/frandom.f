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
* Copyright (C) 1993, Per-Olof Widmark                                 *
*               2010, John Burkardt                                    *
*               2017, Morgane Vacher                                   *
************************************************************************
      Function fRandom_Molcas(iSeed)
************************************************************************
*                                                                      *
*     Generate a random number z, where 0<z<0.5d0**31                  *
*                                                                      *
*     calling arguments:                                               *
*     iSeed  : Type integer, input/output                              *
*              initial value of new chain                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark                                                     *
*     University of Lund, Sweden, 1993                                 *
*     M. Vacher, Uppsala 2017                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Parameter (IM1=134456, IA1=8121, IC1=28411)
      Parameter (IM2=243000, IA2=4561, IC2=51349)
      Parameter (IM3=259200, IA3=7141, IC3=54773)
      real*8 , parameter :: a = 1220703125.0D+00
      real*8 , save :: a1, a2
      integer  i
      integer , save :: ks = 0
      real*8 , save :: r23, r46, t23, t46
      real*8  t1, t2, t3, t4, x, x1, x2, z
      integer  iseed
      Real*8 fRandom_Molcas
      Character*8 rand_info
      Call GetEnvf('MOLCAS_RANDOM',rand_info)
      Call UpCase(rand_info)
      if (rand_info(1:3) .EQ. 'OLD') then
        IX0 = iSeed
        IX1 = MOD(IA1*IX0+IC1,IM1)
        IX2 = MOD(IA2*IX1+IC2,IM2)
        IX3 = MOD(IA3*IX2+IC3,IM3)
        fRandom_Molcas = (DBLE(IX1)+DBLE(IX2)/DBLE(IM2))/DBLE(IM1)
        iSeed = IX3
      else
************************************************************************
*     The number returned is a uniform pseudorandom value in the range *
*     (0,1). The algorithm uses the linear congruential generator:     *
*     X(K+1) = A * X(K)  mod 2^46                                      *
*     This scheme generates 2^44 numbers before repeating.             *
*----------------------------------------------------------------------*
*     Author: John Burkardt                                            *
*     Modified: 08 March 2010                                          *
*----------------------------------------------------------------------*
*     Reference:                                                       *
*     David Bailey, Eric Barszcz, John Barton, D Browning, Robert      *
*     Carter, Leonardo Dagum, Rod Fatoohi, Samuel Fineberg, Paul       *
*     Frederickson, Thomas Lasinski, Robert Schreiber, Horst Simon, V  *
*     Venkatakrishnan, Sisira Weeratunga, The NAS Parallel Benchmarks, *
*     RNR Technical Report RNR-94-007, March 1994.                     *
*     Donald Knuth, The Art of Computer Programming, Volume 2,         *
*     Seminumerical Algorithms, Third Edition, Addison Wesley, 1997,   *
*     ISBN: 0201896842, LC: QA76.6.K64.                                *
*----------------------------------------------------------------------*
*     Parameters:                                                      *
*     Input/output, real ( kind = 8 ) X, the seed.  X should be an     *
*     odd integer such that 1 <= X <= 2^46.                            *
*     Output, real (kind=8) Random_Molcas, the next pseudorandom value.*
************************************************************************
*  If this is the first call, compute
*    R23 = 2 ^ -23,
*    R46 = 2 ^ -46,
*    T23 = 2 ^ 23,
*    T46 = 2 ^ 46.
*  These are computed in loops, rather than by merely using the power operator,
*  in order to insure that the results are exact on all systems.
      x = real(iseed,kind(x))
      if ( ks == 0 ) then
        r23 = 1.0D+00
        r46 = 1.0D+00
        t23 = 1.0D+00
        t46 = 1.0D+00
        do i = 1, 23
          r23 = 0.5D+00 * r23
          t23 = 2.0D+00 * t23
        end do
        do i = 1, 46
          r46 = 0.50D+00 * r46
          t46 = 2.0D+00 * t46
        end do
*  Break A into two parts such that A = 2^23 * A1 + A2.
        t1 = r23 * a
        a1 = real ( int ( t1 ), kind = 8 )
        a2 = a - t23 * a1
        ks = 1
      end if
*  Deal with a 0 input value of X.
      if ( x == 0.0D+00 ) then
        x = 314159265.0D+00
      end if
*  Deal somewhat arbitrarily with negative input X.
      if ( x < 0.0D+00 ) then
        x = - x
      end if
*  Break X into two parts X1 and X2 such that:
*    X = 2^23 * X1 + X2,
*  then compute
*    Z = A1 * X2 + A2 * X1  (mod 2^23)
*    X = 2^23 * Z + A2 * X2  (mod 2^46).
      t1 = r23 * x
      x1 = real ( int ( t1 ), kind = 8 )
      x2 = x - t23 * x1
      t1 = a1 * x2 + a2 * x1
      t2 = real ( int ( r23 * t1 ), kind = 8 )
      z = t1 - t23 * t2
      t3 = t23 * z + a2 * x2
      t4 = real ( int ( r46 * t3 ), kind = 8 )
      x = t3 - t46 * t4
      fRandom_Molcas = r46 * x
      iseed = int(x,kind(iseed))
      endif
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
