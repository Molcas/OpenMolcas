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
C!-----------------------------------------------------------------------!
C!
C!VV Module RandomMod
C!
C! Description:
C!   This random number generator originally appeared in "Toward a Universal
C!   Random Number Generator" by George Marsaglia and Arif Zaman.
C!   Florida State University Report: FSU-SCRI-87-50 (1987)
C!
C!   It was later modified by F. James and published in "A Review of Pseudo-
C!   random Number Generators"
C!
C!   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
C!         (However, a newly discovered technique can yield
C!           a period of 10^600. But that is still in the development stage.)
C!
C!   It passes ALL of the tests for random number generators and has a period
C!   of 2^144, is completely portable (gives bit identical results on all
C!   machines with at least 24-bit mantissas in the floating point
C!   representation).
C!
C!   The algorithm is a combination of a Fibonacci sequence (with lags of 97
C!   and 33, and operation "subtraction plus one, modulo one") and an
C!   "arithmetic sequence" (using subtraction).
C!
C!   On a Vax 11/780, this random number generator can produce a number in
C!   13 microseconds.
C!
C!-----------------------------------------------------------------------!
C!
C!VV Implicit Real*8 ( a-h,o-z )

C!VV Contains

C!-----------------------------------------------------------------------!
C!
      Subroutine Rmarin(ij,kl)
C!
C! This is the initialization routine for the random number generator RANMAR()
C! NOTE: The seed variables can have values between:    0 .le. IJ .le. 31328
C!                                                      0 .le. KL .le. 30081
C! The random number sequences created by these two seeds are of sufficient
C! length to complete an entire calculation with. For example, if sveral
C! different groups are working on different parts of the same calculation,
C! each group could be assigned its own IJ seed. This would leave each group
C! with 30000 choices for the second seed. That is to say, this random
C! number generator can create 900 million different subsequences -- with
C! each subsequence having a length of approximately 10^30.
C!
C! Use IJ = 1802 & KL = 9373 to test the random number generator. The
C! subroutine RANMAR should be used to generate 20000 random numbers.
C! Then display the next six random numbers generated multiplied by 4096*4096
C! If the random number generator is working properly, the random numbers
C! should be:
C!           6533892.0  14220222.0  7275067.0
C!           6172232.0  8354498.0   10633180.0
C!
      Implicit Real*8 ( a-h,o-z )

      Real*8           U(97),c,cd,cm
      Integer           i97,j97
      Logical           test
      common /myrand/ U,c,cd,cm,i97,j97,test
      test = .false.
C!
      If (( ij.lt.0 ).or.( ij.gt.31328 ).or.
     &             ( kl.lt.0 ).or.( kl.gt.30081 )) Then
      write(6, '(a)')
     & 'The first random number seed ',
     & 'must have a value between 0 and 31328'
      write(6, '(a)')
     &  'The second seed must have a value between 0 and 30081'
      call abend()
      End If
C!
      i = mod(ij/177,177)+2
      j = mod(ij    ,177)+2
      k = mod(kl/169,178)+1
      l = mod(kl,    169)
C!
      Do ii = 1,97
      s = 0.0d0
      t = 0.5d0
      Do jj = 1,24
      m = mod(mod(i*j,179)*k,179)
      i = j
      j = k
      k = m
      l = mod(53*l+1,169)
      If (mod(l*m,64).ge.32) then
      s = s+t
      End If
      t = 0.5*t
      End Do
      U(ii) = s
      End Do
C!
      c  = 362436.0d0/16777216.0d0
      cd = 7654321.0d0/16777216.0d0
      cm = 16777213.0d0/16777216.0d0
C!
      i97 = 97
      j97 = 33
C!
      test = .true.
C!
      End
C!
C!-----------------------------------------------------------------------!


C!-----------------------------------------------------------------------!
C!
      Subroutine Ranmar(rvec,len)
C!
C! This is the random number generator proposed by George Marsaglia in
C! Florida State University Report: FSU-SCRI-87-50
C! It was slightly modified by F. James to produce an array of pseudorandom
C! numbers.
C!
      Implicit Real*8 ( a-h,o-z )
      Real*8            U(97),c,cd,cm
      Integer           i97,j97
      Logical           test
      Real*8 rvec(*)
      Real*8         uni
      Integer        ivec
      common /myrand/ U,c,cd,cm,i97,j97,test
C!
      If ( .not.test ) Then
      write(6, '(a)')
     &    ' Call the init routine (RMARIN) before calling RANMAR'
      call abend ()
      End If
C!
      Do ivec = 1,len
      uni = u(i97)-u(j97)
      if( uni.lt.0.0 ) uni = uni+1.0
      u(i97) = uni
      i97 = i97-1
      If ( i97.eq.0 ) i97 = 97
      j97 = j97-1
      If ( j97.eq.0 ) j97 = 97
      c = c-cd
      If ( c.lt.0.0 ) c = c+cm
      uni = uni-c
      If ( uni.lt.0.0 ) uni = uni+1.0
      rvec(ivec) = uni
      End Do
C!
      End
C!
C!-----------------------------------------------------------------------!

C!VV End Module RandomMod
C!
C!-----------------------------------------------------------------------!
