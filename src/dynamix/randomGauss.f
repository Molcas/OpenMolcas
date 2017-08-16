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
C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8


      SUBROUTINE RandomGauss(ValMean,Sigma,iseed,nflag,buffer,Val)
      IMPLICIT REAL*8 (a-h,o-z)

C       ValMean is the mean, and sigma is the standard deviation.
C       nFlag is a binary (0,1) variable for returning the appropiate random value

C When x and y are two variables from [0, 1), uniformly
C distributed, then
C
C    cos(2*pi*x)*sqrt(-2*log(1-y))
C    sin(2*pi*x)*sqrt(-2*log(1-y))
C
C are two *independent* variables with normal distribution
C (mu = 0, sigma = 1).
C (Lambert Meertens)
C (corrected version; bug discovered by Mike Miller, fixed by LM)
      IF (nFlag.eq.0) THEN

          alpha = abs(fRandom_Molcas(iseed))
          beta  = abs(fRandom_Molcas(iseed))

          PI=4.D0*ATAN(1.D0)
          X2pi = alpha * (2.d0*Pi)
          G2rad = sqrt(-2.d0 * log(1.d0 - beta ))

          Z1 = cos(X2pi) * G2rad
          Z2 = sin(X2pi) * G2rad
          Val = ValMean + Z1*Sigma
          buffer = ValMean + Z2*Sigma
          nFlag = 1
      ELSE
         Val = buffer
         nFlag = 0


      END IF
c      write (6,*) 'PI:' , PI
c      write (6,*) 'X2pi:', X2pi
c      write (6,*) 'G2grad: ', G2rad
c      write (6,*) 'buffer:', buffer
c      write (6,*) 'Z1 and Z2', Z1, Z2
c      write (6,*) 'Alpha & Beta', alpha, beta
c      write (6,*) 'VAL into RANDOM:', val
      RETURN

      END

      Subroutine getSeed(iseed)

      External  datimx
      character Line*72,l*1
      integer*8 ihours,minutes,seconds,day

      call datimx(Line)
      read (Line,'(a8,i2,a1,i2,a1,i2,a1,i2)') l,day,l,ihours, l,minutes,
     & l, seconds

      iseed = int((ihours*3600)+(minutes*60)+seconds+(day*86400),
     &            kind(iseed))


      Return

      End
