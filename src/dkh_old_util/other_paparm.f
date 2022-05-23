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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine other_param (ducoeffs,dkhorder,paramtype)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 26.06.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Jena)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,scrmax
      parameter (scrmax=100)
c
      real*8 ducoeffs(maxorder),scr(scrmax)
      character*(3) paramtype
      integer i
      if (maxorder.gt.scrmax) then
        write (stdout,1001) maxorder,scrmax
1001    format(5X,'Parameter  maxorder = ',I3,' is larger than ',I3,'.',
     *         //5X,'Increase size of scratch array "scr" in ',
     *         'subroutine "other_param"!',//5X,'STOP.')
        Call Abend
      endif
      scr(1)=1.0d0
      do 10 i=2,scrmax
        scr(i)=0.0d0
  10  continue
      if (paramtype.eq.'SQR') then
        scr(2) =0.50d0
        scr(4) =-0.125d0
        scr(6) =1.d0/16.d0
        scr(8) =-5.d0/128.d0
        scr(10)=7.d0/256.d0
        scr(12)=-21.d0/1024.d0
        scr(14)=33.d0/2048.d0
        scr(16)=-429.d0/32768.d0
        scr(18)=715.d0/65536.d0
        scr(20)=-2431.d0/262144.d0
        scr(22)=4199.d0/524288.d0
        scr(24)=-29393.d0/4194304.d0
        scr(26)=52003.d0/8388608.d0
        scr(28)=-185725.d0/33554432.d0
        scr(30)=334305.d0/67108864.d0
        scr(32)=-9694845.d0/2147483648.d0
        scr(34)=17678835.d0/4294967296.d0
        scr(36)=-64822395.d0/17179869184.d0
        scr(38)=119409675.d0/34359738368.d0
        scr(40)=-883631595.d0/274877906944.d0
        scr(42)=1641030105.d0/549755813888.d0
        scr(44)=-6116566755.d0/2199023255552.d0
        scr(46)=11435320455.d0/4398046511104.d0
        scr(48)=-171529806825.d0/70368744177664.d0
        if (dkhorder.ge.49) then
          write (stdout,1002)
1002      format (//2X,'Note:  Square-root param. works so far only up',
     *            ' to dkhorder = 48.',/2X)
          Call Abend
        endif
        goto 700
      endif
      if (paramtype.eq.'CAY') then
        scr(2)=0.50d0
        do 20 i=3,scrmax
          scr(i)=scr(i-1)/2.0d0
  20    continue
        goto 700
      endif
      if (paramtype.eq.'MCW') then
        scr(2)=0.50d0
        scr(4)=0.3750d0
        scr(6)=5.0d0/16.0d0
        scr(8)=35.0d0/128.0d0
        scr(10)=63.0d0/256.0d0
        scr(12)=231.0d0/1024.0d0
        scr(14)=429.0d0/2048.0d0
        scr(16)=6435.0d0/32768.0d0
        scr(18)=12155.0d0/65536.0d0
        scr(20)=46189.0d0/262144.0d0
        scr(22)=88179.0d0/524288.0d0
        scr(24)=676039.0d0/4194304.0d0
        scr(26)=1300075.0d0/8388608.0d0
        scr(28)=5014575.0d0/33554432.0d0
        scr(30)=9694845.0d0/67108864.0d0
        scr(32)=300540195.0d0/2147483648.0d0
        scr(34)=583401555.0d0/4294967296.0d0
        scr(36)=2268783825.0d0/17179869184.0d0
        scr(38)=4418157975.0d0/34359738368.0d0
        scr(40)=34461632205.0d0/274877906944.0d0
        scr(42)=67282234305.0d0/549755813888.0d0
        scr(44)=263012370465.0d0/2199023255552.0d0
        scr(46)=514589420475.0d0/4398046511104.0d0
        scr(48)=8061900920775.0d0/70368744177664.0d0
        do 30 i=3,49,2
          scr(i)=scr(i-1)
  30    continue
        if (dkhorder.ge.49) then
          write (stdout,1003)
1003      format (//2X,'Note:  McWeeny param. works so far only up ',
     *            'to dkhorder = 48.',/2X)
          Call Abend
        endif
        goto 700
      endif
      if (paramtype.eq.'OPT') then
        scr(3)=0.25d0*(2.d0-sqrt(2.d0))
        scr(5)=(24.d0-17.d0*sqrt(2.d0))/64.d0
        scr(7)=3.d0*(256.0d0-181.d0*sqrt(2.d0))/2048.0d0
        scr(9)=(12288.0d0-8689.d0*sqrt(2.d0))/32768.0d0
        scr(11)=(786432.d0-556091.d0*sqrt(2.d0))/2097152.0d0
        scr(13)=3.d0*(4194304.d0-2965821.d0*sqrt(2.d0))/33554432.d0
        scr(15)=3.0d0*(134217728.d0-94906265.d0*sqrt(2.d0))
     *             /1073741824.d0
        scr(17)=3.d0*(2147483648.d0-1518500251.d0*sqrt(2.d0))
     *             /17179869184.d0
        scr(19)=3.d0*(274877906944.d0-194368031985.d0*sqrt(2.d0))
     *             /2199023255552.d0
        scr(2) = 0.5d0
        scr(4) = scr(3) - 0.5d0*(scr(2)**2)
        scr(6) = scr(5) - scr(2)*scr(4) + 0.5d0*(scr(3)**2)
        scr(8) = scr(7) - scr(2)*scr(6) + scr(3)*scr(5)
     *             - 0.5d0*(scr(4)**2)
        scr(10)= scr(9) - scr(2)*scr(8) - scr(4)*scr(6) + scr(3)*scr(7)
     *             + 0.5d0*(scr(5)**2)
        scr(12)= scr(11) - scr(2)*scr(10) + scr(3)*scr(9) -scr(4)*scr(8)
     *             + scr(5)*scr(7) - 0.5d0*(scr(6)**2)

        scr(14)= scr(13) -scr(2)*scr(12) +scr(3)*scr(11) -scr(4)*scr(10)
     *             + scr(5)*scr(9) - scr(6)*scr(8) + 0.5d0*(scr(7)**2)
        scr(16)= scr(15) -scr(2)*scr(14) +scr(3)*scr(13) -scr(4)*scr(12)
     *             + scr(5)*scr(11) -scr(6)*scr(10) +scr(7)*scr(9)
     *             - 0.5d0*(scr(8)**2)
        scr(18)= scr(17) -scr(2)*scr(16) +scr(3)*scr(15) -scr(4)*scr(14)
     *             + scr(5)*scr(13) - scr(6)*scr(12) + scr(7)*scr(11)
     *             - scr(8)*scr(10) + 0.5d0*(scr(9)**2)
        scr(20)= scr(19) -scr(2)*scr(18) +scr(3)*scr(17) -scr(4)*scr(16)
     *             + scr(5)*scr(15) - scr(6)*scr(14) + scr(7)*scr(13)
     *             - scr(8)*scr(12) + scr(9)*scr(11) -0.5d0*(scr(10)**2)
        scr(22)= scr(21) -scr(2)*scr(20) +scr(3)*scr(19) -scr(4)*scr(18)
     *             + scr(5)*scr(17) - scr(6)*scr(16) + scr(7)*scr(15)
     *             - scr(8)*scr(14) + scr(9)*scr(13) - scr(10)*scr(12)
     *             + 0.5d0*(scr(11)**2)
        scr(24)= scr(23) -scr(2)*scr(22) +scr(3)*scr(21) -scr(4)*scr(20)
     *             + scr(5)*scr(19) - scr(6)*scr(18) + scr(7)*scr(17)
     *             - scr(8)*scr(16) + scr(9)*scr(15) - scr(10)*scr(14)
     *             + scr(11)*scr(13) - 0.5d0*(scr(12)**2)
c
        if (dkhorder.gt.20) then
          write (stdout,1004)
1004      format (//2X,'Note:  U_opt param. works so far only up ',
     *            'to dkhorder = 20.',/2X)
          Call Abend
        endif
      endif
c
 700  continue
      do 80 i=1,maxorder
        ducoeffs(i)=scr(i)
  80  continue
c
      return
      end
