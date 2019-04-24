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
      subroutine trunc_dkh (nbas,isize,znuc,maxexp,minexp,dkhorder,
     *                      clight,scrno1,scrno2,scrno3,vv,
     *                      nn,e,scr5,scr8,scr9,tran,sinv)
c
************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 25.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer nbas,isize,dkhorder,scrno1,scrno2,scrno3
      real*8 znuc,maxexp,minexp,clight,vv(nbas,nbas),
     *                 nn(nbas,nbas),e(nbas),
     *                 scr5(nbas,nbas,scrno1),scr8(isize,scrno2),
     *                 scr9(nbas,scrno3),tran(nbas,nbas),sinv(nbas,nbas)
c
      real*8 scaling,truncthrsh,ddum1,ddum2,maxew
      parameter (truncthrsh=1.0d-10)
c
C     integer i,j,ij
      integer i,j
      real*8 dum(1)
c
c
      write (stdout,1001) truncthrsh
1001  format (/8X,76('-'),/8X,'|',4X,'Determination of necessary value',
     *        ' of dkhorder for exact decoupling',5X,'|',/8X,76('-'),
     *        //4X,'Truncation threshold (truncthrsh):',1X,D13.6)
c
      scaling=1.0d0/(4.0d0*clight*clight)
c       scaling=1.0d0
c
      write (stdout,1020)
1020  format (/4X,'Calculation of scaled eigenvalues ( scaling=1',
     *        '/(4c^2) ) of',//8X,'the truncation estimate operator',
     *        ' V_k = AVA*(AV~A)^(k-1):',//8X,'k',26X,'Smallest',20X,
     *        'Largest',14X,'Abs_value',/)
c
      call mat_triang (scr8(1,1),nbas,vv)
c
      do 10 i=0,3*maxorder
        dum=1.0
        call diagr (scr8(1,1),nbas,tran,scr9(1,1),sinv,scr5(1,1,3),dum)
        do 20 j=1,nbas
          scr9(j,1)=scr9(j,1)*scaling
  20    continue
        ddum1=abs(scr9(1,1))
        ddum2=abs(scr9(nbas,1))
        if (ddum1.ge.ddum2) maxew=ddum1
        if (ddum2.ge.ddum1) maxew=ddum2
        write (stdout,1030) i+1,i+1,scr9(1,1),scr9(nbas,1),maxew
1030    format (7X,I2,6X,'DKH',I2,3X,F24.12,3X,F24.12,9X,D13.6)
        if (maxew.lt.truncthrsh) then
          write (stdout,1040) i,dkhorder
1040      format (/4X,'For this system (Z, nbas, maxexp, truncthrsh)',
     *            ' DKH',I2,/6X,'should be sufficient for exact ',
     *            'decoupling.',//4X,'Here:  dkhorder = ',I2,'.',
     *            //2X,76('-'))
          goto 150
        endif
        call mat_sq_from_t (scr5(1,1,1),nbas,scr8(1,1))
        call DGEMM_('N','N',nbas,nbas,nbas,1.0d0,scr5(1,1,1),
     *               nbas,nn,nbas,0.0d0,scr5(1,1,2),nbas)
        call mat_copy (scr5(1,1,1),nbas,nbas,scr5(1,1,2))
        call mat_triang (scr8(1,1),nbas,scr5(1,1,1))
c
  10  continue
c
 150  continue
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real(znuc)
        call Unused_real(maxexp)
        call Unused_real(minexp)
        call Unused_real_array(e)
      end if
      end
