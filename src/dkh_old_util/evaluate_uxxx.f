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
      subroutine evaluate_Uxxx (dkhorder,xorder,paramtype,dkhscfflg,
     *                    uused,unumber,utimes,uorder,umult,u,uscrchar,
     *                    uscrleng,utimestot,umulttot,totmult)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 14.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
      character*(3) paramtype
      logical dkhscfflg
c
      integer uused,unumber,utimes(maxunumber),uorder(maxunumber,3),
     *        umult(maxunumber),uscrleng(maxunumber),utimestot,umulttot,
     *        totmult
      character*(4) u(maxunumber)
      character*(3) uscrchar(maxunumber)
c
      integer i
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
c
c-------------------------------------------------------------------------
c
CMR      write (stdout,1007)
CMR1007  format (//5X,'|',70('-'),'|',/5X,'|',6X,'STEP 6: Determine ',
CMR     *        'auxiliary matrices Uxxx ',22X,'|',/5X,'|',
CMR     *        70('-'),'|')
c
      if (dbgflg.ge.1) then
        write (dbgunit,1009)
1009    format (/2X,'Step 6: Enter subroutine "evaluate_Uxxx" ',//2X,
     *          'Determine number of matrix multiplications necessary',
     *          ' for the evaluation of all Uxxx matrices,',/4X,
     *          ' which have been used so far.')
      endif
#ifndef MOLPRO
#ifdef _MOLCAS_
      DKHUnit5=5
      DKHUnit5=IsFreeUnit(DKHUnit5)
      Call Molcas_Open(DKHUnit5,'dkhops.15')
#else
      open (dkhunit5, file='dkhops.15', status='UNKNOWN',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit5)
      write (dkhunit5,1021) dkhorder,paramtype,xorder,dkhscfflg
1021  format (50('-'),/2X,'dkhorder = ',I2,10X,A3,/2X,'xorder   = ',I2,
     *        /2X,'dkhscfflg = ',L1,/'+++',47('-'))
      write (dkhunit5,1042) unumber
1042  format (I3,17X,'order(V)  order(X)  order(tot)')
      do 20 i=1,uused
        if (utimes(i).ge.1) then
          write (dkhunit5,1043) uscrleng(i),u(i),uscrchar(i),
     *                          uorder(i,1),uorder(i,2),uorder(i,3)
1043      format (I1,3X,A4,3X,A3,8X,I2,8X,I2,9X,I2)
        endif
  20  continue
c
#ifdef _MOLCAS_
      close (dkhunit5)
#endif
c
      totmult=totmult+umulttot
      if (dbgflg.ge.1) then
        write (dbgunit,2766) unumber,unumber,umulttot,dkhorder,xorder,
     *                       totmult
2766    format (/2X,'All ',I3,' Uxxx matrices have been written to the',
     *          ' file "dkhops.15".',//2X,'Number of matrix ',
     *          'multiplications necessary to evaluate all ',I3,' Uxxx',
     *          ' expressions:',2X,I5,//2X,'Total number of matrix',
     *          ' multiplications needed for evaluation of all DKH(',
     *          I2,',',I2,')-corrections:',3X,I8)
        write (dbgunit,2769)
2769    format (/2X,'Subroutine evaluate_Uxxx (step6) has now been ',
     *          'completed.',//2X,145('*'))
      endif
CMR      write (stdout,1034) unumber,unumber,umulttot,dkhorder,xorder,
CMR     *                    totmult
CMR1034  format (/2X,'All ',I3,' required Uxxx matrices have been ',
CMR     *        'written to the file "dkhops.15".',
CMR     *        //2X,'Total number of low-level matrix multiplications',
CMR     *        /4X,'necessary for the evaluation of all ',I3,' Uxxx:',
CMR     *        11X,I10,//2X,'Total number of matrix multiplications',
CMR     *        ' for the',/4X,'evaluation of all terms contributing',
CMR     *        ' to DKH(',I2,',',I2,'):',4X,I11,/2X)
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(umult)
        call Unused_integer(utimestot)
      end if
      end
