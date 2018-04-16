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
* Copyright (C) 1989-1992,1998,1999, Roland Lindh                      *
*               1990, IBM                                              *
*               2000-2015, Valera Veryazov                             *
************************************************************************
      subroutine Grid_it_nosupport(iRun,INPORB,ireturn)
c  iRun =1 normal run, 0=trancated from scf
************************************************************************
*                                                                      *
*  Object: Driver for evaluation MO values on a grid.                  *
*                                                                      *
* Called from: None                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              XuFlow (IBM)                                            *
*              Seward_init                                             *
*             *SetUp0                                                  *
*              GetMem                                                  *
*              GetInf                                                  *
*              DrvMO                                                   *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
*          July '89 - May '90                                          *
*                                                                      *
*          Modified to gradient calculations September 1991 -          *
*          February 1992.                                              *
*                                                                      *
*          Modified to evaluating the spin density and spin density    *
*          gradients on a grid                                         *
*                                                                      *
*          Modified to interface with the MSI Cerius 2.                *
*          April 1998                                                  *
*                                                                      *
*          Modified at June- Sept 1999                                 *
*                                                                      *
*   This code was rewritten from scratch by V. Veryazov                *
*          Lund, 2000-2015                                             *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
c      Character*120 Lines(17)
      Character INPORB*(*)
      Logical OldTst, DoRys
#include "grid.nosupport.fh"
#include "warnings.fh"
*
*     Prologue
*
      Call qEnter('GRID')
      levelprint=IPRINTLEVEL(-1)
      if(iRun.eq.0.and.levelprint.lt.3) then
        levelprint=0
        levelprint=IPRINTLEVEL(levelprint)
      endif
      if(iRun.eq.1) then
      iRout=1
      Call SetTim
*
c      Call bXML('GRID_IT')
*
*     Get the work space size
*
      endif
*
c      Call Seward_init
c*
c*     Get the input information as Seward dumped on INFO.
c*
c      nDiff=0
c      DoRys=.False.
c      Call GetInf(Info,nInfo,DoRys,nDiff,idum)
      nDiff=0
      DoRys=.False.
      Call IniSew(Info,DoRys,nDiff)

      OldTst = Test
*
*
*---- Read the input
*
      iReturn=0
      Call Input_Grid_It_nosupport(iRun,INPORB,iReturn)
      if (iReturn.eq._RC_INVOKED_OTHER_MODULE_) then
c* take care to close files and release the potential memory...
c       close(unit=LuOrb)
        close(unit=LuVal)
        if(isUHF.eq.1) close(unit=LuVal_ab)
      goto 999
      endif
*
*
*     Start computing the spin density and spin density gradient
*     at the grid.
*
      Call DrvMO_nosupport(iRun,INPORB)
*
*-----At the end of the calculation free all memory to check for
*     corruption of the memory.
*
999   continue

      Call ClsSew()
*
*     Epilogue
*
      if(iRun.eq.1) then
        Call qStat(' ')
c        Call eXML('GRID_IT')
c      else
c      write(6,*) 'Input file for molcasgv was generated'
      endif
c      ireturn=0


      Call qExit('GRID')
      return
      End
