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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1992, Per Ake Malmqvist                                *
************************************************************************
      SubRoutine SetUpR(nRys)
************************************************************************
*                                                                      *
* Object: to setup the coefficients for the Rys roots and weights.     *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Cehmistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to DaFile February '91                          *
*                                                                      *
*     Added: Call to READAB, P.-A. Malmqvist March 1992, to set up     *
*             the tables needed to calculate large-order roots and     *
*             weights on request.                                      *
************************************************************************
      use Her_RW
      use vRys_RW
      use Leg_RW
      implicit none
#include "real.fh"
#include "stdalloc.fh"
#include "status.fh"
#include "print.fh"
*
      integer :: nRys
*
      integer :: iRys, jRys
      integer :: MemHer, iHer, iOffR
*
      If (Allocated(iHerR2)) Then
         Call WarningMessage(2,
     &          'SetupR: Rys_Status is already active!')
         Call Abend()
      End If
*
#ifdef _RYS_SCRATCH_
      CALL SetAux(1.0D-16)
#endif
*
      CALL Read_ABData
*
      CALL Read_RysRW
*
*     Set up the square of roots and the weights for Hermite polynomials
*     We will only do this for the even numbered polynomials.
*
      MemHer=nRys*(nRys+1)/2
      Call mma_allocate(iHerR2,nRys,label='iHerR2')
      iHerR2(1)=1
      Call mma_allocate(iHerW2,nRys,label='iHerW2')
      iHerW2(1)=1
      Call mma_allocate(HerR2,MemHer,label='HerR2')
      Call mma_allocate(HerW2,MemHer,label='HerW2')
*
      If (2*nRys.gt.MaxHer) Then
         Call WarningMessage(2,'SetupR: 2*nRys>MaxHer')
         Call Abend()
      End If
      Do 110 iRys=1,nRys
         iHer=2*iRys
         iOffR=(iRys*(iRys-1))/2
         iHerR2(iRys) = iHerR2(1) + iOffR
         iHerW2(iRys) = iHerW2(1) + iOffR
         Do 105 jRys=0,iRys-1
            HerR2(iHerR2(iRys)+jRys) = HerR(iHerR(iHer)+iRys+jRys)**2
            HerW2(iHerW2(iRys)+jRys) = HerW(iHerW(iHer)+iRys+jRys)
 105     Continue
 110  Continue
*
*define _DEBUG_
#ifdef _DEBUG_
      Call TriPrt(' Hermite squared roots',' ',HerR2(iHerR2(1)),nRys)
      Call TriPrt(' Hermite weights      ',' ',HerW2(iHerW2(1)),nRys)
#endif
*
      Return
      End
