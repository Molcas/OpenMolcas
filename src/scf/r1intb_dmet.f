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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2016,2017, Roland Lindh                                *
************************************************************************
      Subroutine R1IntB_DMET
************************************************************************
*                                                                      *
*     purpose: Read basis set informations and one electron integrals  *
*              were not needed so far.                                 *
*                                                                      *
*     called from: PrFin                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
*---- Define local variables
      Character*8 Label
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Allocate memory for kinetic energy, mass velocity and darvin
*     integrals
*
      Call mma_allocate(KntE,nBT+4,Label='KntE')
      Call mma_allocate(MssVlc,nBT+4,Label='MssVlc')
      Call mma_allocate(Darwin,nBT+4,Label='Darwin')
*
*---- Read kinetic energy integrals
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='Kinetic '
      Write (6,*) "r1intb0"
#define _DMET_
#ifdef _DMET_
      write(6,*) "comp",comp
      Call PrMtrx(label,1,1,1,KntE)
#endif
*
      Call RdOne(iRc,iOpt,Label,iComp,KntE,iSyLbl)
      If (iRc.ne.0) Go To 777
      Write (6,*) "r1intb1"
      If (iRc.ne.0) Then
        Write (6,*) "r1intb2"
         Write (6,*) 'R1Intb: Error readin ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call QTrace
         Call Abend()
      End If
*
*---- Read mass velocity integrals
      lRel=.False.
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='MassVel '
      Call RdOne(iRc,iOpt,Label,iComp,MssVlc,iSyLbl)
      If (iRc.ne.0) Go To 777
*
*---- Read Darvin integrals
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='Darwin  '
      Call RdOne(iRc,iOpt,Label,iComp,Darwin,iSyLbl)
      If ( iRc.ne.0 ) Go To 777
      lRel=.True.
*
 777  Continue
      If (.Not.lRel) Then
         Call mma_deallocate(MssVlc)
         Call mma_deallocate(Darwin)
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
