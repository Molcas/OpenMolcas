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
      Subroutine Funi_Init
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "nq_info.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Call qEnter('Funi_Init')
*                                                                      *
************************************************************************
*                                                                      *
* Initialize the defaults values of the parameters.
*
*     Default Grid
      Quadrature='MHL '
      nR=75
      L_Quad=29
      Crowding=3.0D0
      Fade=6.0D0
      MBC=' '
*
      ntotgp=0
*
*     Various idefault thresholds for the integral evaluation.
*
      T_X      =1.0D-18
      T_Y      =1.0D-11
      Threshold=1.0D-25
*
      Angular_Prunning = On
      Grid_Type=Moving_Grid
      Rotational_Invariance= On
      NQ_Direct=On
      nGridMax=128
      Packing=On
C     Packing=Off
*
*     Bit 1: set Lobatto, not set Gauss and Gauss-Legendre
*     Bit 2: set scan the whole atomic grid, not set use subset
*     Bit 3: set Lebedev, override bit 1
      iOpt_Angular=4
*                                                                      *
************************************************************************
*                                                                      *
      Flop=Zero
      Do i = 0, LMax_NQ
         R_Max(i)=Zero
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('Funi_Init')
      Return
      End
