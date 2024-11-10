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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine Fmat_m(D1A,FI,FA)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Update the Fock matrix for the active orbitals and transform     *
*     it to MO basis as well as the matrix FI (Fock matrix) for        *
*     frozen and inactive orbitals).                                   *
*                                                                      *
*     calling arguments:                                               *
*     D1A     : array of real*8                                        *
*               active one body density matrix in AO-basis             *
*     FI      : array of real*8                                        *
*               inactive Fock matrix                                   *
*     FA      : array of real*8                                        *
*               active Fock matrix                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     B.O. Roos and P.Aa. Malmqvist                                    *
*     University of Lund, Sweden, 1989                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     - updated for MOLCAS version 2                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1991               *
*     - updated for MOLCAS version 3                                   *
*       M.P. Fuelscher, University of Lund, Sweden, 1993               *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: ECAS, EMY, VIA

      Implicit None

      Real*8 D1A(*) , FI(*) , FA(*)

#include "rasdim.fh"
#include "general.fh"
      Character(LEN=16), Parameter:: ROUTINE='FMAT    '


      Real*8, allocatable:: Tmp1(:)
      Integer iPrLev
      Real*8, External:: DDot_
      real*8 vaa

C Local print level (if any)
      IPRLEV=IPRLOC(4)

!************************************************************
! Here we should start the real work!
!************************************************************
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
!     create FA in AO basis
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
      Call Fold(nSym,nBas,D1A,Tmp1)

      ! Active-Active contribution to ECAS
      VAA = 0.5D0*ddot_(nTot1,FA,1,Tmp1,1)

!     Inactive-active contribution to ECAS
      VIA=dDot_(nTot1,FI,1,Tmp1,1)
      ECAS = EMY + VIA + VAA
      If ( iPrLev.ge.DEBUG ) then
        Write(LF,'(A,ES20.10)') ' Total core energy:            ',EMY
        Write(LF,'(A,ES20.10)') ' inactive-active interaction:  ',VIA
        Write(LF,'(A,ES20.10)') ' active-active interaction:  ',VAA
        Write(LF,'(A,ES20.10)') ' CAS energy (core+interaction):',ECAS
      End If
      Call mma_deallocate(Tmp1)

      End Subroutine Fmat_m
