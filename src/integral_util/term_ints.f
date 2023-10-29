!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1998, Roland Lindh                                     *
!***********************************************************************
      Subroutine Term_Ints()
!***********************************************************************
!                                                                      *
!     Object: to deallocate memory in association with two-electron    *
!             calculations.                                            *
!                                                                      *
!     Author: Roland Lindh, Chemical Physics, University of Lund,      *
!             Sweden. January '98.                                     *
!***********************************************************************
      use k2_arrays, only: FT, Aux, iSOSym, Destroy_BraKet_Base
      use iSD_Data, only: iCntr, iSh2Sh, iSO2Sh, iCntr, nShBf, iShOff
      use stdalloc, only: mma_deallocate
      Implicit None
!                                                                      *
!***********************************************************************
!                                                                      *
!     In case of semi-direct mode the memory is released externally.
!
      Call RlsMem_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Allocated(FT)) Call mma_deallocate(FT)
!
      Call  Destroy_Braket_base()

      If (Allocated(Aux)) Call mma_deallocate(Aux)
!
      If (Allocated(iSOSym)) Call mma_deallocate(iSOSym)
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Allocated(nShBf)) Then
         Call mma_deallocate(nShBF)
         Call mma_deallocate(iShOff)
         Call mma_deallocate(iSh2Sh)
         Call mma_deallocate(iSO2Sh)
         Call mma_deallocate(iCntr)
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Free memory for K2 data
!
      Call FreeK2()
!                                                                      *
!***********************************************************************
!                                                                      *
      End Subroutine Term_Ints
