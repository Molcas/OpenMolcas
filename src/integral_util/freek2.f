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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
      SubRoutine Freek2()
!***********************************************************************
!                                                                      *
!  Object: deallocate memory for pair entities                         *
!                                                                      *
! Called from: ClsSew                                                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, Sweden. November '92                            *
!***********************************************************************
      use stdalloc, only: mma_deallocate
      use k2_structure, only: free_k2data, Indk2, k2_Processed
      Implicit None
!
      If (.Not.Allocated(Indk2)) Return

!     Deallocate k2 entities
!
      Call Free_k2data()

      Call mma_deallocate(Indk2)
      k2_processed=.False.
!
      Return
      End SubRoutine Freek2
