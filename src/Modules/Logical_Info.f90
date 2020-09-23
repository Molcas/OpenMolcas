!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
!     Some incomplete documentation
!
!     Vlct  : logical flag to indicate that the velocity integrals will*
!             be computed                                              *
!     lRel  : logical flag to indicate that the 1-electron Darwin and  *
!             mass-velocity integrals should be computed               *
!     UnNorm: logical flag to indicate that the primitive basis        *
!             functions should not be normalized                       *
!     lSchw : logical flag to indicate that the Cauchy-Schwartz        *
!             inequality should be used in the prescreening            *
!     lAMFI : logical flag for Atomic Mean Field Integrals             *
!     GIAO: integrals for first derivative with respect to B           *
!     Do_FckInt: logical flag for FckInt                               *
!     Do_GuessOrb: logical flag for GuessOrb                           *
!     EMFR: Electromagnetic field radiation                            *
!     FNMC: finite nuclear mass correction
Module Logical_Info
Private
Public :: Vlct, &
          Logical_Info_Dmp, Logical_Info_Get

#include "stdalloc.fh"
Integer i
Logical :: Vlct=.True.
Logical :: lRel=.False.
Logical :: UnNorm=.False.
Logical :: lSchw=.True.
Logical :: lAMFI=.False.
Logical :: NEMO=.False.
Logical :: Do_GuessOrb=.True.
Logical :: Do_FckInt=.True.
Logical :: Align_Only=.False.
Logical :: DoFMM=.False.
Logical :: lRP=.False.
Logical :: lRP_Post=.False.
Logical :: Do_Align=.True.
Logical :: EMFR  =.False.
Logical :: GIAO=.False.
Logical :: lUPONLY=.False.
Logical :: lDOWNONLY=.False.
Logical :: FNMC=.False.
Logical :: lPSOI=.False.

Contains

Subroutine Logical_Info_Dmp()
  Real*8, Allocatable:: iDmp(:)
  Integer:: Len=19

  Call mma_allocate(iDmp,Len,Label='iDmp:Logical')

  i = 0
  If (Vlct) i = 1
  iDmp(01)= i
  i = 0
  If (lRel) i = 1
  iDmp(02)= i
  i = 0
  If (UnNorm) i = 1
  iDmp(03)= i
  i = 0
  If (lSchw) i = 1
  iDmp(04)= i
  i = 0
  If (lAMFI) i = 1
  iDmp(05)= i
  i = 0
  If (NEMO) i = 1
  iDmp(06)= i
  i = 0
  If (Do_GuessOrb) i = 1
  iDmp(07)= i
  i = 0
  If (Do_FckInt) i = 1
  iDmp(08)= i
  i = 0
  If (Align_Only) i = 1
  iDmp(09)= i
  i = 0
  If (DoFMM) i = 1
  iDmp(10)= i
  i = 0
  If (lRP) i = 1
  iDmp(11)= i
  i = 0
  If (lRP_Post) i = 1
  iDmp(12)= i
  i = 0
  If (Do_Align) i = 1
  iDmp(13)= i
  i = 0
  If (EMFR) i = 1
  iDmp(14)= i
  i = 0
  If (GIAO) i = 1
  iDmp(15)= i
  i = 0
  If (lUPONLY) i = 1
  iDmp(16)= i
  i = 0
  If (lDOWNONLY) i = 1
  iDmp(17)= i
  i = 0
  If (FNMC) i = 1
  iDmp(18)= i
  i = 0
  If (lPSOI) i = 1
  iDmp(10)= i

  Call Put_iArray('Logical_Info',iDmp,Len)
  Call mma_deallocate(iDmp)
End Subroutine Logical_Info_Dmp

Subroutine Logical_Info_Get()
  Real*8, Allocatable:: iDmp(:)
  Integer:: Len=19

  Call mma_allocate(iDmp,Len,Label='iDmp:Logical')
  Call Get_iArray('Logical_Info',iDmp,Len)

  Vlct         = iDmp(01).eq.1
  lRel         = iDmp(02).eq.1
  UnNorm       = iDmp(03).eq.1
  lSchw        = iDmp(04).eq.1
  lAMFI        = iDmp(05).eq.1
  NEMO         = iDmp(06).eq.1
  Do_GuessOrb  = iDmp(07).eq.1
  Do_FckInt    = iDmp(08).eq.1
  Align_Only   = iDmp(09).eq.1
  DoFMM        = iDmp(10).eq.1
  lRP          = iDmp(11).eq.1
  lRP_Post     = iDmp(12).eq.1
  Do_Align     = iDmp(13).eq.1
  EMFR         = iDmp(14).eq.1
  GIAO         = iDmp(15).eq.1
  lUPONLY      = iDmp(16).eq.1
  lDOWNONLY    = iDmp(17).eq.1
  FNMC         = iDmp(18).eq.1
  lPSOI        = iDmp(19).eq.1

  Call mma_deallocate(iDmp)


End Subroutine Logical_Info_Get

End Module Logical_Info
