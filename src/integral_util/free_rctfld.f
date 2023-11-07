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
      Subroutine Free_RctFld()
      use PCM_arrays, only: Centr, dCntr, dPnt, dRad, dTes, IntSph,
     &                      MM, NewSph, nVert, PCMDm, PCMiSph, PCMTess,
     &                      PCMSph, PCM_N, PCM_SQ, Vert, SSph
      use Langevin_arrays, only: Cavxyz, Davxyz, Ravxyz, dField, Dip,
     &                           DipEF, Field, Grid, PolEF
      use stdalloc, only: mma_deallocate
      use rctfld_module, only: lLangevin, PCM, DoDeriv
      use External_Centers, only: iXPolType
      Implicit None
!
      If (.NOT.Allocated(MM)) Return

      Call mma_deallocate(MM)

      If (lLangevin .or. (iXPolType.gt.0)) Then
         Call mma_deallocate(Field)
         Call mma_deallocate(dField)
         Call mma_deallocate(Dip)
         Call mma_deallocate(PolEf)
         Call mma_deallocate(DipEf)
         Call mma_deallocate(Grid)
         Call mma_deallocate(davxyz)
         Call mma_deallocate(cavxyz)
         Call mma_deallocate(ravxyz)
      End If

      If (PCM) Then
         Call mma_deallocate(NewSph)
         Call mma_deallocate(IntSph)
         Call mma_deallocate(NVert)
         Call mma_deallocate(PCMiSph)
         Call mma_deallocate(PCM_N)
         Call mma_deallocate(PCMDM)
         Call mma_deallocate(SSph)
         Call mma_deallocate(Centr)
         Call mma_deallocate(Vert)
         Call mma_deallocate(PCMTess)
         Call mma_deallocate(PCMSph)
!
!---- Free the space for geometric derivatives
!
         If (DoDeriv) Then
            Call mma_deallocate(dTes)
            Call mma_deallocate(dPnt)
            Call mma_deallocate(dRad)
            Call mma_deallocate(dCntr)
            Call mma_deallocate(PCM_SQ)
         End If
!
      End If
!
      Return
      End Subroutine Free_RctFld
