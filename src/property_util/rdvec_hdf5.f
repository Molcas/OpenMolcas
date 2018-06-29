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
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
*  RdVec_HDF5
*
*> @brief
*>   Read orbital data from an HDF5 file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Similar to ::RdVec, but for HDF5 files. The \p Label argument can
*>  contain a `B` if one wants the data for beta orbitals.
*>
*> @param[in]  fileid   Identifier of the open HDF5 file
*> @param[in]  Label    Properties to read from the file
*> @param[in]  nSym     Number of irreps
*> @param[in]  nBas     Number of basis functions per irrep
*> @param[out] CMO      Orbital coefficients
*> @param[out] Occ      Orbital occupations
*> @param[out] Ene      Orbital energies
*> @param[out] Ind      Orbital type indices
************************************************************************
      Subroutine RdVec_HDF5(fileid,Label,nSym,nBas,CMO,Occ,Ene,Ind)
      Implicit None
      Integer, Intent(In) :: fileid,nSym,nBas(nSym)
      Character(Len=*), Intent(In) :: Label
      Real*8, Dimension(*) :: CMO,Occ,Ene
      Integer, Dimension(*) :: Ind
#ifdef _HDF5_
#include "mh5.fh"
#include "stdalloc.fh"
      Logical :: Beta
      Integer :: nB
      Character(1), Allocatable :: typestring(:)

      Beta = Index(Label,'B').gt.0

      If (Index(Label,'E').gt.0) Then
        If (Beta) Then
          If (mh5_exists_dset(fileid,'MO_BETA_ENERGIES')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_BETA_ENERGIES',Ene)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'beta MO energies.'
            Call AbEnd()
          End If
        Else
          If (mh5_exists_dset(fileid,'MO_ENERGIES')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_ENERGIES',Ene)
          Else If (mh5_exists_dset(fileid,'MO_ALPHA_ENERGIES')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_ALPHA_ENERGIES',Ene)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'MO energies.'
            Call AbEnd()
          End If
        End If
      End If

      If (Index(Label,'O').gt.0) Then
        If (Beta) Then
          If (mh5_exists_dset(fileid,'MO_BETA_OCCUPATIONS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_BETA_OCCUPATIONS',Occ)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'beta MO occupations.'
            Call AbEnd()
          End If
        Else
          If (mh5_exists_dset(fileid,'MO_OCCUPATIONS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_OCCUPATIONS',Occ)
          Else If (mh5_exists_dset(fileid,'MO_ALPHA_OCCUPATIONS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_ALPHA_OCCUPATIONS',Occ)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'MO occupations.'
            Call AbEnd()
          End If
        End If
      End If

      If (Index(Label,'C').gt.0) Then
        If (Beta) Then
          If (mh5_exists_dset(fileid,'MO_BETA_VECTORS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_BETA_VECTORS',CMO)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'beta MO coefficients.'
            Call AbEnd()
          End If
        Else
          If (mh5_exists_dset(fileid,'MO_VECTORS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_VECTORS',CMO)
          Else If (mh5_exists_dset(fileid,'MO_ALPHA_VECTORS')) Then
            Call mh5_fetch_dset_array_real(fileid,
     &           'MO_ALPHA_VECTORS',CMO)
          Else
            Write(6,'(1X,A)') 'The HDF5 file does not contain '//
     &                        'MO coefficients.'
            Call AbEnd()
          End If
        End If
      End If

      If (Index(Label,'I').gt.0) Then
        nB = Sum(nBas)
        Call mma_allocate(typestring,nB)
        If (Beta) Then
          If (mh5_exists_dset(fileid,'MO_BETA_TYPEINDICES')) Then
            Call mh5_fetch_dset_array_str(fileid,
     &           'MO_BETA_TYPEINDICES',typestring)
            Call tpstr2tpidx(typestring,Ind,nB)
          End If
        Else
          If (mh5_exists_dset(fileid,'MO_TYPEINDICES')) Then
            Call mh5_fetch_dset_array_str(fileid,
     &           'MO_TYPEINDICES',typestring)
            Call tpstr2tpidx(typestring,Ind,nB)
          Else If (mh5_exists_dset(fileid,'MO_ALPHA_TYPEINDICES')) Then
            Call mh5_fetch_dset_array_str(fileid,
     &           'MO_ALPHA_TYPEINDICES',typestring)
            Call tpstr2tpidx(typestring,Ind,nB)
          End If
        End If
        Call mma_deallocate(typestring)
      End If
#else
      Call WarningMessage(2,'Calling RdVec_HDF5, but HDF5 is disabled')
      Call AbEnd()
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(fileid)
         Call Unused_character(Label)
         Call Unused_integer_array(nSym)
         Call Unused_integer_array(nBas)
         Call Unused_real_array(CMO)
         Call Unused_real_array(Occ)
         Call Unused_real_array(Ene)
         Call Unused_integer_array(Ind)
      End If
#endif
      End Subroutine RdVec_HDF5
