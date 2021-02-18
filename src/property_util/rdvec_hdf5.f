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
#ifdef _HDF5_
      Use mh5, Only: mh5_exists_dset, mh5_fetch_dset
#endif
      Implicit None
      Integer, Intent(In) :: fileid,nSym,nBas(nSym)
      Character(Len=*), Intent(In) :: Label
      Real*8, Dimension(*) :: CMO,Occ,Ene
      Integer, Dimension(*) :: Ind
#ifdef _HDF5_
#include "stdalloc.fh"
      Integer :: Beta,nB
      Character(Len=128) :: DataSet,su,sl
      Character(Len=1), Allocatable :: typestring(:)

      Beta=0
      su=''
      sl=''
      If (Index(Label,'A').gt.0) Then
        Beta=-1
        su='ALPHA_'
        sl='alpha '
      End If
      If (Index(Label,'B').gt.0) Then
        If (Beta.ne.0) Then
          Write(6,*)
          Call AbEnd
        End If
        Beta=1
        su='BETA_'
        sl='beta '
      End If

      If (Index(Label,'E').gt.0) Then
        DataSet='MO_'//Trim(su)//'ENERGIES'
        If (mh5_exists_dset(fileid,DataSet)) Then
          Call mh5_fetch_dset(fileid,DataSet,Ene)
        Else
          Write(6,*) 'The HDF5 file does not contain '//
     &               Trim(sl)//'MO energies.'
          Call AbEnd()
        End If
      End If

      If (Index(Label,'O').gt.0) Then
        DataSet='MO_'//Trim(su)//'OCCUPATIONS'
        If (mh5_exists_dset(fileid,DataSet)) Then
          Call mh5_fetch_dset(fileid,DataSet,Occ)
        Else
          Write(6,*) 'The HDF5 file does not contain '//
     &               Trim(sl)//'MO occupations.'
          Call AbEnd()
        End If
      End If

      If (Index(Label,'C').gt.0) Then
        DataSet='MO_'//Trim(su)//'VECTORS'
        If (mh5_exists_dset(fileid,DataSet)) Then
          Call mh5_fetch_dset(fileid,DataSet,CMO)
        Else
          Write(6,*) 'The HDF5 file does not contain '//
     &               Trim(sl)//'MO coefficients.'
          Call AbEnd()
        End If
      End If

      If (Index(Label,'I').gt.0) Then
        nB = Sum(nBas)
        Call mma_allocate(typestring,nB)
        DataSet='MO_'//Trim(su)//'TYPEINDICES'
        If (mh5_exists_dset(fileid,DataSet)) Then
          Call mh5_fetch_dset(fileid,DataSet,typestring)
          Call tpstr2tpidx(typestring,Ind,nB)
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
         Call Unused_integer(nSym)
         Call Unused_integer_array(nBas)
         Call Unused_real_array(CMO)
         Call Unused_real_array(Occ)
         Call Unused_real_array(Ene)
         Call Unused_integer_array(Ind)
      End If
#endif
      End Subroutine RdVec_HDF5
