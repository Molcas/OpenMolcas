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
* Copyright (C) 1997, Per Ake Malmqvist                                *
*               2018, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Print_CI_Mix(EigVec)
      Use RefWfn, Only: refwfn_active, refwfn_is_h5, refwfn_id,
     &                  refwfn_filename, refwfn_close, iadr15
#ifdef _HDF5_
      Use mh5, Only: mh5_open_file_r, mh5_fetch_dset_array_real
#endif
      Implicit None
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
      Real*8 :: EigVec(nState,nState)
      Integer :: iState, iiState, jSNum, iDisk
      Real*8, Allocatable, Dimension(:) :: cCI, mCI
      Logical :: Close_refwfn


      Call mma_allocate(mCI, nConf, Label='MixCICoeff')
      Call mma_allocate(cCI, nConf, Label='CICoeff')

      Close_refwfn = .False.
      If (.Not.refwfn_active) Then
        ! bypass refwfn_open, because we don't want to set global stuff
        If (refwfn_is_h5) Then
#ifdef _HDF5_
          refwfn_id = mh5_open_file_r(refwfn_filename)
#else
* This should never happen
          Call AbEnd()
#endif
        Else
          refwfn_id=15
          Call DAName(refwfn_id,refwfn_filename)
        End If
        Close_refwfn = .True.
      End If

      Call CollapseOutput(1,'Mixed CI coefficients:')

      Write(6,*)
      Write(6,*)' The original CI arrays are now mixed as linear'
      Write(6,*)' combinations, given by the eigenvectors.'
      Write(6,*)

      Do iState=1,nState
        Call FZero(mCI, nConf)
        iDisk=iAdr15(4)
        Do iiState=1,nState
          jSNum=mState(iiState)
          If (refwfn_is_h5) Then
#ifdef _HDF5_
            Call mh5_fetch_dset_array_real(
     &           refwfn_id,'CI_VECTORS',cCI,[nConf,1],[0,jSNum-1])
#else
* This should never happen
            Call AbEnd()
#endif
          Else
            Call dDAFile(refwfn_id,2,cCI,nConf,iDisk)
          End If
          Call daXpY_(nConf,EigVec(iiState,iState),cCI,1,mCI,1)
        End Do
        Write(6,'(1X,A,I3)')
     &     ' The CI coefficients for the MIXED state nr. ',iState
        Call PrWf_CP2(stSym,nConf,mCI,CITHR)
      End Do

      Call CollapseOutput(0,'Mixed CI coefficients:')
      Write(6,*)

      If (Close_refwfn) Call refwfn_close()

      Call mma_deallocate(mCI)
      Call mma_deallocate(cCI)


      End Subroutine Print_CI_Mix
