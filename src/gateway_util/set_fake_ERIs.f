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
      Subroutine Set_Fake_ERIs
      use Basis_Info, only: nBas
      use RICD_Info, only: Do_RI, Cholesky
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "cholesky.fh"
      Character*16 NamRfil
      Integer, Dimension(:), Allocatable :: iSOShl
*
      write(6,*)
      write(6,*)'   *** Skipping anything related to ERIs ***'
      write(6,*)
*
      If (.not.Cholesky .and. .not.Do_RI) Return

      Call Get_NameRun(NamRfil)
      Call NameRun('AUXRFIL')
*
      Call Get_iScalar('ChoVec Address',CHO_ADRVEC)
      nBasT=nBas(0)
      Do i=1,nIrrep-1
         nBasT=nBasT+nBas(i)
      End Do
      Call mma_allocate(iSOShl,nBasT)
      Call Get_dScalar('Cholesky Threshold',THRCOM)
      Call Get_iArray('NumCho',NumCho,nIrrep)
      Call Get_iArray('iSOShl',ISOSHL,NBAST)
*
      Call NameRun(NamRfil)
      CALL Put_iArray('iSOShl',ISOSHL,NBAST)
      Call mma_deallocate(iSOShl)
      CALL Put_iArray('NumCho',NumCho,nIrrep)
      Call Put_iScalar('ChoVec Address',CHO_ADRVEC)
      CALL Put_dScalar('Cholesky Threshold',THRCOM)
*
      Return
      End
