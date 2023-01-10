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
      Subroutine read_mo (CMO,nfro,no,nv,ndel,nbas,nOrb)
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (A-H,O-Z)

*     declaration of calling arguments
      Type (DSBA_Type) CMO
      Integer lthCMO
      integer nfro_scf(8)
      integer nfro
#include "real.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: CMO_t(:,:)

#include "SysDef.fh"

*...  Read nSym, Energy, nBas, nOrb, nOcc, nFro, CMO and orbital energies from COMFILE
*
      Call Get_iArray('nFro',nFro_scf,1)
      If (nFro_scf(1).ne.0) Then
         Write (6,*) 'Some orbitals were frozen in SCF!'
         Call Abend()
      End If
c
      lthCMO=nBas*nBas
      Call mma_allocate(CMO_t,nBas,nBas,Label='CMO_t')
      Call Get_dArray_chk('Last orbitals',CMO_t,lthCMO)
c
c - transpose MO matrix, skip the frozen occupied orbitals
c
      call mo_transp(CMO%A0,CMO_t(:,1+nfro:nOrb),no,nv,ndel,nbas)
c
      Call mma_deallocate(CMO_t)
c
      Return
      End
