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
      Subroutine OFembed_dmat(Dens,nDens)

      use OFembed, only: Do_OFemb
      Implicit Real*8 (a-h,o-z)
      Real*8 Dens(nDens)
#include "stdalloc.fh"
      Character*16 NamRfil
      Real*8, Allocatable :: D_Var(:)

      If (.not.Do_OFemb) Return
!
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name

      Call mma_allocate(D_var,nDens,Label='D_var')
      Call get_dArray('D1aoVar',D_var,nDens)
      Call daxpy_(nDens,One,D_var,1,Dens,1)
      Call mma_deallocate(D_Var)
!
      Call NameRun(NamRfil)   ! switch back to old RUNFILE name
!
      Return
      End
