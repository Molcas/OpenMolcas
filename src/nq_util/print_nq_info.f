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
      Subroutine Print_NQ_Info(iSpin)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "nq_info.fh"
      Logical Reduce_Prt
      External Reduce_Prt
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*                                                                      *
************************************************************************
*                                                                      *
      If (iPL.ge.3) Then
         Call GAIGOP_SCAL(nTotGP,'+')
         Call GADGOP_SCAL(Flop,'+')
         Write (6,*)
         Write (6,'(6X,A,T52,F17.10)')
     &            'Integrated DFT Energy   ',Energy_integrated
         Write (6,'(6X,A,T56,G17.10)')
     &            'Integrated number of electrons',Dens_I
         If (Grad_I.ne.Zero)
     &   Write (6,'(6X,A,T56,G17.10)')
     &            'Integrated |grad|             ',Grad_I
         If (Tau_I .ne.Zero)
     &   Write (6,'(6X,A,T56,G17.10)')
     &            'Integrated tau                ',Tau_I
         Write (6,'(6X,A,T54,I13)')
     &            'Total number of prunned grid points  ',nTotGP
         Write (6,'(6X,A,T52,F17.1)')
     &            'Number of grid points per SO-integral  ',Flop
         Write (6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Add_Info('DFT_Energy',[Energy_integrated],1,6)
      Call Add_Info('NQ_Density',[Dens_I],1,8)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSpin)
      End
