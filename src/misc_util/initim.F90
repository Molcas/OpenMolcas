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
      SubRoutine IniTim()
      Use Para_Info, Only: nProcs
      Implicit None
#include "WrkSpc.fh"
#include "timtra.fh"
!
      if(nfld_tim.eq.0) return
      if(nfld_tim.gt.nfldmax) then
        Call WarningMessage(2,'Too many fields in IniTim')
        Write(6,*) 'nfld_tim:',nfld_tim
        call Abend()
      end if
      Call GetMem('iGATim','Allo','Real',iGATim,nProcs*nfld_tim*2)
      call Fzero(Work(iGATim),nProcs*nfld_tim*2)
      Return
      End
