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
* Copyright (C) 2001-2016, Valera Veryazov                             *
************************************************************************
      subroutine finish(rc)
C     Gracefully shuts down a program module.
C     After everything is closed properly, xquit is
C     called to do the actual termination.
      implicit none
      integer :: rc
#include "WrkSpc.fh"
#include "timtra.fh"
      integer :: idum = 0
      integer :: iwarn

      if(nfld_tim.gt.0) Call GetMem('iGATim','Free','Real',
     &                  iGATim,iDum)
      if(nfld_stat.gt.0) Call GetMem('iGAStat','Free','Real',
     &                   iGAStat,iDum)

      Call GetMem('ip_iDum','Free','Inte',ip_iDummy,1)
      Call GetMem('ip_sDum','Free','SNGL',ip_sDummy, 1)
      Call GetMem('ip_Dum', 'Free','Real',ip_Dummy, 1)
      Call GetMem('Finish','List','Real',iDum,iDum)
      Call GetMem('Finish','Term','Real',iDum,iDum)

      Call fin_run_use()

      Call StatusLine('Happy landing',' ')
      Call WarningCheckOut(iWarn)
      If (iWarn.gt.1) Then
        Call WarningMessage(1,
     &          'There were warnings during the execution;'//
     &          'Please, check the output with care!')
      End If

      Call prgmfree()
      Call AixCheck()
      call xml_close('module')

#ifdef _DELAYED_
      Call close_BLAS()
#endif

      Call xquit(rc)
      End
