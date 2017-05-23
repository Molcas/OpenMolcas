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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine sets the name for the runfile.                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine NameRun(Name)
#include "runinfo.fh"
      Character*(*) Name
*     Integer       iRc
*     Integer       iOpt

      If(Name.eq.'#Pop') Then
         RunName=RnNmStk(1)
         RnNmStk(1)=RnNmStk(2)
         RnNmStk(2)=RnNmStk(3)
         RnNmStk(3)=RnNmStk(4)
      Else
         RnNmStk(4)=RnNmStk(3)
         RnNmStk(3)=RnNmStk(2)
         RnNmStk(2)=RnNmStk(1)
         RnNmStk(1)=RunName
         RunName=Name
      End If

      Call ClrRunCache()
*
* Do not create the run file when naming it, patch 6.7.263
*
*     iRc=0
*     iOpt=1
*     Call MkRun(iRc,iOpt)

      Return
      End

      Subroutine ClrRunCache

      Call ClrRunCacheDS()
      Call ClrRunCacheIS()
      Return

      End

      Subroutine ClrRunCacheDS()
      Implicit None
#include "pg_ds_info.fh"

      Integer  i

      do i=1,num_DS_init
         i_DS_inmem(i)=0.d0
         DS_init(i)=0
         iLbl_DS_inmem(i)=' '
      end do
      num_DS_init=0

      End

      Subroutine ClrRunCacheIS()
      Implicit None
#include "pg_is_info.fh"

      Integer  i

      do i=1,num_IS_init
         i_IS_inmem(i)=0
         IS_init(i)=0
         iLbl_IS_inmem(i)=' '
      end do
      num_IS_init=0

      End
