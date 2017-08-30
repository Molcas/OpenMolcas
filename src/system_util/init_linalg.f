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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine Init_LinAlg()
* Initialization procedure for linear algebra libraries (if needed)
      Implicit None
#ifdef _DELAYED_
      Integer :: iPr
      Character (Len=8) ::  linalg_info
      Character (Len=1024) :: linalg_lib

      linalg_lib='Internal'
      Call GetEnvf('MOLCAS_LINALG',linalg_lib)
      Call GetEnvf('MOLCAS_LINALG_INFO',linalg_info)
      iPr = 0
      Call UpCase(linalg_info)
      linalg_info=AdjustL(linalg_info)
      If ((linalg_info(1:1) .ne. ' ') .and.
     &    (linalg_info(1:3) .ne. 'NO ') .and.
     &    (linalg_info(1:4) .ne. 'OFF ') .and.
     &    (linalg_info(1:1) .ne. '0')) iPr = 1
      Call Initialize_BLAS(linalg_lib,iPr)
#endif
      End Subroutine Init_LinAlg
