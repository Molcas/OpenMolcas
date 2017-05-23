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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
      SUBROUTINE RESTART_CPFMCPF(C,NCONF)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "files_cpf.fh"
      DIMENSION C(*)
      IAD=0
      CALL dDAFILE(Lu_CI,2,C,NCONF,IAD)
      RETURN
      END
