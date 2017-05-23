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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This program generates generally contracted basis sets from multiple *
* atomic densities. The densities are averaged with weights and the    *
* resulting eigenvectors, Natural Orbitals (NO), are used as basis     *
* functions where the corresponding eigenvalues, occupation numbers,   *
* are used as measure of importance.                                   *
* It is also possible to project out the contributions of vectors      *
* from the averaged density before generating the NO's.                *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine GenAno(ireturn)
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
#include "symlab.fh"
      Call MkType
      Call InpCtl_GENANO()
      Do 100 i=1,nSets
         Call RdCmo
         Call UpDens
         If(isUHF.eq.1) Then
            Call RdCmo()
            Call UpDens()
            isUHF=0
         End If
100   Continue
      Call SphAve()
      If(iProj.eq.1) Call Proj1()
      If(iProj.eq.2) Call Proj2()
      Call MkAno()
*
      ireturn=0
      return
      End
