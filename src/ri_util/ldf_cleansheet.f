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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CleanSheet(nShell_Valence,nShell_Auxiliary)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Get a fresh start with Seward.
C
      Implicit None
      Integer nShell_Valence, nShell_Auxiliary
#include "status.fh"

      Logical Verbose, Free_K2
      Logical Indexation, DoFock, DoGrad
      Integer nShell_Tot
      Real*8  ThrAO

      Call Free_iSD()
      ! Essentially close seward
      Verbose=.False.
      Free_K2=.True.
      Call Term_Ints(Verbose,Free_K2)
      Call Free_iSD()

      ! Get number of valence shells
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      Indexation=.False.
      ThrAO=0.0d0
      DoFock=.False.
      DoGrad=.False.
      Call Setup_Ints(nShell_Valence,Indexation,ThrAO,DoFock,DoGrad)
      Verbose=.False.
      Free_K2=.True.
      Call Term_Ints(Verbose,Free_K2)
      Call Free_iSD()

      ! Set up seward for val+aux shells w. indexation
      Call Set_Basis_Mode('WithAuxiliary')
      Call Setup_iSD()
      Indexation=.True.
      ThrAO=0.0d0
      DoFock=.False.
      DoGrad=.False.
      Call Setup_Ints(nShell_Tot,Indexation,ThrAO,DoFock,DoGrad)

      ! Get number of auxiliary shells
      nShell_Auxiliary=nShell_Tot-1-nShell_Valence

      End
