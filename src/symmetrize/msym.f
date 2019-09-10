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
* Copyright (C) 2015, Marcus Johansson                                 *
************************************************************************
      Subroutine fmsym(iReturn)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#ifdef _MSYM_
      Call fmsym_create_context(ctx)
      Call fmsym_set_elements(ctx)
      Call fmsym_find_symmetry(ctx)
      Call fmsym_symmetrize_molecule(ctx)
      Call fmsym_generate_orbital_subspaces(ctx)
      Call fmsym_symmetrize_orb_file(ctx,'INPORB')
      Call fmsym_release_context(ctx)
#endif
      iReturn = 0
      Return
      End
