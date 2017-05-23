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
      SubRoutine PriTim(TimDat,nFld_par,nProcs)
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Real*8 TimDat(2*nFld_par,nProcs)
*
*     Call QEnter('PriTim')
      Write (6,'(1X,A)') 'Timing statistics of individual servers:'
      Write (6,'(1X,A5,6X,9A13)') ' node',
     &                           '   wait Dens.',
     &                           '     k2 stuff',
     &                           '     k4 stuff',
     &                           '  reduce Grad',
     &                           '    T O T A L',
     &                           '      # tasks',
     &                           '  # shl quads'
      Do inode = 1, nProcs
        Write (6,'(1X,I5,A6,5F13.2,2F13.0)')
     &           inode,' CPU  ',(TimDat(iFld,inode),iFld=1,nFld_par-1),
     &                           TimDat(nFld_par,inode),
     &                           TimDat(2*nFld_par,inode)
        Write (6,'(1X,5X,A6,9F13.2)')
     &                 ' Wall ',
     &                 (TimDat(nFld_par+iFld,inode),iFld=1,nFld_par-1)
      End Do
*     Call QExit('PriTim')
      Return
      End
