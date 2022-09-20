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

subroutine WRIDIR(G,IDIMENS,IFILE,IAS)

implicit none
integer IDIMENS, IFILE, IAS
real*8 G(IDIMENS)

write(IFILE,rec=IAS) G

return

entry READIR(G,IDIMENS,IFILE,IAS)

read(IFILE,rec=IAS) G

return

end subroutine WRIDIR
