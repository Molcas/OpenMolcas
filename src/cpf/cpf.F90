!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine CPF(IRETURN)
!***********************************************************************
!                                                                      *
! PER SIEGBAHN                                                         *
! MARGARETA BLOMBERG                                                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY                                  *
! UNIVERSITY OF LUND                                                   *
! SWEDEN                                                               *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                                C P F                                 *
! MODIFIED TO IBM BY ROLAND LINDH 02/17/88                             *
! MODIFIED TO MOLCAS-2 BY ROLAND LINDH 03/26/91                        *
! MODIFIED TO MOLCAS-3 BY M.P. FUELSCHER 08/31/93                      *
! MODIFIED TO MOLCAS-4 BY P.A. MALMQVIST AND N.W. MORTIARTY 10/25/96   *
! MODIFIED TO MOLCAS 4.1 BY R. LINDH 02/24/98 (Multi fileing)          *
! MODIFIED TO MODERN FORTRAN BY I. FDEZ. GALVAN 2022                   *
!***********************************************************************

use cpf_global, only: ICASE, INDX, ISAB, JSY, Lu_25, Lu_27, Lu_30, Lu_CI, Lu_CIGuga, Lu_CPFORB, Lu_TiABCD, Lu_TiABCI, Lu_TiABIJ, &
                      Lu_TraInt, Lu_TraOne
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: IRETURN
integer(kind=iwp) :: MEMORY

call mma_maxdble(MEMORY)
MEMORY = int(MEMORY*0.8_wp)

! Open files

Lu_CIGuga = 10
call DANAME(Lu_CIGuga,'CIGUGA')
Lu_TraInt = 50
call DANAME_MF(Lu_TraInt,'TRAINT')
Lu_TraOne = 17
call DANAME(Lu_TraOne,'TRAONE')
Lu_CI = 26
call DANAME(Lu_CI,'CPFVECT')
Lu_CPFORB = 19
! Temporaries:
Lu_TiABIJ = 60
call DANAME_MF(Lu_TiABIJ,'TIABIJ')
Lu_TiABCI = 70
call DANAME_MF(Lu_TiABCI,'TIABCI')
Lu_TiABCD = 80
call DANAME_MF(Lu_TiABCD,'TIABCD')
Lu_25 = 25
call DANAME(Lu_25,'FT25F001')
Lu_27 = 27
call DANAME(Lu_27,'FT27F001')
Lu_30 = 30
call DANAME(Lu_30,'FT30F001')

! Body

call SDCI_CPF(MEMORY)

! Deallocate the workspace

call mma_deallocate(ICASE)
call mma_deallocate(JSY)
call mma_deallocate(INDX)
call mma_deallocate(ISAB)

!                                                                      *
!***********************************************************************
!                                                                      *
! Close open dafiles

call DACLOS(Lu_CIGuga)
call DACLOS(Lu_TraInt)
call DACLOS(Lu_TraOne)
call DACLOS(Lu_CI)
call DACLOS(Lu_TiABIJ)
call DACLOS(Lu_TiABCI)
call DACLOS(Lu_TiABCD)
call DACLOS(Lu_25)
call DACLOS(Lu_27)
call DACLOS(Lu_30)
!                                                                      *
!***********************************************************************
!                                                                      *
call FASTIO('STATUS')
IRETURN = 0

return

end subroutine CPF
