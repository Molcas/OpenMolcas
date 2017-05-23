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
        subroutine ccsdt(ireturn)
        integer ireturn
        logical run_triples
        call reorg(run_triples,ireturn)
        Call Disable_Spool()
        call ccsd(ireturn,run_triples)
        if(run_triples) call cct3(ireturn)
        return
        end
