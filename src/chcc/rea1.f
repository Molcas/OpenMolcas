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
        subroutine rea1 (lun,length,A)
c
c       nacitane bloku dat z Sekvencneho suboru, ak koniec tak
c       rewind a citanie znova od zaciatku
c       Toto je ozaj ditry
c
        implicit none
        integer lun,length
        real*8 A(1:length)
c
        read (lun,end=99) A
        return
c
99      rewind(lun)
        read (lun) A
c
        end
