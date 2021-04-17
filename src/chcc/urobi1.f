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
        subroutine UrobI1 (I1,NaGrp,LunAux)
c
c       vyraba fily so simulovanymi I1       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,LunAux
        real*8 I1(1)
c
c       help variables
        integer aGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
c
c1.1      def length
          len=no*DimGrpv(aGrp)*no*(no+1)/2
c
c1.2      full I1 with random numbers
          schem=1.0d-2
          call RNFill (len,I1(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=I1Name(aGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I1Name(aGrp))
c
c1.4      write I1 into proper file
          write (6,*) aGrp,len
          call wri_chcc (LunAux,len,I1(1))
c
        close (LunAux)
c
        end do
c
c
        return
        end
