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
        subroutine UrobI3 (I3,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi I3       vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 I3(1)
c
c       help variables
        integer aGrp,beGrp,len
        real*8 schem
c
c1      cycle over a,be Groups
c
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
c
c1.1      def length
          if (aGrp.eq.beGrp) then
          len=no*(no+1)*DimGrpv(aGrp)*(DimGrpv(beGrp)+1)/4
          else
          len=no*(no+1)*DimGrpv(aGrp)*DimGrpv(beGrp)/2
          end if
c
c1.2      full I3 with random numbers
          schem=1.0d-2
          call RNFill (len,I3(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=I3Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I3Name(aGrp,beGrp))
c
c1.4      write I3 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,I3(1))
c
          close (LunAux)
c
        end do
        end do
c
c
        return
        end
