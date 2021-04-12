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
        subroutine UrobChV (L2,NaGrp,NbeGrp,LunAux)
c
c       vyraba fily so simulovanymi Cholesky vektormi
c       so spravnou strukturou
c
        implicit none
#include "chcc1.fh"
c@@     include 'o2v4.fh' stare
#include "o3v3.fh"
#include "chcc_files.fh"
c
        integer NaGrp,NbeGrp,LunAux
        real*8 L2(1)
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
          len=nc*DimGrpv(aGrp)*DimGrpv(beGrp)
c
c1.2      full L2 with random numbers
          schem=1.0d-2
          call RNFill (len,L2(1),schem)
c
c1.3      open proper file
*         open (unit=LunAux,file=L2Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,L2Name(aGrp,beGrp))
c
c1.4      write L2 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,L2(1))
c
        end do
        end do
c
c
        return
        end
