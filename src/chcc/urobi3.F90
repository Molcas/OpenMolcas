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
        subroutine UrobI3 (I3,NaGrp,NbeGrp,LunAux)
!
!       vyraba fily so simulovanymi I3       vektormi
!       so spravnou strukturou
!
        implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
!
        integer NaGrp,NbeGrp,LunAux
        real*8 I3(1)
!
!       help variables
        integer aGrp,beGrp,len
        real*8 schem
!
!1      cycle over a,be Groups
!
        do aGrp=1,NaGrp
        do beGrp=1,NbeGrp
!
!1.1      def length
          if (aGrp.eq.beGrp) then
          len=no*(no+1)*DimGrpv(aGrp)*(DimGrpv(beGrp)+1)/4
          else
          len=no*(no+1)*DimGrpv(aGrp)*DimGrpv(beGrp)/2
          end if
!
!1.2      full I3 with random numbers
          schem=1.0d-2
          call RNFill (len,I3(1),schem)
!
!1.3      open proper file
!         open (unit=LunAux,file=I3Name(aGrp,beGrp),form='unformatted')
          Call MOLCAS_BinaryOpen_Vanilla(LunAux,I3Name(aGrp,beGrp))
!
!1.4      write I3 into proper file
          write (6,*) aGrp,beGrp,len
          call wri_chcc (LunAux,len,I3(1))
!
          close (LunAux)
!
        end do
        end do
!
!
        return
        end
