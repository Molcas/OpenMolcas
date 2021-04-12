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
        subroutine GetTau (Tau,T1,
     c                     aGrp,bGrp,dima,dimb,adda,addb,lunTau)
c
c       this routine do:
c       1) read the block of T2((a,b)',ij) from T2Name(aGrp,bGrp)
c        2) Make Tau ((a,b)',ij) in T2((a,b)',ij) array
c
c       I/O parameter description:
c       Tau      - array for Tau((a,b),ij) (O)
c       aGrp     - group of a index
c       bGrp     - group of b index
c       dima     - dimension group of a' index
c       dimb     - dimension group of b' index
c       adda     - shift of a' in full virtual space
c       addb     - shift of b' in full virtual space
c       lunTau   - Lun of opened file, where Tau is stored
c
c
        implicit none
#include "chcc1.fh"
#include "chcc_files.fh"
c
        real*8 Tau(1)
        real*8 T1(1)
        integer aGrp,bGrp,dima,dimb,adda,addb,lunTau
c
c       help variables
        integer length
c
c
c1      def legth
c
        if (aGrp.eq.bGrp) then
c       groups of a and b are equal, reading for a'>=b'
          length=no*no*Dima*(Dima+1)/2
        else
c       aGrp>bGrp, reading for a',b' in given groups
          length=no*no*Dima*Dimb
        end if
c
c
c2      read block of T2 amplitudes
c
        call GetX (Tau(1),length,LunTau,T2Name(aGrp,bGrp),1,1)
c
c
c3        make Tau
c
        if (aGrp.ne.bGrp) then
           call GetTauHlp1 (Tau(1),T1(1),dima,dimb,adda,addb,no,nv)
        else
           call GetTauHlp2 (Tau(1),T1(1),dima,adda,no,nv)
        end if


        return
        end
