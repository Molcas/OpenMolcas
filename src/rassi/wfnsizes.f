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
      subroutine wfnsizes_rassi()
      use rasdef, only: NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T
      use rassi_aux, only : nasht_save
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NISHT,NASHT,NOSHT,NSSHT,NBST,
     &                      NAES,NASH,NBASF,NISH,NOSH,NSSH
      implicit none
      integer :: isym
C The structure of the orbital space:
      NISHT=0
      NASHT=0
      NRS1T=0
      NRS2T=0
      NRS3T=0
      NOSHT=0
      NSSHT=0
      NBST=0
      DO ISYM=1,NSYM
        NAES(ISYM)=NASHT
        NISHT=NISHT+NISH(ISYM)
        NASHT=NASHT+NASH(ISYM)
        NRS1T=NRS1T+NRS1(ISYM)
        NRS2T=NRS2T+NRS2(ISYM)
        NRS3T=NRS3T+NRS3(ISYM)
        NOSHT=NOSHT+NOSH(ISYM)
        NSSHT=NSSHT+NSSH(ISYM)
        NBST=NBST+NBASF(ISYM)
      END DO
      nasht_save=nasht
      end subroutine wfnsizes_rassi
