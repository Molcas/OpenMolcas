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
        subroutine ccsd_exc (key)
c
c       check, if there is atleast one determinant in CCSD expansion
c       key=0 - no determinant in expansion
c           1 - only monoexcitations in expansion
c           2 - both mono and biexcitations in expansion
c
        implicit none
#include "ccsd1.fh"
        integer key
c
c       help variables
        integer isym,jsym,ijsym,asym,bsym,nij,nab
        integer naa,nbb,naaaa,nbbbb,nabab
c
c
c1.1    calc # of monoexcitations
c       taking into account also symmetry

        naa=0
        nbb=0
        do isym=1,nsym
          asym=isym
          naa=naa+noa(isym)*nva(asym)
          nbb=nbb+nob(isym)*nvb(asym)
        end do
c
c1.2    calc # of biexcitation
c       taking into account also symmetry
c
        naaaa=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        if (isym.eq.jsym) then
          nij=noa(isym)*(noa(isym)-1)/2
        else
          nij=noa(isym)*noa(jsym)
        end if
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          if (bsym.lt.asym) then
            nab=nva(asym)*nva(bsym)
          else if (bsym.eq.asym) then
            nab=nva(asym)*(nva(asym)-1)/2
          else
            nab=0
          end if
          naaaa=naaaa+nij*nab
          end do
        end do
        end do
c
        nbbbb=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        if (isym.eq.jsym) then
          nij=nob(isym)*(nob(isym)-1)/2
        else
          nij=nob(isym)*nob(jsym)
        end if
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          if (bsym.lt.asym) then
            nab=nvb(asym)*nvb(bsym)
          else if (bsym.eq.asym) then
            nab=nvb(asym)*(nvb(asym)-1)/2
          else
            nab=0
          end if
          nbbbb=nbbbb+nij*nab
          end do
        end do
        end do
c
        nabab=0
        do isym=1,nsym
        do jsym=1,isym
        ijsym=mmul(isym,jsym)
        nij=noa(isym)*nob(jsym)
          do asym=1,nsym
          bsym=mmul(ijsym,asym)
          nab=nva(asym)*nvb(bsym)
          nabab=nabab+nij*nab
          end do
        end do
        end do
c
c
c2      set key
c
        if ((naaaa+nbbbb+nabab).eq.0) then
          if ((naa+nbb).eq.0) then
            key=0
          else
            key=1
          end if
        else
          key=2
        end if
c
c
        return
        end
