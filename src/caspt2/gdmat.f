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
      SUBROUTINE GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,OCC,GDAO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CNAT(*),OCC(*)
      DIMENSION GDAO(*)
      DIMENSION ISTART(NSYM),NUSE(NSYM),NBAS(NSYM)
#include "WrkSpc.fh"

* General density matrix, in the sense of using a specified
* but arbitrary range of orbitals in each symmetry, and an
* occupation number.

* The input arguments:
*    NSYM=Nr of symmetry blocks,
*    NBAS(i), i=1..NSYM: Number of basis functions, equal
*       also to total number of orbitals, in each symmetry.
*    CNAT(i),i=1..sum(NBAS(i)**2,i=1..NSYM): The CMO
*       coefficients, stored as square symmetry blocks CN(a,q)
*       where a is basis function, q is orbital index.
*    ISTART(i), i=1..NSYM: Orbital number to start using
*        in each symmetry block, numbered 1,2,..NBAS(i)
*    NUSE(i), i=1..NSYM: How many orbitals to use,
*    OCC(i),i=1..sum(NBAS(i),i=1..NSYM): The natural
*        occupation number of each orbital.
* Output: Array GDAO.

* Computes the AO density matrix, using formula
*  D(a,b) = sum( CN(a,p)*xn(p)*CN(b,p), p=1..n)
* where a,b are basis function indices, p is an MO index,
*  CN() and xn() are CMO coefficients and occupation
* numbers of natural orbitals, and the indices and matrices
* are defined within symmetry blocks.
* The symmetry blocks of D are then stored triangularly
* after each other in the array GDAO.
      IOEND=0
      ICEND=0
      IDAB=0
      DO ISYM=1,NSYM
       NB=NBAS(ISYM)
       IF (NB.GT.0) THEN
         CALL DCOPY_( (NB*(NB+1))/2,0.0D0,0,GDAO(IDAB+1),1)
         NW=NUSE(ISYM)
         IF(NW.GT.0) THEN
          IW1=ISTART(ISYM)
          IW2=IW1-1+NW
          DO IA=1,NB
           DO IB=1,IA
            IDAB=IDAB+1
            DAB=GDAO(IDAB)
            DO IW=IW1,IW2
              DAB=DAB+OCC(IOEND+IW)*CNAT(ICEND+IA+NB*(IW-1))*
     &                             CNAT(ICEND+IB+NB*(IW-1))
            END DO
            GDAO(IDAB)=DAB
           END DO
          END DO
         ELSE
          IDAB=IDAB+(NB*(NB+1))/2
         END IF
         IOEND=IOEND+NB
         ICEND=ICEND+NB**2
       END IF
      END DO

      RETURN
      END
