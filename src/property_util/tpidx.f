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
      SUBROUTINE tpidx2orb(NSYM,NB,
     $        TYPEINDEX,
     $        NF,NI,N1,N2,N3,NS,ND)
      IMPLICIT NONE
      integer :: NSYM, NB(NSYM)
      integer :: TYPEINDEX(*)
      integer :: NF(*), NI(*), N1(*), N2(*), N3(*), NS(*), ND(*)
      integer :: iSym, iStart
      iStart=1
      do isym=1,nsym
        call tpidx2orb_sym(
     $          TYPEINDEX(iStart),NB(iSym),
     $          NF(iSym),NI(iSym),
     $          N1(iSym),N2(iSym),N3(iSym),
     $          NS(iSym),ND(iSym))
        iStart=iStart+nB(iSym)
      end do
      END

      SUBROUTINE orb2tpidx(NSYM,NB,
     $        NF,NI,N1,N2,N3,NS,ND,
     $        TYPEINDEX)
      IMPLICIT NONE
      integer :: NSYM, NB(NSYM)
      integer :: NF(*), NI(*), N1(*), N2(*), N3(*), NS(*), ND(*)
      integer :: TYPEINDEX(*)
      integer :: iSym, iStart
      iStart=1
      do isym=1,nsym
        call orb2tpidx_sym(
     $          NF(iSym),NI(iSym),
     $          N1(iSym),N2(iSym),N3(iSym),
     $          NS(iSym),ND(iSym),
     $          TYPEINDEX(iStart))
        iStart=iStart+nB(iSym)
      end do
      END

      SUBROUTINE tpstr2tpidx(TYPESTRING,TYPEINDEX,NB)
      IMPLICIT NONE
      character :: TYPESTRING(*)
      integer :: TYPEINDEX(*),NB,i
      do i=1,NB
        select case (TYPESTRING(i))
          case ('F','f')
            TYPEINDEX(i)=1
          case ('I','i')
            TYPEINDEX(i)=2
          case ('1')
            TYPEINDEX(i)=3
          case ('2')
            TYPEINDEX(i)=4
          case ('3')
            TYPEINDEX(i)=5
          case ('S','s')
            TYPEINDEX(i)=6
          case ('D','d')
            TYPEINDEX(i)=7
        end select
      end do
      END

      SUBROUTINE tpidx2tpstr(TYPEINDEX,TYPESTRING,NB)
      IMPLICIT NONE
      character :: TYPESTRING(*)
      integer :: TYPEINDEX(*),NB,i
      do i=1,NB
        select case (TYPEINDEX(i))
          case (1)
            TYPESTRING(i)='F'
          case (2)
            TYPESTRING(i)='I'
          case (3)
            TYPESTRING(i)='1'
          case (4)
            TYPESTRING(i)='2'
          case (5)
            TYPESTRING(i)='3'
          case (6)
            TYPESTRING(i)='S'
          case (7)
            TYPESTRING(i)='D'
        end select
      end do
      END

      SUBROUTINE tpstr2orb(NSYM,NB,
     $        TYPESTRING,
     $        NF,NI,N1,N2,N3,NS,ND)
      IMPLICIT NONE
      integer :: NSYM, NB(NSYM)
      character :: typestring(*)
      integer :: NF(*), NI(*), N1(*), N2(*), N3(*), NS(*), ND(*)
      integer :: iSym, iStart
      iStart=1
      do isym=1,nsym
        call tpstr2orb_sym(
     $          TYPESTRING(iStart),NB(iSym),
     $          NF(iSym),NI(iSym),
     $          N1(iSym),N2(iSym),N3(iSym),
     $          NS(iSym),ND(iSym))
        iStart=iStart+nB(iSym)
      end do
      END

      SUBROUTINE orb2tpstr(NSYM,NB,
     $        NF,NI,N1,N2,N3,NS,ND,
     $        TYPESTRING)
      IMPLICIT NONE
      integer :: NSYM, NB(NSYM)
      integer :: NF(*), NI(*), N1(*), N2(*), N3(*), NS(*), ND(*)
      character :: typestring(*)
      integer :: iSym, iStart
      iStart=1
      do isym=1,nsym
        call orb2tpstr_sym(
     $          NF(iSym),NI(iSym),
     $          N1(iSym),N2(iSym),N3(iSym),
     $          NS(iSym),ND(iSym),
     $          TYPESTRING(iStart))
        iStart=iStart+nB(iSym)
      end do
      END

      SUBROUTINE tpidx2orb_sym(TYPEINDEX,NB,NF,NI,N1,N2,N3,NS,ND)
*SVC: read orbital partition info from a typeindex array
*     corresponding to a _specific_symmetry_ (so these variables are
*     scalars!!
*     A typeindex array consists of integers with 7 possible values
*     corresponding to the types 'fi123sd' -> '1234567'.
      implicit none
      integer :: typeindex(*), NB
      integer :: NF, NI, N1, N2, N3, NS, ND

      integer :: IB, ITYPE
*
      NF=0
      NI=0
      N1=0
      N2=0
      N3=0
      NS=0
      ND=0
* switch typeindex:
      DO IB=1,NB
        ITYPE=TYPEINDEX(IB)
        IF     (ITYPE.EQ.1) THEN
          NF = NF + 1
        ELSE IF(ITYPE.EQ.2) THEN
          NI = NI + 1
        ELSE IF(ITYPE.EQ.3) THEN
          N1 = N1 + 1
        ELSE IF(ITYPE.EQ.4) THEN
          N2 = N2 + 1
        ELSE IF(ITYPE.EQ.5) THEN
          N3 = N3 + 1
        ELSE IF(ITYPE.EQ.6) THEN
          NS = NS + 1
        ELSE IF(ITYPE.EQ.7) THEN
          ND = ND + 1
        ELSE
          WRITE (6,*) 'TPIDX2ORB_SYM: unknown type index number'
          call AbEnd
        END IF
      END DO

      End

      SUBROUTINE orb2tpidx_sym(NF,NI,N1,N2,N3,NS,ND,TYPEINDEX)
*SVC: convert orbital partition info to a typeindex array
*     corresponding to a _specific_symmetry_ (so these variables are
*     scalars!!
*     A typeindex array consists of integers with 7 possible values
*     corresponding to the types 'fi123sd' -> '1234567'.
      implicit none
      integer :: typeindex(*)
      integer :: NF, NI, N1, N2, N3, NS, ND

      integer :: IB, iOff

      iOff=0
*
      DO IB=1,NF
        TYPEINDEX(iOff+IB)=1
      END DO
      iOff=iOff+NF

      DO IB=1,NI
        TYPEINDEX(iOff+IB)=2
      END DO
      iOff=iOff+NI

      DO IB=1,N1
        TYPEINDEX(iOff+IB)=3
      END DO
      iOff=iOff+N1

      DO IB=1,N2
        TYPEINDEX(iOff+IB)=4
      END DO
      iOff=iOff+N2

      DO IB=1,N3
        TYPEINDEX(iOff+IB)=5
      END DO
      iOff=iOff+N3

      DO IB=1,NS
        TYPEINDEX(iOff+IB)=6
      END DO
      iOff=iOff+NS

      DO IB=1,ND
        TYPEINDEX(iOff+IB)=7
      END DO

      End

      SUBROUTINE tpstr2orb_sym(TYPESTRING,NB,NF,NI,N1,N2,N3,NS,ND)
*SVC: read orbital partition info from a typestring array
*     corresponding to a _specific_symmetry_ (so these variables are
*     scalars!!
*     A typestring array consists of characters of 'fi123sd'
      implicit none
      character :: typestring(*)
      integer :: NB
      integer :: NF, NI, N1, N2, N3, NS, ND

      integer :: IB
      character :: CTYPE
*
      NF=0
      NI=0
      N1=0
      N2=0
      N3=0
      NS=0
      ND=0
* switch typeindex:
      DO IB=1,NB
        CTYPE=TYPESTRING(IB)
        CALL UPCASE(CTYPE)
        IF     (CTYPE.EQ.'F') THEN
          NF = NF + 1
        ELSE IF(CTYPE.EQ.'I') THEN
          NI = NI + 1
        ELSE IF(CTYPE.EQ.'1') THEN
          N1 = N1 + 1
        ELSE IF(CTYPE.EQ.'2') THEN
          N2 = N2 + 1
        ELSE IF(CTYPE.EQ.'3') THEN
          N3 = N3 + 1
        ELSE IF(CTYPE.EQ.'S') THEN
          NS = NS + 1
        ELSE IF(CTYPE.EQ.'D') THEN
          ND = ND + 1
        ELSE
          WRITE (6,*) 'TPSTR2ORB_SYM: unknown type index character '
     $            //CTYPE
          call AbEnd
        END IF
      END DO

      End

      SUBROUTINE orb2tpstr_sym(NF,NI,N1,N2,N3,NS,ND,TYPESTRING)
*SVC: convert orbital partition info to a typestring array
*     corresponding to a _specific_symmetry_ (so these variables are
*     scalars!!
*     A typestring array consists of characters of 'fi123sd'
      implicit none
      character :: typestring(*)
      integer :: NF, NI, N1, N2, N3, NS, ND

      integer :: IB, iOff

      iOff=0
*
      DO IB=1,NF
        TYPESTRING(iOff+IB)='F'
      END DO
      iOff=iOff+NF

      DO IB=1,NI
        TYPESTRING(iOff+IB)='I'
      END DO
      iOff=iOff+NI

      DO IB=1,N1
        TYPESTRING(iOff+IB)='1'
      END DO
      iOff=iOff+N1

      DO IB=1,N2
        TYPESTRING(iOff+IB)='2'
      END DO
      iOff=iOff+N2

      DO IB=1,N3
        TYPESTRING(iOff+IB)='3'
      END DO
      iOff=iOff+N3

      DO IB=1,NS
        TYPESTRING(iOff+IB)='S'
      END DO
      iOff=iOff+NS

      DO IB=1,ND
        TYPESTRING(iOff+IB)='D'
      END DO

      End
