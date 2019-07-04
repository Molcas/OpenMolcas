!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_utility_routines

 public dgetsp_util
 public pretty_print_util
 public is_nan
 public lower_to_upper
 public prepare_local_input
 public find_qcmaquis_keyword
 public get_state_tag

contains
      subroutine dgetsp_util(n,age,asp)
!
!     routine taken from Dalton -> dgetsp and adapted
!     originally written on the 8-Feb-1987 by Hans Joergen Aa. Jensen
!
!     Purpose: Transform from GE format to SP format, that is:
!              extract symmetric part of general matrix AGE
!              to symmetric, packed matrix ASP.
!
      implicit none
      real*8 , intent(in)    :: AGE(N,*)
      integer, intent(in)    :: n
      real*8 , intent(inout) :: ASP(*)
      real*8 , parameter     :: DP5 = 0.5D0

      integer                :: i, j, joff
      do j = 1,n
         joff = (j*j-j)/2
         do i = 1,j
            asp(joff+i) = dp5 * (age(i,j) + age(j,i))
         end do
      end do

      end subroutine dgetsp_util

      SUBROUTINE pretty_print_util(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,NCTL,LUPRI)
!.......................................................................
! Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
!         16-Jun-1986 hjaaj ( removed Hollerith )
!
! OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;
!
!        AMATRX(',').........MATRIX TO BE OUTPUT
!
!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
!
!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
!
!        ROWDIM..............ROW DIMENSION OF AMATRX(',')
!
!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
!
!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!                            hjaaj: negative for 132 col width
!
! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.
!
! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971
!
!.......................................................................
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      LOGICAL, external :: IS_NAN
      PARAMETER (ZERO=0.D00, KCOLP=5, KCOLN=8)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
!
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
!
      AMAX = ZERO
      N_NAN = 0
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
!           IF ( IS_NAN(AMATRX(I,J),AMATRX(I,J)) ) THEN
!              N_NAN = N_NAN + 1
!           ELSE
               AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
!           END IF
   10 CONTINUE
      IF (N_NAN .GT. 0) WRITE (LUPRI,'(/T6,A,I10,A)') 'WARNING: matrix contains',N_NAN,' NaN.'
      IF (AMAX <= 1.0d-20) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
!        use F output format
         PFMT = '(A1,I7,2X,8F18.11)'
         thrpri = 0.5D-10
      ELSE
!        use 1PE output format
         PFMT = '(A1,I7,2X,1P,8E15.6)'
         thrpri = 1.0D-10*AMAX
      END IF
!
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
!
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (abs(AMATRX(K,I)).gt.thrpri) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 WRITE(LUPRI,'(A)') '    ==== End of matrix output ===='
      RETURN
 1000 FORMAT (/10X,8(5X,A6,I4))

      end subroutine pretty_print_util

      LOGICAL FUNCTION IS_NAN(XA,XB)
!
!     May 2010, Hans Joergen Aa. Jensen
!     Purpose: IS_NAN(X,X) is true iff X is NAN
!
      REAL*8 XA, XB
      IS_NAN = XA .NE. XB
      END FUNCTION IS_NAN

      subroutine lower_to_upper(str)

        character(*), intent(inout) :: str

        integer                     :: i
        do i = 1, len(str)
          select case(str(i:i))
            case("a":"z")
              str(i:i) = achar(iachar(str(i:i))-32)
          end select
        end do
      end subroutine lower_to_upper

      subroutine prepare_local_input(funit,input_string,strdim,initial_guess,do_restart,iter,&
                                     E_threshold,norb,fiedler_or_cideas)

        integer,            intent(in)               :: strdim
        character(len=500), intent(in), dimension(*) :: input_string
        character(len=2),   intent(inout)            :: initial_guess
        integer,            intent(in)               :: funit
        integer,            intent(in)               :: iter
        logical,            intent(inout)            :: do_restart
        real*8,             intent(in)               :: E_threshold
        integer,            intent(in)               :: norb
        logical,            intent(in)               :: fiedler_or_cideas

        integer                                      :: i, j, nsweeps_tmp, k, local_m, fiedler_sweeps
        integer, save                                :: nsweeps_new
        character(len=500)                           :: string
        character(len=500)                           :: string2
        character(len=500)                           :: mod_string
        character(len=500)                           :: keyword

        !> in general: if a keyword occurs several times in the input, the last occurence counts in QCMaquis

        do i = 1, strdim,2
          string(1:500) = " "
          string        = trim(input_string(i))
          call lower_to_upper(string)
          !> test for exceptions - for example
!            - MOLCAS cuts out the " " for the actual input that is processed
!            - ???
          if(trim(string) == "SWEEP_BOND_DIMENSIONS")then
            write(funit,'(a,a,a,a)') trim(input_string(i)),' = "',trim(input_string(i+1)),'"'
            cycle
          else if(trim(string) == "ORBITAL_ORDER")then
            write(funit,'(a,a,a,a)') trim(input_string(i)),' = "',trim(input_string(i+1)),'"'
            cycle
          else if(trim(string) == "SYMMETRY")then
            write(funit,'(a,a,a,a)') trim(input_string(i)),' = "',trim(input_string(i+1)),'"'
            cycle
          else if(trim(string) == "RESTART")then
            !> deprecated keyword
            cycle
          else if(trim(string) == "INIT_STATE")then
            string2(1:500)     = " "
            string2            = trim(input_string(i+1))
            initial_guess(1:2) = string2(1:2)
            call lower_to_upper(string2)
            if(string2(1:2) == "CI") cycle
          end if
          write(funit,'(a,a,a)') trim(input_string(i)),' = ',trim(input_string(i+1))

        end do

        !> restart
        j = 0; keyword = 'NSWEEPS'
        call find_qcmaquis_keyword(input_string,strdim,keyword,j)
        if(j > 0)then
          read(input_string(j+1),'(i20)') nsweeps_tmp
        else
          print *, 'nsweeps is a mandatory keyword in QCMaquis but not given, I quit...'
          stop 'mandatory QCMaquis input keyword missing'
        end if
        do_restart                = iter > 1
        if(iter == 1) nsweeps_new = 0
        nsweeps_new = nsweeps_new + nsweeps_tmp
        mod_string(1:500) = ' '
        write(mod_string,'(i20)') nsweeps_new
        if(do_restart)then
          write(funit,'(a,a)') 'nsweeps = ',trim(mod_string)
        end if

        !> For the threshold, if not set explicitly by the user default to the parent (QC host program) threshold
        j = 0; keyword = 'CONV_THRESH'
        call find_qcmaquis_keyword(input_string,strdim,keyword,j)
        if(j <= 0)then
          write(funit,'(a,a,E12.6)') 'conv_thresh ',' = ',E_threshold
        end if
        j = 0; keyword = 'TRUNCATION_FINAL'
        call find_qcmaquis_keyword(input_string,strdim,keyword,j)
        if(j <= 0)then
          write(funit,'(a,a,E12.6)') 'truncation_final ',' = ',E_threshold*0.001
        end if
        j = 0; keyword = 'IETL_JCD_TOL'
        call find_qcmaquis_keyword(input_string,strdim,keyword,j)
        if(j <= 0)then
          write(funit,'(a,a,E12.6)') 'ietl_jcd_tol ',' = ',E_threshold*0.001
        end if

        !> orbital order
        j = 0; keyword = 'ORBITAL_ORDER'
        call find_qcmaquis_keyword(input_string,strdim,keyword,j)
        if(j <= 0)then
          write(funit,'(a)',advance='no') 'orbital_order="'
          do i = 1, norb-1
            if(i < 10)then
              write(funit,"(I1,A1)",advance='no') i,','
            else if(i < 100)then
              write(funit,"(I2,A1)",advance='no') i,','
            else if(i < 1000)then
              write(funit,"(I3,A1)",advance='no') i,','
            else
              stop 'more than 999 orbitals in orbital order'
            end if
          end do
          if(i < 10)then
            write(funit,"(I1,A1)",advance='yes') norb,'"'
          else if(i < 100)then
            write(funit,"(I2,A1)",advance='yes') norb,'"'
          else if(i < 1000)then
            write(funit,"(I3,A1)",advance='yes') norb,'"'
          else
            stop 'more than 999 orbitals in orbital order'
          end if
        end if

        !> init state
        if(initial_guess == "hf")then
          j = 0; keyword = "INIT_STATE"
          call find_qcmaquis_keyword(input_string,strdim,keyword,j)
          if(j <= 0)then
            write(funit,'(a,a,a,a)') 'init_state=','"',trim(initial_guess),'"'
          end if
        end if

        !> prepare input for initial fiedler vector ordering determination
        if(fiedler_or_cideas)then
           !> m
                         local_m = 128
           if(norb > 24) local_m = 256
           mod_string(1:500) = ' '
           write(mod_string,'(i20)') local_m
           write(funit,'(a,a)') 'max_bond_dimension = ',trim(mod_string)

           !> nsweeps
           fiedler_sweeps = 4
           mod_string(1:500) = ' '
           write(mod_string,'(i20)') fiedler_sweeps
           write(funit,'(a,a)') 'nsweeps = ',trim(mod_string)

           !> init_state
           write(funit,'(a,a,a,a)') 'init_state=','"','default','"'
        end if

      end subroutine prepare_local_input

      subroutine find_qcmaquis_keyword(input_string,strdim,keyword,location)

        integer,            intent(in)                       :: strdim
        character(len=500), intent(in), dimension(*)         :: input_string
        character(len=500), intent(in)                       :: keyword
        integer,            intent(inout)                    :: location

        integer                                              :: i
        character(len=500)                                   :: string

        location = -1

        do i = 1, strdim,2
          string(1:500) = " "
          string        = trim(input_string(i))
          call lower_to_upper(string)
          if(trim(string) == trim(keyword))then
            location = i
            exit
          end if
        end do

      end subroutine find_qcmaquis_keyword

      subroutine get_state_tag(state,state_tag,offset)

        integer,            intent(in)                       :: state
        character(len=5  ), intent(inout)                    :: state_tag
        integer,            intent(in)                       :: offset

        integer                                              :: irootm1

        state_tag(1:5) = ' '
        irootm1        = state-1+offset
        if(irootm1 < 10)then
          write(state_tag,'(i1)') irootm1
        else if(irootm1 < 100)then
          write(state_tag,'(i2)') irootm1
        else if(irootm1 < 1000)then
          write(state_tag,'(i3)') irootm1
        else if(irootm1 < 10000)then
          write(state_tag,'(i4)') irootm1
        else if(irootm1 < 100000)then
          write(state_tag,'(i5)') irootm1
        end if

      end subroutine get_state_tag

end module qcmaquis_interface_utility_routines

