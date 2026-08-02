!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Dmitri Laikov                                          *
!               2010, John Burkardt                                    *
!               2024, Ignacio Fdez. Galvan                             *
!***********************************************************************

! Adapted for OpenMolcas in Feb. 2024

module Lebedev_quadrature

use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp), parameter :: rule_max = 65, &
                                ! availability table
                                available(rule_max) = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0, &
                                                       1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1], &
                                ! precision (degree of accuracy) is just 2*n+1, so odd numbers from 3 to 131
                                ! order (number of points) table
                                order(rule_max,0:4) = reshape([ &
                                                        ! sym = 0  N = 14(+12)+24*(N1+N2)+48*N3
                                                        6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,386,434,482,530,590, &
                                                        650,698,770,830,890,974,1046,1118,1202,1274,1358,1454,1538,1622,1730,1814, &
                                                        1910,2030,2126,2222,2354,2450,2558,2702,2810,2930,3074,3182,3314,3470, &
                                                        3590,3722,3890,4010,4154,4334,4466,4610,4802,4934,5090,5294,5438,5606, &
                                                        5810, &
                                                        ! sym = 1  N = 9(+8)+12*N1+16*N2+24*N3
                                                        5, 9,17,25,29,45,49, 61, 77, 93,105,125,141,161,185,209,229,261,281,309, &
                                                        341,373,401,437,469,505, 545, 585, 621, 669, 705, 749, 797, 845, 889, 941, &
                                                         989,1041,1097,1153,1205,1269,1321,1381,1445,1509,1569,1637,1701,1769, &
                                                        1841,1913,1981,2061,2129,2205,2285,2365,2441,2525,2605,2689,2777,2865, &
                                                        2949, &
                                                        ! sym = 2  N = 6(+5)+6*N1+10*N2+12*N3
                                                        4, 6,11,16,17,27,28, 34, 41, 51, 57, 68, 75, 86, 98,113,121,141,149,162, &
                                                        179,199,209,230,247,262, 284, 306, 321, 351, 366, 386, 413, 440, 457, 488, &
                                                         512, 534, 566, 598, 617, 657, 682, 706, 743, 777, 801, 842, 873, 902, &
                                                         944, 983,1009,1059,1091,1122,1169,1213,1241,1292,1333,1366,1418,1464, &
                                                        1497, &
                                                        ! sym = 3  N = 4(+3)+3*N1+6*(N2+N3)
                                                        3, 4, 7,10,10,16,16, 19, 22, 28, 31, 37, 40, 46, 52, 61, 64, 76, 79, 85, &
                                                         94,106,109,121,130,136, 148, 160, 166, 184, 190, 199, 214, 229, 235, 253, &
                                                         265, 274, 292, 310, 316, 340, 352, 361, 382, 400, 409, 433, 448, 460, &
                                                         484, 505, 514, 544, 559, 571, 598, 622, 631, 661, 682, 694, 724, 748, &
                                                         760, &
                                                        ! sym = 4  N = 7(+6)+12*(N1+N2)+24*N3
                                                        3, 7,13,19,25,37,43, 55, 73, 85, 97,115,133,151,175,193,217,241,265,295, &
                                                        325,349,385,415,445,487, 523, 559, 601, 637, 679, 727, 769, 811, 865, 907, &
                                                         955,1015,1063,1111,1177,1225,1279,1351,1405,1465,1537,1591,1657,1735, &
                                                        1795,1861,1945,2005,2077,2167,2233,2305,2401,2467,2545,2647,2719,2803, &
                                                        2905 &
                                                      ],[rule_max,5])

! AVAILABLE_TABLE returns the availability of a Lebedev rule.
! ORDER_TABLE     returns the order of a Lebedev rule.
! PRECISION_TABLE returns the precision of a Lebedev rule.
! LD_BY_RULE      returns a Lebedev angular grid given its rule number.

public :: available_table, ld_by_rule, order_table, rule_max

contains

function available_table(rule)

  !*********************************************************************
  !
  !! available_table() returns the availability of a Lebedev rule.
  !
  !  Modified:
  !
  !    12 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Input, integer RULE, the index of the rule, between 1 and 65.
  !
  !    Output, integer AVAILABLE_TABLE, the availability of the rule.
  !    * -1, there is no such rule;
  !    *  0, there is such a rule, but it is not available in this library.
  !    *  1, the rule is available in this library.

  integer(kind=iwp) available_table
  integer(kind=iwp), intent(in) :: rule

  if (rule < 1) then
    available_table = -1
  else if (rule_max < rule) then
    available_table = -1
  else
    available_table = available(rule)
  end if

  return

end function available_table

function order_table(sym,rule)

  !*********************************************************************
  !
  !! ORDER_TABLE returns the order of a Lebedev rule.
  !
  !    SYM specifies the type of symmetry:
  !     0: no symmetry
  !     1: z symmetry
  !     2: x,z symmetry
  !     3: x,y,z symmetry
  !     4: inversion symmetry
  !
  !  Modified:
  !
  !    11 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Input, integer RULE, the index of the rule, between 1 and 65.
  !
  !    Output, integer ORDER_TABLE, the order of the rule.

  integer(kind=iwp) :: order_table
  integer(kind=iwp), intent(in) :: sym, rule

  if ((sym < 0) .or. (sym > 4)) then
    write(u6,'(a)') ' '
    write(u6,'(a)') 'ORDER_TABLE - Fatal error!'
    write(u6,'(a)') '  Illegal value of SYM.'
    call Abend()
  end if

  if (rule < 1) then
    write(u6,'(a)') ' '
    write(u6,'(a)') 'ORDER_TABLE - Fatal error!'
    write(u6,'(a)') '  RULE < 1.'
    call Abend()
  else if (rule_max < rule) then
    write(u6,'(a)') ' '
    write(u6,'(a)') 'ORDER_TABLE - Fatal error!'
    write(u6,'(a)') '  RULE_MAX < RULE.'
    call Abend()
  end if

  order_table = order(rule,sym)

  return

end function order_table

subroutine ld_by_rule(sym,rule,x,y,z,w)

  !*********************************************************************
  !
  !! LD_BY_RULE returns a Lebedev angular grid given its rule number.
  !
  !  Discussion:
  !
  !    Only a certain set of such rules are available through this function.
  !
  !    SYM specifies the type of symmetry:
  !     0: no symmetry
  !     1: z symmetry
  !     2: x,z symmetry
  !     3: x,y,z symmetry
  !     4: inversion symmetry
  !
  !    Changed from LD_BY_ORDER, because with symmetry some orders are duplicated
  !
  !  Modified:
  !
  !    13 September 2010
  !    27 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Input, integer RULE, the rule number.
  !
  !    Output, real ( kind = wp ) X(order(rule,sym)), Y(order(rule,sym)), Z(order(rule,sym)), W(order(rule,sym)),
  !    the coordinates and weights of the points.

  integer(kind=iwp), intent(in) :: sym, rule
  real(kind=wp), intent(out) :: x(order(rule,sym)), y(order(rule,sym)), z(order(rule,sym)), w(order(rule,sym))

  if ((sym < 0) .or. (sym > 4)) then
    write(u6,'(a)') ' '
    write(u6,'(a)') 'LD_BY_RULE - Fatal error!'
    write(u6,'(a)') '  Illegal value of SYM.'
    call Abend()
  end if

  select case (rule)
    case (1)
      call ld0006(sym,x,y,z,w)
    case (2)
      call ld0014(sym,x,y,z,w)
    case (3)
      call ld0026(sym,x,y,z,w)
    case (4)
      call ld0038(sym,x,y,z,w)
    case (5)
      call ld0050(sym,x,y,z,w)
    case (6)
      call ld0074(sym,x,y,z,w)
    case (7)
      call ld0086(sym,x,y,z,w)
    case (8)
      call ld0110(sym,x,y,z,w)
    case (9)
      call ld0146(sym,x,y,z,w)
    case (10)
      call ld0170(sym,x,y,z,w)
    case (11)
      call ld0194(sym,x,y,z,w)
    case (12)
      call ld0230(sym,x,y,z,w)
    case (13)
      call ld0266(sym,x,y,z,w)
    case (14)
      call ld0302(sym,x,y,z,w)
    case (15)
      call ld0350(sym,x,y,z,w)
    case (16) ! not available
      call ld0386(sym,x,y,z,w)
    case (17)
      call ld0434(sym,x,y,z,w)
    case (18) ! not available
      call ld0482(sym,x,y,z,w)
    case (19) ! not available
      call ld0530(sym,x,y,z,w)
    case (20)
      call ld0590(sym,x,y,z,w)
    case (21) ! not available
      call ld0650(sym,x,y,z,w)
    case (22) ! not available
      call ld0698(sym,x,y,z,w)
    case (23)
      call ld0770(sym,x,y,z,w)
    case (24) ! not available
      call ld0830(sym,x,y,z,w)
    case (25) ! not available
      call ld0890(sym,x,y,z,w)
    case (26)
      call ld0974(sym,x,y,z,w)
    case (27) ! not available
      call ld1046(sym,x,y,z,w)
    case (28) ! not available
      call ld1118(sym,x,y,z,w)
    case (29)
      call ld1202(sym,x,y,z,w)
    case (30) ! not available
      call ld1274(sym,x,y,z,w)
    case (31) ! not available
      call ld1358(sym,x,y,z,w)
    case (32)
      call ld1454(sym,x,y,z,w)
    case (33) ! not available
      call ld1538(sym,x,y,z,w)
    case (34) ! not available
      call ld1622(sym,x,y,z,w)
    case (35)
      call ld1730(sym,x,y,z,w)
    case (36) ! not available
      call ld1814(sym,x,y,z,w)
    case (37) ! not available
      call ld1910(sym,x,y,z,w)
    case (38)
      call ld2030(sym,x,y,z,w)
    case (39) ! not available
      call ld2126(sym,x,y,z,w)
    case (40) ! not available
      call ld2222(sym,x,y,z,w)
    case (41)
      call ld2354(sym,x,y,z,w)
    case (42) ! not available
      call ld2450(sym,x,y,z,w)
    case (43) ! not available
      call ld2558(sym,x,y,z,w)
    case (44)
      call ld2702(sym,x,y,z,w)
    case (45) ! not available
      call ld2810(sym,x,y,z,w)
    case (46) ! not available
      call ld2930(sym,x,y,z,w)
    case (47)
      call ld3074(sym,x,y,z,w)
    case (48) ! not available
      call ld3182(sym,x,y,z,w)
    case (49) ! not available
      call ld3314(sym,x,y,z,w)
    case (50)
      call ld3470(sym,x,y,z,w)
    case (51) ! not available
      call ld3590(sym,x,y,z,w)
    case (52) ! not available
      call ld3722(sym,x,y,z,w)
    case (53)
      call ld3890(sym,x,y,z,w)
    case (54) ! not available
      call ld4010(sym,x,y,z,w)
    case (55) ! not available
      call ld4154(sym,x,y,z,w)
    case (56)
      call ld4334(sym,x,y,z,w)
    case (57) ! not available
      call ld4466(sym,x,y,z,w)
    case (58) ! not available
      call ld4610(sym,x,y,z,w)
    case (59)
      call ld4802(sym,x,y,z,w)
    case (60) ! not available
      call ld4934(sym,x,y,z,w)
    case (61) ! not available
      call ld5090(sym,x,y,z,w)
    case (62)
      call ld5294(sym,x,y,z,w)
    case (63) ! not available
      call ld5438(sym,x,y,z,w)
    case (64) ! not available
      call ld5606(sym,x,y,z,w)
    case (65)
      call ld5810(sym,x,y,z,w)
    case default
      write(u6,'(a)') ' '
      write(u6,'(a)') 'LD_BY_RULE - Fatal error!'
      write(u6,'(a,i4)') '  Unexpected value of RULE.'
      call Abend()
  end select

  return

end subroutine ld_by_rule

subroutine gen_oh(sym,code,num,a_in,b_in,v,x,y,z,w)

  !*********************************************************************
  !
  !! GEN_OH generates points under OH symmetry.
  !
  !  Discussion:
  !
  !    Given a point on a sphere, specified by A and B, this routine generates
  !    all the equivalent points under OH symmetry, making grid points with
  !    weight V.
  !
  !    The variable NUM is increased by the number of different points
  !    generated.
  !
  !    Depending on CODE, there are from 6 to 48 different but equivalent
  !    points that are generated:
  !
  !      CODE=1:   (0,0,1) etc                                (  6 points)
  !      CODE=2:   (0,A,A) etc, A=1/sqrt(2)                   ( 12 points)
  !      CODE=3:   (A,A,A) etc, A=1/sqrt(3)                   (  8 points)
  !      CODE=4:   (A,A,B) etc, B=sqrt(1-2 A^2)               ( 24 points)
  !      CODE=5:   (A,B,0) etc, B=sqrt(1-A^2), A input        ( 24 points)
  !      CODE=6:   (A,B,C) etc, C=sqrt(1-A^2-B^2), A, B input ( 48 points)
  !
  !    SYM specifies the type of symmetry:
  !     0: no symmetry
  !     1: z symmetry
  !     2: x,z symmetry
  !     3: x,y,z symmetry
  !     4: inversion symmetry
  !
  !  Modified:
  !
  !    11 September 2010
  !    28 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Input, integer CODE, selects the symmetry group.
  !
  !    Input/output, integer NUM, presumably a counter for the
  !    total number of points.  It is incremented by the number of points
  !    generated on this call.
  !
  !    Input, real ( kind = wp ) A, B, information that may be needed to
  !    generate the coordinates of the points (for code = 5 or 6 only).
  !
  !    Input, real ( kind = wp ) V, the weight to be assigned the points.
  !
  !    Output, real ( kind = wp ) X(NUM), Y(NUM), Z(NUM), W(NUM), the coordinates
  !    and weights of the symmetric points generated on this call.

  use Constants, only: Zero, One, Two, Three, Four, Eight, Half

  integer(kind=iwp), parameter :: points(6,0:4) =  reshape([6,12,8,24,24,48, &
                                                            5, 8,4,12,16,24, &
                                                            4, 5,2, 6,10,12, &
                                                            3, 3,1, 3, 6, 6, &
                                                            3, 6,4,12,12,24],[6,5])
  integer(kind=iwp), intent(in) :: sym, code
  integer(kind=iwp), intent(inout) :: num
  real(kind=wp), intent(in) :: a_in, b_in, v
  real(kind=wp), intent(out) :: x(points(code,sym)), y(points(code,sym)), z(points(code,sym)), w(points(code,sym))
  real(kind=wp) :: a, b, c

  if ((sym < 0) .or. (sym > 4)) then
    write(u6,'(a)') ' '
    write(u6,'(a)') 'GEN_OH - Fatal error!'
    write(u6,'(a)') '  Illegal value of SYM.'
    call Abend()
  end if

  select case (code)

    case (1)

      a = One
      select case (sym)
        case (0)
          x(:) = [   a,  -a,Zero,Zero,Zero,Zero]
          y(:) = [Zero,Zero,   a,  -a,Zero,Zero]
          z(:) = [Zero,Zero,Zero,Zero,   a,  -a]
          w(:) = [ One, One, One, One, One, One]*v
        case (1)
          x(:) = [   a,  -a,Zero,Zero,Zero]
          y(:) = [Zero,Zero,   a,  -a,Zero]
          z(:) = [Zero,Zero,Zero,Zero,   a]
          w(:) = [ One, One, One, One, Two]*v
        case (2)
          x(:) = [   a,Zero,Zero,Zero]
          y(:) = [Zero,   a,  -a,Zero]
          z(:) = [Zero,Zero,Zero,   a]
          w(:) = [ Two, One, One, Two]*v
        case (3,4)
          x(:) = [   a,Zero,Zero]
          y(:) = [Zero,   a,Zero]
          z(:) = [Zero,Zero,   a]
          w(:) = [ Two, Two, Two]*v
      end select

    case (2)

      a = sqrt(Half)
      select case (sym)
        case (0)
          x(:) = [Zero,Zero,Zero,Zero,   a,  -a,   a,  -a,   a,  -a,   a,  -a]
          y(:) = [   a,  -a,   a,  -a,Zero,Zero,Zero,Zero,   a,   a,  -a,  -a]
          z(:) = [   a,   a,  -a,  -a,   a,   a,  -a,  -a,Zero,Zero,Zero,Zero]
          w(:) = [ One, One, One, One, One, One, One, One, One, One, One, One]*v
        case (1)
          x(:) = [Zero,Zero,   a,  -a,   a,  -a,   a,  -a]
          y(:) = [   a,  -a,Zero,Zero,   a,   a,  -a,  -a]
          z(:) = [   a,   a,   a,   a,Zero,Zero,Zero,Zero]
          w(:) = [ Two, Two, Two, Two, One, One, One, One]*v
        case (2)
          x(:) = [Zero,Zero,   a,   a,   a]
          y(:) = [   a,  -a,Zero,   a,  -a]
          z(:) = [   a,   a,   a,Zero,Zero]
          w(:) = [ Two, Two,Four, Two, Two]*v
        case (3)
          x(:) = [Zero,   a,   a]
          y(:) = [   a,Zero,   a]
          z(:) = [   a,   a,Zero]
          w(:) = [Four,Four,Four]*v
        case (4)
          x(:) = [Zero,Zero,   a,  -a,   a,  -a]
          y(:) = [   a,  -a,Zero,Zero,   a,   a]
          z(:) = [   a,   a,   a,   a,Zero,Zero]
          w(:) = [ Two, Two, Two, Two, Two, Two]*v
      end select

    case (3)

      a = sqrt(One/Three)
      select case (sym)
        case (0)
          x(:) = [  a, -a,  a, -a,  a, -a,  a, -a]
          y(:) = [  a,  a, -a, -a,  a,  a, -a, -a]
          z(:) = [  a,  a,  a,  a, -a, -a, -a, -a]
          w(:) = [One,One,One,One,One,One,One,One]*v
        case (1,4)
          x(:) = [  a, -a,  a, -a]
          y(:) = [  a,  a, -a, -a]
          z(:) = [  a,  a,  a,  a]
          w(:) = [Two,Two,Two,Two]*v
        case (2)
          x(:) = [   a,   a]
          y(:) = [   a,  -a]
          z(:) = [   a,   a]
          w(:) = [Four,Four]*v
        case (3)
          x(:) = [    a]
          y(:) = [    a]
          z(:) = [    a]
          w(:) = [Eight]*v
      end select

    case (4)

      a = a_in
      b = sqrt(One-Two*a_in*a_in)
      select case (sym)
        case (0)
          x(:) = [  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  b, -b,  b, -b,  b, -b,  b, -b]
          y(:) = [  a,  a, -a, -a,  a,  a, -a, -a,  b,  b, -b, -b,  b,  b, -b, -b,  a,  a, -a, -a,  a,  a, -a, -a]
          z(:) = [  b,  b,  b,  b, -b, -b, -b, -b,  a,  a,  a,  a, -a, -a, -a, -a,  a,  a,  a,  a, -a, -a, -a, -a]
          w(:) = [One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One]*v
        case (1,4)
          x(:) = [  a, -a,  a, -a,  a, -a,  a, -a,  b, -b,  b, -b]
          y(:) = [  a,  a, -a, -a,  b,  b, -b, -b,  a,  a, -a, -a]
          z(:) = [  b,  b,  b,  b,  a,  a,  a,  a,  a,  a,  a,  a]
          w(:) = [Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two]*v
        case (2)
          x(:) = [   a,   a,   a,   a,   b,   b]
          y(:) = [   a,  -a,   b,  -b,   a,  -a]
          z(:) = [   b,   b,   a,   a,   a,   a]
          w(:) = [Four,Four,Four,Four,Four,Four]*v
        case (3)
          x(:) = [    a,    a,    b]
          y(:) = [    a,    b,    a]
          z(:) = [    b,    a,    a]
          w(:) = [Eight,Eight,Eight]*v
      end select

    case (5)

      a = a_in
      b = sqrt(One-a_in*a_in)
      select case (sym)
        case (0)
          x(:) = [   a,  -a,   a,  -a,   b,  -b,   b,  -b,   a,  -a,   a,  -a,   b,  -b,   b,  -b,Zero,Zero,Zero,Zero,Zero,Zero, &
                  Zero,Zero]
          y(:) = [   b,   b,  -b,  -b,   a,   a,  -a,  -a,Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,   a,  -a,   a,  -a,   b,  -b, &
                     b,  -b]
          z(:) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,   b,   b,  -b,  -b,   a,   a,  -a,  -a,   b,   b,  -b,  -b,   a,   a, &
                    -a,  -a]
          w(:) = [ One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, One, &
                   One, One]*v
        case (1)
          x(:) = [   a,  -a,   a,  -a,   b,  -b,   b,  -b,   a,  -a,   b,  -b,Zero,Zero,Zero,Zero]
          y(:) = [   b,   b,  -b,  -b,   a,   a,  -a,  -a,Zero,Zero,Zero,Zero,   a,  -a,   b,  -b]
          z(:) = [Zero,Zero,Zero,Zero,Zero,Zero,Zero,Zero,   b,   b,   a,   a,   b,   b,   a,   a]
          w(:) = [ One, One, One, One, One, One, One, One, Two, Two, Two, Two, Two, Two, Two, Two]*v
        case (2)
          x(:) = [   a,   a,   b,   b,   a,   b,Zero,Zero,Zero,Zero]
          y(:) = [   b,  -b,   a,  -a,Zero,Zero,   a,  -a,   b,  -b]
          z(:) = [Zero,Zero,Zero,Zero,   b,   a,   b,   b,   a,   a]
          w(:) = [ Two, Two, Two, Two,Four,Four, Two, Two, Two, Two]*v
        case (3)
          x(:) = [   a,   b,   a,   b,Zero,Zero]
          y(:) = [   b,   a,Zero,Zero,   a,   b]
          z(:) = [Zero,Zero,   b,   a,   b,   a]
          w(:) = [Four,Four,Four,Four,Four,Four]*v
        case (4)
          x(:) = [   a,  -a,   b,  -b,   a,  -a,   b,  -b,Zero,Zero,Zero,Zero]
          y(:) = [   b,   b,   a,   a,Zero,Zero,Zero,Zero,   a,  -a,   b,  -b]
          z(:) = [Zero,Zero,Zero,Zero,   b,   b,   a,   a,   b,   b,   a,   a]
          w(:) = [ Two, Two, Two, Two, Two, Two, Two, Two, Two, Two, Two, Two]*v
      end select

    case (6)

      a = a_in
      b = b_in
      c = sqrt(One-a_in*a_in-b_in*b_in)
      select case (sym)
        case (0)
          x(:) = [  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  a, -a,  b, -b,  b, -b,  b, -b,  b, -b,  b, -b,  b, -b, &
                    b, -b,  b, -b,  c, -c,  c, -c,  c, -c,  c, -c,  c, -c,  c, -c,  c, -c,  c, -c]
          y(:) = [  b,  b, -b, -b,  b,  b, -b, -b,  c,  c, -c, -c,  c,  c, -c, -c,  a,  a, -a, -a,  a,  a, -a, -a,  c,  c, -c, -c, &
                    c,  c, -c, -c,  a,  a, -a, -a,  a,  a, -a, -a,  b,  b, -b, -b,  b,  b, -b, -b]
          z(:) = [  c,  c,  c,  c, -c, -c, -c, -c,  b,  b,  b,  b, -b, -b, -b, -b,  c,  c,  c,  c, -c, -c, -c, -c,  a,  a,  a,  a, &
                   -a, -a, -a, -a,  b,  b,  b,  b, -b, -b, -b, -b,  a,  a,  a,  a, -a, -a, -a, -a]
          w(:) = [One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One, &
                  One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One,One]*v
        case (1,4)
          x(:) = [  a, -a,  a, -a,  a, -a,  a, -a,  b, -b,  b, -b,  b, -b,  b, -b,  c, -c,  c, -c,  c, -c,  c, -c]
          y(:) = [  b,  b, -b, -b,  c,  c, -c, -c,  a,  a, -a, -a,  c,  c, -c, -c,  a,  a, -a, -a,  b,  b, -b, -b]
          z(:) = [  c,  c,  c,  c,  b,  b,  b,  b,  c,  c,  c,  c,  a,  a,  a,  a,  b,  b,  b,  b,  a,  a,  a,  a]
          w(:) = [Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two,Two]*v
        case (2)
          x(:) = [   a,   a,   a,   a,   b,   b,   b,   b,   c,   c,   c,   c]
          y(:) = [   b,  -b,   c,  -c,   a,  -a,   c,  -c,   a,  -a,   b,  -b]
          z(:) = [   c,   c,   b,   b,   c,   c,   a,   a,   b,   b,   a,   a]
          w(:) = [Four,Four,Four,Four,Four,Four,Four,Four,Four,Four,Four,Four]*v
        case (3)
          x(:) = [    a,    a,    b,    b,    c,    c]
          y(:) = [    b,    c,    a,    c,    a,    b]
          z(:) = [    c,    b,    c,    a,    b,    a]
          w(:) = [Eight,Eight,Eight,Eight,Eight,Eight]*v
      end select

    case default

      write(u6,'(a)') ' '
      write(u6,'(a)') 'GEN_OH - Fatal error!'
      write(u6,'(a)') '  Illegal value of CODE.'
      call Abend()

  end select

  num = num+points(code,sym)

  return

end subroutine gen_oh

subroutine ld_all(sym,N,N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  !*********************************************************************
  !
  !! LD_ALL computes a generic Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    28 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  integer(kind=iwp), intent(in) :: N, N1, N2, N3
  real(kind=wp), intent(in) :: a4(N1), a5(N2), a6(N3), b6(N3), v1, v2, v3, v4(N1), v5(N2), v6(N3)
  real(kind=wp), intent(out) :: x(N), y(N), z(N), w(N)
  integer(kind=iwp) :: i, m

  m = 1
  if (v1 /= Zero) call gen_oh(sym,1,m,Zero,Zero,v1,x(m),y(m),z(m),w(m))
  if (v2 /= Zero) call gen_oh(sym,2,m,Zero,Zero,v2,x(m),y(m),z(m),w(m))
  if (v3 /= Zero) call gen_oh(sym,3,m,Zero,Zero,v3,x(m),y(m),z(m),w(m))
  do i=1,N1
    call gen_oh(sym,4,m,a4(i),Zero,v4(i),x(m),y(m),z(m),w(m))
  end do
  do i=1,N2
    call gen_oh(sym,5,m,a5(i),Zero,v5(i),x(m),y(m),z(m),w(m))
  end do
  do i=1,N3
    call gen_oh(sym,6,m,a6(i),b6(i),v6(i),x(m),y(m),z(m),w(m))
  end do

  return

end subroutine ld_all

subroutine ld0006(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0006 computes the 6 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero, One, Six

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(1,sym)), y(order(1,sym)), z(order(1,sym)), w(order(1,sym))
  integer(kind=iwp), parameter :: N1 = 0, N2 = 0, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [real(kind=wp)::], &
                              a5(N2) = [real(kind=wp)::], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = One/Six, v2 = Zero, v3 = Zero, &
                              v4(N1) = [real(kind=wp)::], &
                              v5(N2) = [real(kind=wp)::], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(1,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0006

subroutine ld0014(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0014 computes the 14 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero, One, Three

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(2,sym)), y(order(2,sym)), z(order(2,sym)), w(order(2,sym))
  integer(kind=iwp), parameter :: N1 = 0, N2 = 0, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [real(kind=wp)::], &
                              a5(N2) = [real(kind=wp)::], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = One/15.0_wp, v2 = Zero, v3 = Three/40.0_wp, &
                              v4(N1) = [real(kind=wp)::], &
                              v5(N2) = [real(kind=wp)::], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(2,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0014

subroutine ld0026(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0026 computes the 26 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: One, Four, Nine

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(3,sym)), y(order(3,sym)), z(order(3,sym)), w(order(3,sym))
  integer(kind=iwp), parameter :: N1 = 0, N2 = 0, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [real(kind=wp)::], &
                              a5(N2) = [real(kind=wp)::], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = One/21.0_wp, v2 = Four/105.0_wp, v3 = Nine/280.0_wp, &
                              v4(N1) = [real(kind=wp)::], &
                              v5(N2) = [real(kind=wp)::], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(3,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0026

subroutine ld0038(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0038 computes the 38 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(4,sym)), y(order(4,sym)), z(order(4,sym)), w(order(4,sym))
  integer(kind=iwp), parameter :: N1 = 0, N2 = 1, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [real(kind=wp)::], &
                              a5(N2) = [0.4597008433809831_wp], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = 0.9523809523809524e-2_wp, v2 = Zero, v3 = 0.3214285714285714e-1_wp, &
                              v4(N1) = [real(kind=wp)::], &
                              v5(N2) = [0.2857142857142857e-1_wp], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(4,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0038

subroutine ld0050(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0050 computes the 50 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(5,sym)), y(order(5,sym)), z(order(5,sym)), w(order(5,sym))
  integer(kind=iwp), parameter :: N1 = 1, N2 = 0, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [0.3015113445777636_wp], &
                              a5(N2) = [real(kind=wp)::], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = 0.1269841269841270e-1_wp, v2 = 0.2257495590828924e-1_wp, v3 = 0.2109375000000000e-1_wp, &
                              v4(N1) = [0.2017333553791887e-1_wp], &
                              v5(N2) = [real(kind=wp)::], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(5,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0050

subroutine ld0074(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0074 computes the 74 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(6,sym)), y(order(6,sym)), z(order(6,sym)), w(order(6,sym))
  integer(kind=iwp), parameter :: N1 = 1, N2 = 1, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [0.4803844614152614_wp], &
                              a5(N2) = [0.2657620708215946e-1_wp], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = 0.5130671797338464e-3_wp, v2 = 0.1660406956574204e-1_wp, v3 = -0.2958603896103896e-1_wp, &
                              v4(N1) = [0.3207726489807764_wp], &
                              v5(N2) = [0.1652217099371571e-1_wp], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(6,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0074

subroutine ld0086(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0086 computes the 86 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(7,sym)), y(order(7,sym)), z(order(7,sym)), w(order(7,sym))
  integer(kind=iwp), parameter :: N1 = 2, N2 = 1, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [0.3696028464541502_wp,0.6943540066026664_wp], &
                              a5(N2) = [0.3742430390903412_wp], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = 0.1154401154401154e-1_wp, v2 = Zero, v3 = 0.1194390908585628e-1_wp, &
                              v4(N1) = [0.1111055571060340e-1_wp,0.1187650129453714e-1_wp], &
                              v5(N2) = [0.1181230374690448e-1_wp], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(7,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0086

subroutine ld0110(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0110 computes the 110 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(8,sym)), y(order(8,sym)), z(order(8,sym)), w(order(8,sym))
  integer(kind=iwp), parameter :: N1 = 3, N2 = 1, N3 = 0
  real(kind=wp), parameter :: a4(N1) = [0.1851156353447362_wp,0.6904210483822922_wp,0.3956894730559419_wp], &
                              a5(N2) = [0.4783690288121502_wp], &
                              a6(N3) = [real(kind=wp)::], &
                              b6(N3) = [real(kind=wp)::], &
                              v1 = 0.3828270494937162e-2_wp, v2 = Zero, v3 = 0.9793737512487512e-2_wp, &
                              v4(N1) = [0.8211737283191111e-2_wp,0.9942814891178103e-2_wp,0.9595471336070963e-2_wp], &
                              v5(N2) = [0.9694996361663028e-2_wp], &
                              v6(N3) = [real(kind=wp)::]

  call ld_all(sym,order(8,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0110

subroutine ld0146(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0146 computes the 146 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(9,sym)), y(order(9,sym)), z(order(9,sym)), w(order(9,sym))
  integer(kind=iwp), parameter :: N1 = 3, N2 = 0, N3 = 1
  real(kind=wp), parameter :: a4(N1) = [0.6764410400114264_wp,0.4174961227965453_wp,0.1574676672039082_wp], &
                              a5(N2) = [real(kind=wp)::], &
                              a6(N3) = [0.1403553811713183_wp], &
                              b6(N3) = [0.4493328323269557_wp], &
                              v1 = 0.5996313688621381e-3_wp, v2 = 0.7372999718620756e-2_wp, v3 = 0.7210515360144488e-2_wp, &
                              v4(N1) = [0.7116355493117555e-2_wp,0.6753829486314477e-2_wp,0.7574394159054034e-2_wp], &
                              v5(N2) = [real(kind=wp)::], &
                              v6(N3) = [0.6991087353303262e-2_wp]

  call ld_all(sym,order(9,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0146

subroutine ld0170(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0170 computes the 170 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(10,sym)), y(order(10,sym)), z(order(10,sym)), w(order(10,sym))
  integer(kind=iwp), parameter :: N1 = 3, N2 = 1, N3 = 1
  real(kind=wp), parameter :: a4(N1) = [0.2551252621114134_wp,0.6743601460362766_wp,0.4318910696719410_wp], &
                              a5(N2) = [0.2613931360335988_wp], &
                              a6(N3) = [0.4990453161796037_wp], &
                              b6(N3) = [0.1446630744325115_wp], &
                              v1 = 0.5544842902037365e-2_wp, v2 = 0.6071332770670752e-2_wp, v3 = 0.6383674773515093e-2_wp, &
                              v4(N1) = [0.5183387587747790e-2_wp,0.6317929009813725e-2_wp,0.6201670006589077e-2_wp], &
                              v5(N2) = [0.5477143385137348e-2_wp], &
                              v6(N3) = [0.5968383987681156e-2_wp]

  call ld_all(sym,order(10,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0170

subroutine ld0194(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0194 computes the 194 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(11,sym)), y(order(11,sym)), z(order(11,sym)), w(order(11,sym))
  integer(kind=iwp), parameter :: N1 = 4, N2 = 1, N3 = 1
  real(kind=wp), parameter :: a4(N1) = [0.6712973442695226_wp,0.2892465627575439_wp,0.4446933178717437_wp,0.1299335447650067_wp], &
                              a5(N2) = [0.3457702197611283_wp], &
                              a6(N3) = [0.1590417105383530_wp], &
                              b6(N3) = [0.8360360154824589_wp], &
                              v1 = 0.1782340447244611e-2_wp, v2 = 0.5716905949977102e-2_wp, v3 = 0.5573383178848738e-2_wp, &
                              v4(N1) = [0.5608704082587997e-2_wp,0.5158237711805383e-2_wp,0.5518771467273614e-2_wp, &
                                        0.4106777028169394e-2_wp], &
                              v5(N2) = [0.5051846064614808e-2_wp], &
                              v6(N3) = [0.5530248916233094e-2_wp]

  call ld_all(sym,order(11,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0194

subroutine ld0230(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0230 computes the 230 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(12,sym)), y(order(12,sym)), z(order(12,sym)), w(order(12,sym))
  integer(kind=iwp), parameter :: N1 = 5, N2 = 2, N3 = 1
  real(kind=wp), parameter :: a4(N1) = [0.4492044687397611_wp,0.2520419490210201_wp,0.6981906658447242_wp,0.6587405243460960_wp, &
                                        0.4038544050097660e-1_wp], &
                              a5(N2) = [0.5823842309715585_wp,0.3545877390518688_wp], &
                              a6(N3) = [0.2272181808998187_wp], &
                              b6(N3) = [0.4864661535886647_wp], &
                              v1 = -0.5522639919727325e-1_wp, v2 = Zero, v3 = 0.4450274607445226e-2_wp, &
                              v4(N1) = [0.4496841067921404e-2_wp,0.5049153450478750e-2_wp,0.3976408018051883e-2_wp, &
                                        0.4401400650381014e-2_wp,0.1724544350544401e-1_wp], &
                              v5(N2) = [0.4231083095357343e-2_wp,0.5198069864064399e-2_wp], &
                              v6(N3) = [0.4695720972568883e-2_wp]

  call ld_all(sym,order(12,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0230

subroutine ld0266(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0266 computes the 266 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(13,sym)), y(order(13,sym)), z(order(13,sym)), w(order(13,sym))
  integer(kind=iwp), parameter :: N1 = 5, N2 = 1, N3 = 2
  real(kind=wp), parameter :: a4(N1) = [0.7039373391585475_wp,0.1012526248572414_wp,0.4647448726420539_wp,0.3277420654971629_wp, &
                                        0.6620338663699974_wp], &
                              a5(N2) = [0.8506508083520399_wp], &
                              a6(N3) = [0.3233484542692899_wp,0.2314790158712601_wp], &
                              b6(N3) = [0.1153112011009701_wp,0.5244939240922365_wp], &
                              v1 = -0.1313769127326952e-2_wp, v2 = -0.2522728704859336e-2_wp, v3 = 0.4186853881700583e-2_wp, &
                              v4(N1) = [0.5315167977810885e-2_wp,0.4047142377086219e-2_wp,0.4112482394406990e-2_wp, &
                                        0.3595584899758782e-2_wp,0.4256131351428158e-2_wp], &
                              v5(N2) = [0.4229582700647240e-2_wp], &
                              v6(N3) = [0.4080914225780505e-2_wp,0.4071467593830964e-2_wp]

  call ld_all(sym,order(13,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0266

subroutine ld0302(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0302 computes the 302 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(14,sym)), y(order(14,sym)), z(order(14,sym)), w(order(14,sym))
  integer(kind=iwp), parameter :: N1 = 6, N2 = 2, N3 = 2
  real(kind=wp), parameter :: a4(N1) = [0.3515640345570105_wp,0.6566329410219612_wp,0.4729054132581005_wp, &
                                        0.9618308522614784e-1_wp,0.2219645236294178_wp,0.7011766416089545_wp], &
                              a5(N2) = [0.2644152887060663_wp,0.5718955891878961_wp], &
                              a6(N3) = [0.2510034751770465_wp,0.1233548532583327_wp], &
                              b6(N3) = [0.8000727494073952_wp,0.4127724083168531_wp], &
                              v1 = 0.8545911725128148e-3_wp, v2 = Zero, v3 = 0.3599119285025571e-2_wp, &
                              v4(N1) = [0.3449788424305883e-2_wp,0.3604822601419882e-2_wp,0.3576729661743367e-2_wp, &
                                        0.2352101413689164e-2_wp,0.3108953122413675e-2_wp,0.3650045807677255e-2_wp], &
                              v5(N2) = [0.2982344963171804e-2_wp,0.3600820932216460e-2_wp], &
                              v6(N3) = [0.3571540554273387e-2_wp,0.3392312205006170e-2_wp]

  call ld_all(sym,order(14,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0302

subroutine ld0350(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0350 computes the 350 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(15,sym)), y(order(15,sym)), z(order(15,sym)), w(order(15,sym))
  integer(kind=iwp), parameter :: N1 = 6, N2 = 2, N3 = 3
  real(kind=wp), parameter :: a4(N1) = [0.7068965463912316_wp,0.4794682625712025_wp,0.1927533154878019_wp,0.6930357961327123_wp, &
                                        0.3608302115520091_wp,0.6498486161496169_wp], &
                              a5(N2) = [0.1932945013230339_wp,0.3800494919899303_wp], &
                              a6(N3) = [0.2899558825499574_wp,0.9684121455103957e-1_wp,0.1833434647041659_wp], &
                              b6(N3) = [0.7934537856582316_wp,0.8280801506686862_wp,0.9074658265305127_wp], &
                              v1 = 0.3006796749453936e-2_wp, v2 = Zero, v3 = 0.3050627745650771e-2_wp, &
                              v4(N1) = [0.1621104600288991e-2_wp,0.3005701484901752e-2_wp,0.2990992529653774e-2_wp, &
                                        0.2982170644107595e-2_wp,0.2721564237310992e-2_wp,0.3033513795811141e-2_wp], &
                              v5(N2) = [0.3007949555218533e-2_wp,0.2881964603055307e-2_wp], &
                              v6(N3) = [0.2958357626535696e-2_wp,0.3036020026407088e-2_wp,0.2832187403926303e-2_wp]

  call ld_all(sym,order(15,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0350

subroutine ld0386(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(16,sym)), y(order(16,sym)), z(order(16,sym)), w(order(16,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 6, N2 = 3, N3 = 3
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(16,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0386

subroutine ld0434(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0434 computes the 434 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(17,sym)), y(order(17,sym)), z(order(17,sym)), w(order(17,sym))
  integer(kind=iwp), parameter :: N1 = 7, N2 = 2, N3 = 4
  real(kind=wp), parameter :: a4(N1) = [0.6909346307509111_wp,0.1774836054609158_wp,0.4914342637784746_wp,0.6456664707424256_wp, &
                                        0.2861289010307638_wp,0.7568084367178018e-1_wp,0.3927259763368002_wp], &
                              a5(N2) = [0.8818132877794288_wp,0.9776428111182649_wp], &
                              a6(N3) = [0.2054823696403044_wp,0.5905157048925271_wp,0.5550152361076807_wp,0.9371809858553722_wp], &
                              b6(N3) = [0.8689460322872412_wp,0.7999278543857286_wp,0.7717462626915901_wp,0.3344363145343455_wp], &
                              v1 = 0.5265897968224436e-3_wp, v2 = 0.2548219972002607e-2_wp, v3 = 0.2512317418927307e-2_wp, &
                              v4(N1) = [0.2530403801186355e-2_wp,0.2014279020918528e-2_wp,0.2501725168402936e-2_wp, &
                                        0.2513267174597564e-2_wp,0.2302694782227416e-2_wp,0.1462495621594614e-2_wp, &
                                        0.2445373437312980e-2_wp], &
                              v5(N2) = [0.2417442375638981e-2_wp,0.1910951282179532e-2_wp], &
                              v6(N3) = [0.2416930044324775e-2_wp,0.2512236854563495e-2_wp,0.2496644054553086e-2_wp, &
                                        0.2236607760437849e-2_wp]

  call ld_all(sym,order(17,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0434

subroutine ld0482(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(18,sym)), y(order(18,sym)), z(order(18,sym)), w(order(18,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 7, N2 = 4, N3 = 4
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(18,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0482

subroutine ld0530(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(19,sym)), y(order(19,sym)), z(order(19,sym)), w(order(19,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 8, N2 = 3, N3 = 5
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(19,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0530

subroutine ld0590(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0590 computes the 590 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(20,sym)), y(order(20,sym)), z(order(20,sym)), w(order(20,sym))
  integer(kind=iwp), parameter :: N1 = 9, N2 = 3, N3 = 6
  real(kind=wp), parameter :: a4(N1) = [0.7040954938227469_wp,0.6807744066455243_wp,0.6372546939258752_wp,0.5044419707800358_wp, &
                                        0.4215761784010967_wp,0.3317920736472123_wp,0.2384736701421887_wp,0.1459036449157763_wp, &
                                        0.6095034115507196e-1_wp], &
                              a5(N2) = [0.6116843442009876_wp,0.3964755348199858_wp,0.1724782009907724_wp], &
                              a6(N3) = [0.5610263808622060_wp,0.4742392842551980_wp,0.5984126497885380_wp,0.3791035407695563_wp, &
                                        0.2778673190586244_wp,0.5033564271075117_wp], &
                              b6(N3) = [0.3518280927733519_wp,0.2634716655937950_wp,0.1816640840360209_wp,0.1720795225656878_wp, &
                                        0.8213021581932511e-1_wp,0.8999205842074875e-1_wp], &
                              v1 = 0.3095121295306187e-3_wp, v2 = Zero, v3 = 0.1852379698597489e-2_wp, &
                              v4(N1) = [0.1871790639277744e-2_wp,0.1858812585438317e-2_wp,0.1852028828296213e-2_wp, &
                                        0.1846715956151242e-2_wp,0.1818471778162769e-2_wp,0.1749564657281154e-2_wp, &
                                        0.1617210647254411e-2_wp,0.1384737234851692e-2_wp,0.9764331165051050e-3_wp], &
                              v5(N2) = [0.1857161196774078e-2_wp,0.1705153996395864e-2_wp,0.1300321685886048e-2_wp], &
                              v6(N3) = [0.1842866472905286e-2_wp,0.1802658934377451e-2_wp,0.1849830560443660e-2_wp, &
                                        0.1713904507106709e-2_wp,0.1555213603396808e-2_wp,0.1802239128008525e-2_wp]

  call ld_all(sym,order(20,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0590

subroutine ld0650(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(21,sym)), y(order(21,sym)), z(order(21,sym)), w(order(21,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 9, N2 = 3, N3 = 7
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(21,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0650

subroutine ld0698(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(22,sym)), y(order(22,sym)), z(order(22,sym)), w(order(22,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 9, N2 = 5, N3 = 7
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(22,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0698

subroutine ld0770(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0770 computes the 770 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(23,sym)), y(order(23,sym)), z(order(23,sym)), w(order(23,sym))
  integer(kind=iwp), parameter :: N1 = 10, N2 = 3, N3 = 9
  real(kind=wp), parameter :: a4(N1) = [0.5087204410502360e-1_wp,0.1228198790178831_wp,0.2026890814408786_wp, &
                                        0.2847745156464294_wp,0.3656719078978026_wp,0.4428264886713469_wp,0.5140619627249735_wp, &
                                        0.6306401219166803_wp,0.6716883332022612_wp,0.6979792685336881_wp], &
                              a5(N2) = [0.1446865674195309_wp,0.3390263475411216_wp,0.5335804651263506_wp], &
                              a6(N3) = [0.6944024393349413e-1_wp,0.2269004109529460_wp,0.8025574607775339e-1_wp, &
                                        0.1467999527896572_wp,0.1571507769824727_wp,0.2365702993157246_wp, &
                                        0.7714815866765732e-1_wp,0.3062936666210730_wp,0.3822477379524787_wp], &
                              b6(N3) = [0.2355187894242326_wp,0.4102182474045730_wp,0.6214302417481605_wp,0.3245284345717394_wp, &
                                        0.5224482189696630_wp,0.6017546634089558_wp,0.4346575516141163_wp,0.4908826589037616_wp, &
                                        0.5648768149099500_wp], &
                              v1 = 0.2192942088181184e-3_wp, v2 = 0.1436433617319080e-2_wp, v3 = 0.1421940344335877e-2_wp, &
                              v4(N1) = [0.6798123511050502e-3_wp,0.9913184235294912e-3_wp,0.1180207833238949e-2_wp, &
                                        0.1296599602080921e-2_wp,0.1365871427428316e-2_wp,0.1402988604775325e-2_wp, &
                                        0.1418645563595609e-2_wp,0.1421376741851662e-2_wp,0.1423996475490962e-2_wp, &
                                        0.1431554042178567e-2_wp], &
                              v5(N2) = [0.9254401499865368e-3_wp,0.1250239995053509e-2_wp,0.1394365843329230e-2_wp], &
                              v6(N3) = [0.1127089094671749e-2_wp,0.1345753760910670e-2_wp,0.1424957283316783e-2_wp, &
                                        0.1261523341237750e-2_wp,0.1392547106052696e-2_wp,0.1418761677877656e-2_wp, &
                                        0.1338366684479554e-2_wp,0.1393700862676131e-2_wp,0.1415914757466932e-2_wp]

  call ld_all(sym,order(23,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0770

subroutine ld0830(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(24,sym)), y(order(24,sym)), z(order(24,sym)), w(order(24,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 11, N2 = 5, N3 = 9
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(24,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0830

subroutine ld0890(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(25,sym)), y(order(25,sym)), z(order(25,sym)), w(order(25,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 11, N2 = 5, N3 = 10
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(25,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0890

subroutine ld0974(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD0974 computes the 974 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(26,sym)), y(order(26,sym)), z(order(26,sym)), w(order(26,sym))
  integer(kind=iwp), parameter :: N1 = 12, N2 = 4, N3 = 12
  real(kind=wp), parameter :: a4(N1) = [0.4292963545341347e-1_wp,0.1051426854086404_wp,0.1750024867623087_wp, &
                                        0.2477653379650257_wp,0.3206567123955957_wp,0.3916520749849983_wp,0.4590825874187624_wp, &
                                        0.5214563888415861_wp,0.6253170244654199_wp,0.6637926744523170_wp,0.6910410398498301_wp, &
                                        0.7052907007457760_wp], &
                              a5(N2) = [0.1236686762657990_wp,0.2940777114468387_wp,0.4697753849207649_wp,0.6334563241139567_wp], &
                              a6(N3) = [0.5974048614181342e-1_wp,0.1375760408473636_wp,0.3391016526336286_wp, &
                                        0.1271675191439820_wp,0.2693120740413512_wp,0.1419786452601918_wp, &
                                        0.6709284600738255e-1_wp,0.7057738183256172e-1_wp,0.2783888477882155_wp, &
                                        0.1979578938917407_wp,0.2087307061103274_wp,0.4055122137872836_wp], &
                              b6(N3) = [0.2029128752777523_wp,0.4602621942484054_wp,0.5030673999662036_wp,0.2817606422442134_wp, &
                                        0.4331561291720157_wp,0.6256167358580814_wp,0.3798395216859157_wp,0.5517505421423520_wp, &
                                        0.6029619156159187_wp,0.3589606329589096_wp,0.5348666438135476_wp,0.5674997546074373_wp], &
                              v1 = 0.1438294190527431e-3_wp, v2 = Zero, v3 = 0.1125772288287004e-2_wp, &
                              v4(N1) = [0.4948029341949241e-3_wp,0.7357990109125470e-3_wp,0.8889132771304384e-3_wp, &
                                        0.9888347838921435e-3_wp,0.1053299681709471e-2_wp,0.1092778807014578e-2_wp, &
                                        0.1114389394063227e-2_wp,0.1123724788051555e-2_wp,0.1125239325243814e-2_wp, &
                                        0.1126153271815905e-2_wp,0.1130286931123841e-2_wp,0.1134986534363955e-2_wp], &
                              v5(N2) = [0.6823367927109931e-3_wp,0.9454158160447096e-3_wp,0.1074429975385679e-2_wp, &
                                        0.1129300086569132e-2_wp], &
                              v6(N3) = [0.8436884500901954e-3_wp,0.1075255720448885e-2_wp,0.1108577236864462e-2_wp, &
                                        0.9566475323783357e-3_wp,0.1080663250717391e-2_wp,0.1126797131196295e-2_wp, &
                                        0.1022568715358061e-2_wp,0.1108960267713108e-2_wp,0.1122790653435766e-2_wp, &
                                        0.1032401847117460e-2_wp,0.1107249382283854e-2_wp,0.1121780048519972e-2_wp]

  call ld_all(sym,order(26,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld0974

subroutine ld1046(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(27,sym)), y(order(27,sym)), z(order(27,sym)), w(order(27,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 12, N2 = 5, N3 = 13
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(27,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1046

subroutine ld1118(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(28,sym)), y(order(28,sym)), z(order(28,sym)), w(order(28,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 12, N2 = 6, N3 = 14
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(28,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1118

subroutine ld1202(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD1202 computes the 1202 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(29,sym)), y(order(29,sym)), z(order(29,sym)), w(order(29,sym))
  integer(kind=iwp), parameter :: N1 = 13, N2 = 4, N3 = 16
  real(kind=wp), parameter :: a4(N1) = [0.3712636449657089e-1_wp,0.9140060412262223e-1_wp,0.1531077852469906_wp, &
                                        0.2180928891660612_wp,0.2839874532200175_wp,0.3491177600963764_wp,0.4121431461444309_wp, &
                                        0.4718993627149127_wp,0.5273145452842337_wp,0.6209475332444019_wp,0.6569722711857291_wp, &
                                        0.6841788309070143_wp,0.7012604330123631_wp], &
                              a5(N2) = [0.1072382215478166_wp,0.2582068959496968_wp,0.4172752955306717_wp,0.5700366911792503_wp], &
                              a6(N3) = [0.9827986018263947_wp,0.9624249230326228_wp,0.9402007994128811_wp,0.9320822040143202_wp, &
                                        0.9043674199393299_wp,0.8912407560074747_wp,0.8676435628462708_wp,0.8581979986041619_wp, &
                                        0.8396753624049856_wp,0.8165288564022188_wp,0.8015469370783529_wp,0.7773563069070351_wp, &
                                        0.7661621213900394_wp,0.7553584143533510_wp,0.7344305757559503_wp,0.7043837184021765_wp], &
                              b6(N3) = [0.1771774022615325_wp,0.2475716463426288_wp,0.3354616289066489_wp,0.3173615246611977_wp, &
                                        0.4090268427085357_wp,0.3854291150669224_wp,0.4932221184851285_wp,0.4785320675922435_wp, &
                                        0.4507422593157064_wp,0.5632123020762100_wp,0.5434303569693900_wp,0.5123518486419871_wp, &
                                        0.6394279634749102_wp,0.6269805509024392_wp,0.6031161693096310_wp,0.5693702498468441_wp], &
                              v1 = 0.1105189233267572e-3_wp, v2 = 0.9205232738090741e-3_wp, v3 = 0.9133159786443561e-3_wp, &
                              v4(N1) = [0.3690421898017899e-3_wp,0.5603990928680660e-3_wp,0.6865297629282609e-3_wp, &
                                        0.7720338551145630e-3_wp,0.8301545958894795e-3_wp,0.8686692550179628e-3_wp, &
                                        0.8927076285846890e-3_wp,0.9060820238568219e-3_wp,0.9119777254940867e-3_wp, &
                                        0.9128720138604181e-3_wp,0.9130714935691735e-3_wp,0.9152873784554116e-3_wp, &
                                        0.9187436274321654e-3_wp], &
                              v5(N2) = [0.5176977312965694e-3_wp,0.7331143682101417e-3_wp,0.8463232836379928e-3_wp, &
                                        0.9031122694253992e-3_wp], &
                              v6(N3) = [0.6485778453163257e-3_wp,0.7435030910982369e-3_wp,0.7998527891839054e-3_wp, &
                                        0.8101731497468018e-3_wp,0.8483389574594331e-3_wp,0.8556299257311812e-3_wp, &
                                        0.8803208679738260e-3_wp,0.8811048182425720e-3_wp,0.8850282341265444e-3_wp, &
                                        0.9021342299040653e-3_wp,0.9010091677105086e-3_wp,0.9022692938426915e-3_wp, &
                                        0.9158016174693465e-3_wp,0.9131578003189435e-3_wp,0.9107813579482705e-3_wp, &
                                        0.9105760258970126e-3_wp]

  call ld_all(sym,order(29,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1202

subroutine ld1274(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(30,sym)), y(order(30,sym)), z(order(30,sym)), w(order(30,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 13, N2 = 7, N3 = 16
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(30,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1274

subroutine ld1358(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(31,sym)), y(order(31,sym)), z(order(31,sym)), w(order(31,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 14, N2 = 6, N3 = 18
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(31,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1358

subroutine ld1454(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD1454 computes the 1454 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(32,sym)), y(order(32,sym)), z(order(32,sym)), w(order(32,sym))
  integer(kind=iwp), parameter :: N1 = 15, N2 = 5, N3 = 20
  real(kind=wp), parameter :: a4(N1) = [0.3229290663413854e-1_wp,0.8036733271462222e-1_wp,0.1354289960531653_wp, &
                                        0.1938963861114426_wp,0.2537343715011275_wp,0.3135251434752570_wp,0.3721558339375338_wp, &
                                        0.4286809575195696_wp,0.4822510128282994_wp,0.5320679333566263_wp,0.6172998195394274_wp, &
                                        0.6510679849127481_wp,0.6777315251687360_wp,0.6963109410648741_wp,0.7058935009831749_wp], &
                              a5(N2) = [0.9955546194091857_wp,0.9734115901794209_wp,0.9275693732388626_wp,0.8568022422795103_wp, &
                                        0.7623495553719372_wp], &
                              a6(N3) = [0.5707522908892223_wp,0.5196463388403083_wp,0.4646337531215351_wp,0.4063901697557691_wp, &
                                        0.3456329466643087_wp,0.2831395121050332_wp,0.2197682022925330_wp,0.1564696098650355_wp, &
                                        0.6027356673721295_wp,0.5496032320255096_wp,0.4921707755234567_wp,0.4309422998598483_wp, &
                                        0.3664108182313672_wp,0.2990189057758436_wp,0.6268724013144998_wp,0.5707324144834607_wp, &
                                        0.5096360901960365_wp,0.4438729938312456_wp,0.6419978471082389_wp,0.5817218061802611_wp], &
                              b6(N3) = [0.4387028039889501_wp,0.3858908414762617_wp,0.3301937372343854_wp,0.2725423573563777_wp, &
                                        0.2139510237495250_wp,0.1555922309786647_wp,0.9892878979686097e-1_wp, &
                                        0.4598642910675510e-1_wp,0.3376625140173426_wp,0.2822301309727988_wp, &
                                        0.2248632342592540_wp,0.1666224723456479_wp,0.1086964901822169_wp, &
                                        0.5251989784120085e-1_wp,0.2297523657550023_wp,0.1723080607093800_wp, &
                                        0.1140238465390513_wp,0.5611522095882537e-1_wp,0.1164174423140873_wp, &
                                        0.5797589531445219e-1_wp], &
                              v1 = 0.7777160743261247e-4_wp, v2 = Zero, v3 = 0.7557646413004701e-3_wp, &
                              v4(N1) = [0.2841633806090617e-3_wp,0.4374419127053555e-3_wp,0.5417174740872172e-3_wp, &
                                        0.6148000891358593e-3_wp,0.6664394485800705e-3_wp,0.7025039356923220e-3_wp, &
                                        0.7268511789249627e-3_wp,0.7422637534208629e-3_wp,0.7509545035841214e-3_wp, &
                                        0.7548535057718401e-3_wp,0.7554088969774001e-3_wp,0.7553147174442808e-3_wp, &
                                        0.7564767653292297e-3_wp,0.7587991808518730e-3_wp,0.7608261832033027e-3_wp], &
                              v5(N2) = [0.4021680447874916e-3_wp,0.5804871793945964e-3_wp,0.6792151955945159e-3_wp, &
                                        0.7336741211286294e-3_wp,0.7581866300989608e-3_wp], &
                              v6(N3) = [0.7538257859800743e-3_wp,0.7483517247053123e-3_wp,0.7371763661112059e-3_wp, &
                                        0.7183448895756934e-3_wp,0.6895815529822191e-3_wp,0.6480105801792886e-3_wp, &
                                        0.5897558896594636e-3_wp,0.5095708849247346e-3_wp,0.7536906428909755e-3_wp, &
                                        0.7472505965575118e-3_wp,0.7343017132279698e-3_wp,0.7130871582177445e-3_wp, &
                                        0.6817022032112776e-3_wp,0.6380941145604121e-3_wp,0.7550381377920310e-3_wp, &
                                        0.7478646640144802e-3_wp,0.7335918720601220e-3_wp,0.7110120527658118e-3_wp, &
                                        0.7571363978689501e-3_wp,0.7489908329079234e-3_wp]

  call ld_all(sym,order(32,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1454

subroutine ld1538(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(33,sym)), y(order(33,sym)), z(order(33,sym)), w(order(33,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 15, N2 = 6, N3 = 21
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(33,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1538

subroutine ld1622(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(34,sym)), y(order(34,sym)), z(order(34,sym)), w(order(34,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 15, N2 = 8, N3 = 22
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(34,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1622

subroutine ld1730(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD1730 computes the 1730 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(35,sym)), y(order(35,sym)), z(order(35,sym)), w(order(35,sym))
  integer(kind=iwp), parameter :: N1 = 16, N2 = 5, N3 = 25
  real(kind=wp), parameter :: a4(N1) = [0.2860923126194662e-1_wp,0.7142556767711522e-1_wp,0.1209199540995559_wp, &
                                        0.1738673106594379_wp,0.2284645438467734_wp,0.2834807671701512_wp,0.3379680145467339_wp, &
                                        0.3911355454819537_wp,0.4422860353001403_wp,0.4907781568726057_wp,0.5360006153211468_wp, &
                                        0.6142105973596603_wp,0.6459300387977504_wp,0.6718056125089225_wp,0.6910888533186254_wp, &
                                        0.7030467416823252_wp], &
                              a5(N2) = [0.8354951166354646e-1_wp,0.2050143009099486_wp,0.3370208290706637_wp, &
                                        0.4689051484233963_wp,0.5939400424557334_wp], &
                              a6(N3) = [0.1394983311832261_wp,0.1967999180485014_wp,0.2546183732548967_wp,0.3121281074713875_wp, &
                                        0.3685981078502492_wp,0.4233760321547856_wp,0.4758671236059246_wp,0.5255178579796463_wp, &
                                        0.5718025633734589_wp,0.2686927772723415_wp,0.3306006819904809_wp,0.3904906850594983_wp, &
                                        0.4479957951904390_wp,0.5027076848919780_wp,0.5542087392260217_wp,0.6020850887375187_wp, &
                                        0.4019851409179594_wp,0.4635614567449800_wp,0.5215860931591575_wp,0.5758202499099271_wp, &
                                        0.6259893683876795_wp,0.5313795124811891_wp,0.5893317955931995_wp,0.6426246321215801_wp, &
                                        0.6511904367376113_wp], &
                              b6(N3) = [0.4097581162050343e-1_wp,0.8851987391293348e-1_wp,0.1397680182969819_wp, &
                                        0.1929452542226526_wp,0.2467898337061562_wp,0.3003104124785409_wp,0.3526684328175033_wp, &
                                        0.4031134861145713_wp,0.4509426448342351_wp,0.4711322502423248e-1_wp, &
                                        0.9784487303942695e-1_wp,0.1505395810025273_wp,0.2039728156296050_wp, &
                                        0.2571529941121107_wp,0.3092191375815670_wp,0.3593807506130276_wp, &
                                        0.5063389934378671e-1_wp,0.1032422269160612_wp,0.1566322094006254_wp, &
                                        0.2098082827491099_wp,0.2618824114553391_wp,0.5263245019338556e-1_wp, &
                                        0.1061059730982005_wp,0.1594171564034221_wp,0.5354789536565540e-1_wp], &
                              v1 = 0.6309049437420976e-4_wp, v2 = 0.6398287705571748e-3_wp, v3 = 0.6357185073530720e-3_wp, &
                              v4(N1) = [0.2221207162188168e-3_wp,0.3475784022286848e-3_wp,0.4350742443589804e-3_wp, &
                                        0.4978569136522127e-3_wp,0.5435036221998053e-3_wp,0.5765913388219542e-3_wp, &
                                        0.6001200359226003e-3_wp,0.6162178172717512e-3_wp,0.6265218152438485e-3_wp, &
                                        0.6323987160974212e-3_wp,0.6350767851540569e-3_wp,0.6354362775297107e-3_wp, &
                                        0.6352302462706235e-3_wp,0.6358117881417972e-3_wp,0.6373101590310117e-3_wp, &
                                        0.6390428961368665e-3_wp], &
                              v5(N2) = [0.3186913449946576e-3_wp,0.4678028558591711e-3_wp,0.5538829697598626e-3_wp, &
                                        0.6044475907190476e-3_wp,0.6313575103509012e-3_wp], &
                              v6(N3) = [0.4078626431855630e-3_wp,0.4759933057812725e-3_wp,0.5268151186413440e-3_wp, &
                                        0.5643048560507316e-3_wp,0.5914501076613073e-3_wp,0.6104561257874195e-3_wp, &
                                        0.6230252860707806e-3_wp,0.6305618761760796e-3_wp,0.6343092767597889e-3_wp, &
                                        0.5176268945737826e-3_wp,0.5564840313313692e-3_wp,0.5856426671038980e-3_wp, &
                                        0.6066386925777091e-3_wp,0.6208824962234458e-3_wp,0.6296314297822907e-3_wp, &
                                        0.6340423756791859e-3_wp,0.5829627677107342e-3_wp,0.6048693376081110e-3_wp, &
                                        0.6202362317732461e-3_wp,0.6299005328403779e-3_wp,0.6347722390609353e-3_wp, &
                                        0.6203778981238834e-3_wp,0.6308414671239979e-3_wp,0.6362706466959498e-3_wp, &
                                        0.6375414170333233e-3_wp]

  call ld_all(sym,order(35,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1730

subroutine ld1814(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(36,sym)), y(order(36,sym)), z(order(36,sym)), w(order(36,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 17, N2 = 8, N3 = 25
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(36,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1814

subroutine ld1910(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(37,sym)), y(order(37,sym)), z(order(37,sym)), w(order(37,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 17, N2 = 8, N3 = 27
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(37,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld1910

subroutine ld2030(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD2030 computes the 2030 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(38,sym)), y(order(38,sym)), z(order(38,sym)), w(order(38,sym))
  integer(kind=iwp), parameter :: N1 = 18, N2 = 6, N3 = 30
  real(kind=wp), parameter :: a4(N1) = [0.2540835336814348e-1_wp,0.6399322800504915e-1_wp,0.1088269469804125_wp, &
                                        0.1570670798818287_wp,0.2071163932282514_wp,0.2578914044450844_wp,0.3085687558169623_wp, &
                                        0.3584719706267024_wp,0.4070135594428709_wp,0.4536618626222638_wp,0.4979195686463577_wp, &
                                        0.5393075111126999_wp,0.6115617676843916_wp,0.6414308435160159_wp,0.6664099412721607_wp, &
                                        0.6859161771214913_wp,0.6993625593503890_wp,0.7062393387719380_wp], &
                              a5(N2) = [0.7479028168349763e-1_wp,0.1848951153969366_wp,0.3059529066581305_wp, &
                                        0.4285556101021362_wp,0.5468758653496526_wp,0.6565821978343439_wp], &
                              a6(N3) = [0.1253901572367117_wp,0.1775721510383941_wp,0.2305693358216114_wp,0.2836502845992063_wp, &
                                        0.3361794746232590_wp,0.3875979172264824_wp,0.4374019316999074_wp,0.4851275843340022_wp, &
                                        0.5303391803806868_wp,0.5726197380596287_wp,0.2431520732564863_wp,0.3002096800895869_wp, &
                                        0.3558554457457432_wp,0.4097782537048887_wp,0.4616337666067458_wp,0.5110707008417874_wp, &
                                        0.5577415286163795_wp,0.6013060431366950_wp,0.3661596767261781_wp,0.4237633153506581_wp, &
                                        0.4786328454658452_wp,0.5305702076789774_wp,0.5793436224231788_wp,0.6247069017094747_wp, &
                                        0.4874315552535204_wp,0.5427337322059053_wp,0.5943493747246700_wp,0.6421314033564943_wp, &
                                        0.6020628374713980_wp,0.6529222529856881_wp], &
                              b6(N3) = [0.3681917226439641e-1_wp,0.7982487607213301e-1_wp,0.1264640966592335_wp, &
                                        0.1751585683418957_wp,0.2247995907632670_wp,0.2745299257422246_wp,0.3236373482441118_wp, &
                                        0.3714967859436741_wp,0.4175353646321745_wp,0.4612084406355461_wp, &
                                        0.4258040133043952e-1_wp,0.8869424306722721e-1_wp,0.1368811706510655_wp, &
                                        0.1860739985015033_wp,0.2354235077395853_wp,0.2842074921347011_wp,0.3317784414984102_wp, &
                                        0.3775299002040700_wp,0.4599367887164592e-1_wp,0.9404893773654421e-1_wp, &
                                        0.1431377109091971_wp,0.1924186388843570_wp,0.2411590944775190_wp,0.2886871491583605_wp, &
                                        0.4804978774953206e-1_wp,0.9716857199366665e-1_wp,0.1465205839795055_wp, &
                                        0.1953579449803574_wp,0.4916375015738108e-1_wp,0.9861621540127005e-1_wp], &
                              v1 = 0.4656031899197431e-4_wp, v2 = Zero, v3 = 0.5421549195295507e-3_wp, &
                              v4(N1) = [0.1778522133346553e-3_wp,0.2811325405682796e-3_wp,0.3548896312631459e-3_wp, &
                                        0.4090310897173364e-3_wp,0.4493286134169965e-3_wp,0.4793728447962723e-3_wp, &
                                        0.5015415319164265e-3_wp,0.5175127372677937e-3_wp,0.5285522262081019e-3_wp, &
                                        0.5356832703713962e-3_wp,0.5397914736175170e-3_wp,0.5416899441599930e-3_wp, &
                                        0.5419308476889938e-3_wp,0.5416936902030596e-3_wp,0.5419544338703164e-3_wp, &
                                        0.5428983656630975e-3_wp,0.5442286500098193e-3_wp,0.5452250345057301e-3_wp], &
                              v5(N2) = [0.2568002497728530e-3_wp,0.3827211700292145e-3_wp,0.4579491561917824e-3_wp, &
                                        0.5042003969083574e-3_wp,0.5312708889976025e-3_wp,0.5438401790747117e-3_wp], &
                              v6(N3) = [0.3316041873197344e-3_wp,0.3899113567153771e-3_wp,0.4343343327201309e-3_wp, &
                                        0.4679415262318919e-3_wp,0.4930847981631031e-3_wp,0.5115031867540091e-3_wp, &
                                        0.5245217148457367e-3_wp,0.5332041499895321e-3_wp,0.5384583126021542e-3_wp, &
                                        0.5411067210798852e-3_wp,0.4259797391468714e-3_wp,0.4604931368460021e-3_wp, &
                                        0.4871814878255202e-3_wp,0.5072242910074885e-3_wp,0.5217069845235350e-3_wp, &
                                        0.5315785966280310e-3_wp,0.5376833708758905e-3_wp,0.5408032092069521e-3_wp, &
                                        0.4842744917904866e-3_wp,0.5048926076188130e-3_wp,0.5202607980478373e-3_wp, &
                                        0.5309932388325743e-3_wp,0.5377419770895208e-3_wp,0.5411696331677717e-3_wp, &
                                        0.5197996293282420e-3_wp,0.5311120836622945e-3_wp,0.5384309319956951e-3_wp, &
                                        0.5421859504051886e-3_wp,0.5390948355046314e-3_wp,0.5433312705027845e-3_wp]

  call ld_all(sym,order(38,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2030

subroutine ld2126(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(39,sym)), y(order(39,sym)), z(order(39,sym)), w(order(39,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 18, N2 = 8, N3 = 31
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(39,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2126

subroutine ld2222(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(40,sym)), y(order(40,sym)), z(order(40,sym)), w(order(40,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 18, N2 = 10, N3 = 32
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(40,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2222

subroutine ld2354(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD2354 computes the 2354 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(41,sym)), y(order(41,sym)), z(order(41,sym)), w(order(41,sym))
  integer(kind=iwp), parameter :: N1 = 19, N2 = 6, N3 = 36
  real(kind=wp), parameter :: a4(N1) = [0.2290024646530589e-1_wp,0.5779086652271284e-1_wp,0.9863103576375984e-1_wp, &
                                        0.1428155792982185_wp,0.1888978116601463_wp,0.2359091682970210_wp,0.2831228833706171_wp, &
                                        0.3299495857966693_wp,0.3758840802660796_wp,0.4204751831009480_wp,0.4633068518751051_wp, &
                                        0.5039849474507313_wp,0.5421265793440747_wp,0.6092660230557310_wp,0.6374654204984869_wp, &
                                        0.6615136472609892_wp,0.6809487285958127_wp,0.6952980021665196_wp,0.7041245497695400_wp], &
                              a5(N2) = [0.6744033088306065e-1_wp,0.1678684485334166_wp,0.2793559049539613_wp, &
                                        0.3935264218057639_wp,0.5052629268232558_wp,0.6107905315437531_wp], &
                              a6(N3) = [0.1135081039843524_wp,0.1612866626099378_wp,0.2100786550168205_wp,0.2592282009459942_wp, &
                                        0.3081740561320203_wp,0.3564289781578164_wp,0.4035587288240703_wp,0.4491671196373903_wp, &
                                        0.4928854782917489_wp,0.5343646791958988_wp,0.5732683216530990_wp,0.2214131583218986_wp, &
                                        0.2741796504750071_wp,0.3259797439149485_wp,0.3765441148826891_wp,0.4255773574530558_wp, &
                                        0.4727795117058430_wp,0.5178546895819012_wp,0.5605141192097460_wp,0.6004763319352512_wp, &
                                        0.3352842634946949_wp,0.3891971629814670_wp,0.4409875565542281_wp,0.4904893058592484_wp, &
                                        0.5375056138769549_wp,0.5818255708669969_wp,0.6232334858144959_wp,0.4489485354492058_wp, &
                                        0.5015136875933150_wp,0.5511300550512623_wp,0.5976720409858000_wp,0.6409956378989354_wp, &
                                        0.5581222330827514_wp,0.6074705984161695_wp,0.6532272537379033_wp,0.6594761494500487_wp], &
                              b6(N3) = [0.3331954884662588e-1_wp,0.7247167465436538e-1_wp,0.1151539110849745_wp, &
                                        0.1599491097143677_wp,0.2058699956028027_wp,0.2521624953502911_wp,0.2982090785797674_wp, &
                                        0.3434762087235733_wp,0.3874831357203437_wp,0.4297814821746926_wp,0.4699402260943537_wp, &
                                        0.3873602040643895e-1_wp,0.8089496256902013e-1_wp,0.1251732177620872_wp, &
                                        0.1706260286403185_wp,0.2165115147300408_wp,0.2622089812225259_wp,0.3071721431296201_wp, &
                                        0.3508998998801138_wp,0.3929160876166931_wp,0.4202563457288019e-1_wp, &
                                        0.8614309758870850e-1_wp,0.1314500879380001_wp,0.1772189657383859_wp, &
                                        0.2228277110050294_wp,0.2677179935014386_wp,0.3113675035544165_wp, &
                                        0.4409162378368174e-1_wp,0.8939009917748489e-1_wp,0.1351806029383365_wp, &
                                        0.1808370355053196_wp,0.2257852192301602_wp,0.4532173421637160e-1_wp, &
                                        0.9117488031840314e-1_wp,0.1369294213140155_wp,0.4589901487275583e-1_wp], &
                              v1 = 0.3922616270665292e-4_wp, v2 = 0.4703831750854424e-3_wp, v3 = 0.4678202801282136e-3_wp, &
                              v4(N1) = [0.1437832228979900e-3_wp,0.2303572493577644e-3_wp,0.2933110752447454e-3_wp, &
                                        0.3402905998359838e-3_wp,0.3759138466870372e-3_wp,0.4030638447899798e-3_wp, &
                                        0.4236591432242211e-3_wp,0.4390522656946746e-3_wp,0.4502523466626247e-3_wp, &
                                        0.4580577727783541e-3_wp,0.4631391616615899e-3_wp,0.4660928953698676e-3_wp, &
                                        0.4674751807936953e-3_wp,0.4676414903932920e-3_wp,0.4674086492347870e-3_wp, &
                                        0.4674928539483207e-3_wp,0.4680748979686447e-3_wp,0.4690449806389040e-3_wp, &
                                        0.4699877075860818e-3_wp], &
                              v5(N2) = [0.2099942281069176e-3_wp,0.3172269150712804e-3_wp,0.3832051358546523e-3_wp, &
                                        0.4252193818146985e-3_wp,0.4513807963755000e-3_wp,0.4657797469114178e-3_wp], &
                              v6(N3) = [0.2733362800522836e-3_wp,0.3235485368463559e-3_wp,0.3624908726013453e-3_wp, &
                                        0.3925540070712828e-3_wp,0.4156129781116235e-3_wp,0.4330644984623263e-3_wp, &
                                        0.4459677725921312e-3_wp,0.4551593004456795e-3_wp,0.4613341462749918e-3_wp, &
                                        0.4651019618269806e-3_wp,0.4670249536100625e-3_wp,0.3549555576441708e-3_wp, &
                                        0.3856108245249010e-3_wp,0.4098622845756882e-3_wp,0.4286328604268950e-3_wp, &
                                        0.4427802198993945e-3_wp,0.4530473511488561e-3_wp,0.4600805475703138e-3_wp, &
                                        0.4644599059958017e-3_wp,0.4667274455712508e-3_wp,0.4069360518020356e-3_wp, &
                                        0.4260442819919195e-3_wp,0.4408678508029063e-3_wp,0.4518748115548597e-3_wp, &
                                        0.4595564875375116e-3_wp,0.4643988774315846e-3_wp,0.4668827491646946e-3_wp, &
                                        0.4400541823741973e-3_wp,0.4514512890193797e-3_wp,0.4596198627347549e-3_wp, &
                                        0.4648659016801781e-3_wp,0.4675502017157673e-3_wp,0.4598494476455523e-3_wp, &
                                        0.4654916955152048e-3_wp,0.4684709779505137e-3_wp,0.4691445539106986e-3_wp]

  call ld_all(sym,order(41,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2354

subroutine ld2450(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(42,sym)), y(order(42,sym)), z(order(42,sym)), w(order(42,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 10, N2 = 10, N3 = 36
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(42,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2450

subroutine ld2558(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(43,sym)), y(order(43,sym)), z(order(43,sym)), w(order(43,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 20, N2 = 10, N3 = 38
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(43,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2558

subroutine ld2702(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD2702 computes the 2702 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(44,sym)), y(order(44,sym)), z(order(44,sym)), w(order(44,sym))
  integer(kind=iwp), parameter :: N1 = 21, N2 = 7, N3 = 42
  real(kind=wp), parameter :: a4(N1) = [0.2065562538818703e-1_wp,0.5250918173022379e-1_wp,0.8993480082038376e-1_wp, &
                                        0.1306023924436019_wp,0.1732060388531418_wp,0.2168727084820249_wp,0.2609528309173586_wp, &
                                        0.3049252927938952_wp,0.3483484138084404_wp,0.3908321549106406_wp,0.4320210071894814_wp, &
                                        0.4715824795890053_wp,0.5091984794078453_wp,0.5445580145650803_wp,0.6072575796841768_wp, &
                                        0.6339484505755803_wp,0.6570718257486958_wp,0.6762557330090709_wp,0.6911161696923790_wp, &
                                        0.7012841911659961_wp,0.7064559272410020_wp], &
                              a5(N2) = [0.6123554989894765e-1_wp,0.1533070348312393_wp,0.2563902605244206_wp, &
                                        0.3629346991663361_wp,0.4683949968987538_wp,0.5694479240657952_wp,0.6634465430993955_wp], &
                              a6(N3) = [0.1033958573552305_wp,0.1473521412414395_wp,0.1924552158705967_wp,0.2381094362890328_wp, &
                                        0.2838121707936760_wp,0.3291323133373415_wp,0.3736896978741460_wp,0.4171406040760013_wp, &
                                        0.4591677985256915_wp,0.4994733831718418_wp,0.5377731830445096_wp,0.5737917830001331_wp, &
                                        0.2027323586271389_wp,0.2516942375187273_wp,0.3000227995257181_wp,0.3474806691046342_wp, &
                                        0.3938103180359209_wp,0.4387519590455703_wp,0.4820503960077787_wp,0.5234573778475101_wp, &
                                        0.5627318647235282_wp,0.5996390607156954_wp,0.3084780753791947_wp,0.3589988275920223_wp, &
                                        0.4078628415881973_wp,0.4549287258889735_wp,0.5000278512957279_wp,0.5429785044928199_wp, &
                                        0.5835939850491711_wp,0.6216870353444856_wp,0.4151104662709091_wp,0.4649804275009218_wp, &
                                        0.5124695757009662_wp,0.5574711100606224_wp,0.5998597333287227_wp,0.6395007148516600_wp, &
                                        0.5188456224746252_wp,0.5664190707942778_wp,0.6110464353283153_wp,0.6526430302051563_wp, &
                                        0.6167551880377548_wp,0.6607195418355383_wp], &
                              b6(N3) = [0.3034544009063584e-1_wp,0.6618803044247135e-1_wp,0.1054431128987715_wp, &
                                        0.1468263551238858_wp,0.1894486108187886_wp,0.2326374238761579_wp,0.2758485808485768_wp, &
                                        0.3186179331996921_wp,0.3605329796303794_wp,0.4012147253586509_wp,0.4403050025570692_wp, &
                                        0.4774565904277483_wp,0.3544122504976147e-1_wp,0.7418304388646328e-1_wp, &
                                        0.1150502745727186_wp,0.1571963371209364_wp,0.1999631877247100_wp,0.2428073457846535_wp, &
                                        0.2852575132906155_wp,0.3268884208674639_wp,0.3673033321675939_wp,0.4061211551830290_wp, &
                                        0.3860125523100059e-1_wp,0.7928938987104867e-1_wp,0.1212614643030087_wp, &
                                        0.1638770827382693_wp,0.2065965798260176_wp,0.2489436378852235_wp,0.2904811368946891_wp, &
                                        0.3307941957666609_wp,0.4064829146052554e-1_wp,0.8258424547294755e-1_wp, &
                                        0.1251841962027289_wp,0.1679107505976331_wp,0.2102805057358715_wp,0.2518418087774107_wp, &
                                        0.4194321676077518e-1_wp,0.8457661551921499e-1_wp,0.1273652932519396_wp, &
                                        0.1698173239076354_wp,0.4266398851548864e-1_wp,0.8551925814238349e-1_wp], &
                              v1 = 0.2998675149888161e-4_wp, v2 = Zero, v3 = 0.4077860529495355e-3_wp, &
                              v4(N1) = [0.1185349192520667e-3_wp,0.1913408643425751e-3_wp,0.2452886577209897e-3_wp, &
                                        0.2862408183288702e-3_wp,0.3178032258257357e-3_wp,0.3422945667633690e-3_wp, &
                                        0.3612790520235922e-3_wp,0.3758638229818521e-3_wp,0.3868711798859953e-3_wp, &
                                        0.3949429933189938e-3_wp,0.4006068107541156e-3_wp,0.4043192149672723e-3_wp, &
                                        0.4064947495808078e-3_wp,0.4075245619813152e-3_wp,0.4076423540893566e-3_wp, &
                                        0.4074280862251555e-3_wp,0.4074163756012244e-3_wp,0.4077647795071246e-3_wp, &
                                        0.4084517552782530e-3_wp,0.4092468459224052e-3_wp,0.4097872687240906e-3_wp], &
                              v5(N2) = [0.1738986811745028e-3_wp,0.2659616045280191e-3_wp,0.3240596008171533e-3_wp, &
                                        0.3621195964432943e-3_wp,0.3868838330760539e-3_wp,0.4018911532693111e-3_wp, &
                                        0.4089929432983252e-3_wp], &
                              v6(N3) = [0.2279907527706409e-3_wp,0.2715205490578897e-3_wp,0.3057917896703976e-3_wp, &
                                        0.3326913052452555e-3_wp,0.3537334711890037e-3_wp,0.3700567500783129e-3_wp, &
                                        0.3825245372589122e-3_wp,0.3918125171518296e-3_wp,0.3984720419937579e-3_wp, &
                                        0.4029746003338211e-3_wp,0.4057428632156627e-3_wp,0.4071719274114857e-3_wp, &
                                        0.2990236950664119e-3_wp,0.3262951734212878e-3_wp,0.3482634608242413e-3_wp, &
                                        0.3656596681700892e-3_wp,0.3791740467794218e-3_wp,0.3894034450156905e-3_wp, &
                                        0.3968600245508371e-3_wp,0.4019931351420050e-3_wp,0.4052108801278599e-3_wp, &
                                        0.4068978613940934e-3_wp,0.3454275351319704e-3_wp,0.3629963537007920e-3_wp, &
                                        0.3770187233889873e-3_wp,0.3878608613694378e-3_wp,0.3959065270221274e-3_wp, &
                                        0.4015286975463570e-3_wp,0.4050866785614717e-3_wp,0.4069320185051913e-3_wp, &
                                        0.3760120964062763e-3_wp,0.3870969564418064e-3_wp,0.3955287790534055e-3_wp, &
                                        0.4015361911302668e-3_wp,0.4053836986719548e-3_wp,0.4073578673299117e-3_wp, &
                                        0.3954628379231406e-3_wp,0.4017645508847530e-3_wp,0.4059030348651293e-3_wp, &
                                        0.4080565809484880e-3_wp,0.4063018753664651e-3_wp,0.4087191292799671e-3_wp]

  call ld_all(sym,order(44,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2702

subroutine ld2810(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(45,sym)), y(order(45,sym)), z(order(45,sym)), w(order(45,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 21, N2 = 9, N3 = 43
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(45,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2810

subroutine ld2930(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(46,sym)), y(order(46,sym)), z(order(46,sym)), w(order(46,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 21, N2 = 10, N3 = 45
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(46,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld2930

subroutine ld3074(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD3074 computes the 3074 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(47,sym)), y(order(47,sym)), z(order(47,sym)), w(order(47,sym))
  integer(kind=iwp), parameter :: N1 = 22, N2 = 7, N3 = 49
  real(kind=wp), parameter :: a4(N1) = [0.1886108518723392e-1_wp,0.4800217244625303e-1_wp,0.8244922058397242e-1_wp, &
                                        0.1200408362484023_wp,0.1595773530809965_wp,0.2002635973434064_wp,0.2415127590139982_wp, &
                                        0.2828584158458477_wp,0.3239091015338138_wp,0.3643225097962194_wp,0.4037897083691802_wp, &
                                        0.4420247515194127_wp,0.4787572538464938_wp,0.5137265251275234_wp,0.5466764056654611_wp, &
                                        0.6054859420813535_wp,0.6308106701764562_wp,0.6530369230179584_wp,0.6718609524611158_wp, &
                                        0.6869676499894013_wp,0.6980467077240748_wp,0.7048241721250522_wp], &
                              a5(N2) = [0.5591105222058232e-1_wp,0.1407384078513916_wp,0.2364035438976309_wp, &
                                        0.3360602737818170_wp,0.4356292630054665_wp,0.5321569415256174_wp,0.6232956305040554_wp], &
                              a6(N3) = [0.9469870086838469e-1_wp,0.1353170300568141_wp,0.1771679481726077_wp, &
                                        0.2197066664231751_wp,0.2624783557374927_wp,0.3050969521214442_wp,0.3472252637196021_wp, &
                                        0.3885610219026360_wp,0.4288273776062765_wp,0.4677662471302948_wp,0.5051333589553359_wp, &
                                        0.5406942145810492_wp,0.5742204122576457_wp,0.1865407027225188_wp,0.2321186453689432_wp, &
                                        0.2773159142523882_wp,0.3219200192237254_wp,0.3657032593944029_wp,0.4084376778363622_wp, &
                                        0.4499004945751427_wp,0.4898758141326335_wp,0.5281547442266309_wp,0.5645346989813992_wp, &
                                        0.5988181252159848_wp,0.2850425424471603_wp,0.3324619433027876_wp,0.3785848333076282_wp, &
                                        0.4232891028562115_wp,0.4664287050829722_wp,0.5078458493735726_wp,0.5473779816204180_wp, &
                                        0.5848617133811376_wp,0.6201348281584888_wp,0.3852191185387871_wp,0.4325025061073423_wp, &
                                        0.4778486229734490_wp,0.5211663693009000_wp,0.5623469504853703_wp,0.6012718188659246_wp, &
                                        0.6378179206390117_wp,0.4836936460214534_wp,0.5293792562683797_wp,0.5726281253100033_wp, &
                                        0.6133658776169068_wp,0.6515085491865307_wp,0.5778692716064976_wp,0.6207904288086192_wp, &
                                        0.6608688171046802_wp,0.6656263089489130_wp], &
                              b6(N3) = [0.2778748387309470e-1_wp,0.6076569878628364e-1_wp,0.9703072762711040e-1_wp, &
                                        0.1354112458524762_wp,0.1750996479744100_wp,0.2154896907449802_wp,0.2560954625740152_wp, &
                                        0.2965070050624096_wp,0.3363641488734497_wp,0.3753400029836788_wp,0.4131297522144286_wp, &
                                        0.4494423776081795_wp,0.4839938958841502_wp,0.3259144851070796e-1_wp, &
                                        0.6835679505297343e-1_wp,0.1062284864451989_wp,0.1454404409323047_wp, &
                                        0.1854018282582510_wp,0.2256297412014750_wp,0.2657104425000896_wp,0.3052755487631557_wp, &
                                        0.3439863920645423_wp,0.3815229456121914_wp,0.4175752420966734_wp, &
                                        0.3562149509862536e-1_wp,0.7330318886871096e-1_wp,0.1123226296008472_wp, &
                                        0.1521084193337708_wp,0.1921844459223610_wp,0.2321360989678303_wp,0.2715886486360520_wp, &
                                        0.3101924707571355_wp,0.3476121052890973_wp,0.3763224880035108e-1_wp, &
                                        0.7659581935637135e-1_wp,0.1163381306083900_wp,0.1563890598752899_wp, &
                                        0.1963320810149200_wp,0.2357847407258738_wp,0.2743846121244060_wp, &
                                        0.3895902610739024e-1_wp,0.7871246819312640e-1_wp,0.1187963808202981_wp, &
                                        0.1587914708061787_wp,0.1983058575227646_wp,0.3977209689791542e-1_wp, &
                                        0.7990157592981152e-1_wp,0.1199671308754309_wp,0.4015955957805969e-1_wp], &
                              v1 = 0.2599095953754734e-4_wp, v2 = 0.3603134089687541e-3_wp, v3 = 0.3586067974412447e-3_wp, &
                              v4(N1) = [0.9831528474385880e-4_wp,0.1605023107954450e-3_wp,0.2072200131464099e-3_wp, &
                                        0.2431297618814187e-3_wp,0.2711819064496707e-3_wp,0.2932762038321116e-3_wp, &
                                        0.3107032514197368e-3_wp,0.3243808058921213e-3_wp,0.3349899091374030e-3_wp, &
                                        0.3430580688505218e-3_wp,0.3490124109290343e-3_wp,0.3532148948561955e-3_wp, &
                                        0.3559862669062833e-3_wp,0.3576224317551411e-3_wp,0.3584050533086076e-3_wp, &
                                        0.3584903581373224e-3_wp,0.3582991879040586e-3_wp,0.3582371187963125e-3_wp, &
                                        0.3584353631122350e-3_wp,0.3589120166517785e-3_wp,0.3595445704531601e-3_wp, &
                                        0.3600943557111074e-3_wp], &
                              v5(N2) = [0.1456447096742039e-3_wp,0.2252370188283782e-3_wp,0.2766135443474897e-3_wp, &
                                        0.3110729491500851e-3_wp,0.3342506712303391e-3_wp,0.3491981834026860e-3_wp, &
                                        0.3576003604348932e-3_wp], &
                              v6(N3) = [0.1921921305788564e-3_wp,0.2301458216495632e-3_wp,0.2604248549522893e-3_wp, &
                                        0.2845275425870697e-3_wp,0.3036870897974840e-3_wp,0.3188414832298066e-3_wp, &
                                        0.3307046414722089e-3_wp,0.3398330969031360e-3_wp,0.3466757899705373e-3_wp, &
                                        0.3516095923230054e-3_wp,0.3549645184048486e-3_wp,0.3570415969441392e-3_wp, &
                                        0.3581251798496118e-3_wp,0.2543491329913348e-3_wp,0.2786711051330776e-3_wp, &
                                        0.2985552361083679e-3_wp,0.3145867929154039e-3_wp,0.3273290662067609e-3_wp, &
                                        0.3372705511943501e-3_wp,0.3448274437851510e-3_wp,0.3503592783048583e-3_wp, &
                                        0.3541854792663162e-3_wp,0.3565995517909428e-3_wp,0.3578802078302898e-3_wp, &
                                        0.2958644592860982e-3_wp,0.3119548129116835e-3_wp,0.3250745225005984e-3_wp, &
                                        0.3355153415935208e-3_wp,0.3435847568549328e-3_wp,0.3495786831622488e-3_wp, &
                                        0.3537767805534621e-3_wp,0.3564459815421428e-3_wp,0.3578464061225468e-3_wp, &
                                        0.3239748762836212e-3_wp,0.3345491784174287e-3_wp,0.3429126177301782e-3_wp, &
                                        0.3492420343097421e-3_wp,0.3537399050235257e-3_wp,0.3566209152659172e-3_wp, &
                                        0.3581084321919782e-3_wp,0.3426522117591512e-3_wp,0.3491848770121379e-3_wp, &
                                        0.3539318235231476e-3_wp,0.3570231438458694e-3_wp,0.3586207335051714e-3_wp, &
                                        0.3541196205164025e-3_wp,0.3574296911573953e-3_wp,0.3591993279818963e-3_wp, &
                                        0.3595855034661997e-3_wp]

  call ld_all(sym,order(47,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3074

subroutine ld3182(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(48,sym)), y(order(48,sym)), z(order(48,sym)), w(order(48,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 23, N2 = 11, N3 = 49
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(48,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3182

subroutine ld3314(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(49,sym)), y(order(49,sym)), z(order(49,sym)), w(order(49,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 23, N2 = 10, N3 = 52
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(49,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3314

subroutine ld3470(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD3470 computes the 3470 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(50,sym)), y(order(50,sym)), z(order(50,sym)), w(order(50,sym))
  integer(kind=iwp), parameter :: N1 = 24, N2 = 8, N3 = 56
  real(kind=wp), parameter :: a4(N1) = [0.1721420832906233e-1_wp,0.4408875374981770e-1_wp,0.7594680813878681e-1_wp, &
                                        0.1108335359204799_wp,0.1476517054388567_wp,0.1856731870860615_wp,0.2243634099428821_wp, &
                                        0.2633006881662727_wp,0.3021340904916283_wp,0.3405594048030089_wp,0.3783044434007372_wp, &
                                        0.4151194767407910_wp,0.4507705766443257_wp,0.4850346056573187_wp,0.5176950817792470_wp, &
                                        0.5485384240820989_wp,0.6039117238943308_wp,0.6279956655573113_wp,0.6493636169568952_wp, &
                                        0.6677644117704504_wp,0.6829368572115624_wp,0.6946195818184121_wp,0.7025711542057026_wp, &
                                        0.7066004767140119_wp], &
                              a5(N2) = [0.5132537689946062e-1_wp,0.1297994661331225_wp,0.2188852049401307_wp, &
                                        0.3123174824903457_wp,0.4064037620738195_wp,0.4984958396944782_wp,0.5864975046021365_wp, &
                                        0.6686711634580175_wp], &
                              a6(N3) = [0.8715738780835950e-1_wp,0.1248383123134007_wp,0.1638062693383378_wp, &
                                        0.2035586203373176_wp,0.2436798975293774_wp,0.2838207507773806_wp,0.3236787502217692_wp, &
                                        0.3629849554840691_wp,0.4014948081992087_wp,0.4389818379260225_wp,0.4752331143674377_wp, &
                                        0.5100457318374018_wp,0.5432238388954868_wp,0.5745758685072442_wp,0.1723981437592809_wp, &
                                        0.2149553257844597_wp,0.2573256081247422_wp,0.2993163751238106_wp,0.3407238005148000_wp, &
                                        0.3813454978483264_wp,0.4209848104423343_wp,0.4594519699996300_wp,0.4965640166185930_wp, &
                                        0.5321441655571562_wp,0.5660208438582166_wp,0.5980264315964364_wp,0.2644215852350733_wp, &
                                        0.3090113743443063_wp,0.3525871079197808_wp,0.3950418005354029_wp,0.4362475663430163_wp, &
                                        0.4760661812145854_wp,0.5143551042512103_wp,0.5509709026935597_wp,0.5857711030329428_wp, &
                                        0.6186149917404392_wp,0.3586894569557064_wp,0.4035266610019441_wp,0.4467775312332510_wp, &
                                        0.4883638346608543_wp,0.5281908348434601_wp,0.5661542687149311_wp,0.6021450102031452_wp, &
                                        0.6360520783610050_wp,0.4521611065087196_wp,0.4959365651560963_wp,0.5376815804038283_wp, &
                                        0.5773314480243768_wp,0.6148113245575056_wp,0.6500407462842380_wp,0.5425151448707213_wp, &
                                        0.5841860556907931_wp,0.6234632186851500_wp,0.6602934551848843_wp,0.6278573968375105_wp, &
                                        0.6665611711264577_wp], &
                              b6(N3) = [0.2557175233367578e-1_wp,0.5604823383376681e-1_wp,0.8968568601900765e-1_wp, &
                                        0.1254086651976279_wp,0.1624780150162012_wp,0.2003422342683208_wp,0.2385628026255263_wp, &
                                        0.2767731148783578_wp,0.3146542308245309_wp,0.3519196415895088_wp,0.3883050984023654_wp, &
                                        0.4235613423908649_wp,0.4574484717196220_wp,0.4897311639255524_wp, &
                                        0.3010630597881105e-1_wp,0.6326031554204694e-1_wp,0.9848566980258631e-1_wp, &
                                        0.1350835952384266_wp,0.1725184055442181_wp,0.2103559279730725_wp,0.2482278774554860_wp, &
                                        0.2858099509982883_wp,0.3228075659915428_wp,0.3589459907204151_wp,0.3939630088864310_wp, &
                                        0.4276029922949089_wp,0.3300939429072552e-1_wp,0.6803887650078501e-1_wp, &
                                        0.1044326136206709_wp,0.1416751597517679_wp,0.1793408610504821_wp,0.2170630750175722_wp, &
                                        0.2545145157815807_wp,0.2913940101706601_wp,0.3274169910910705_wp,0.3623081329317265_wp, &
                                        0.3497354386450040e-1_wp,0.7129736739757095e-1_wp,0.1084758620193165_wp, &
                                        0.1460915689241772_wp,0.1837790832369980_wp,0.2212075390874021_wp,0.2580682841160985_wp, &
                                        0.2940656362094121_wp,0.3631055365867002e-1_wp,0.7348318468484350e-1_wp, &
                                        0.1111087643812648_wp,0.1488226085145408_wp,0.1862892274135151_wp,0.2231909701714456_wp, &
                                        0.3718201306118944e-1_wp,0.7483616335067346e-1_wp,0.1125990834266120_wp, &
                                        0.1501303813157619_wp,0.3767559930245720e-1_wp,0.7548443301360158e-1_wp], &
                              v1 = 0.2040382730826330e-4_wp, v2 = Zero, v3 = 0.3178149703889544e-3_wp, &
                              v4(N1) = [0.8288115128076110e-4_wp,0.1360883192522954e-3_wp,0.1766854454542662e-3_wp, &
                                        0.2083153161230153e-3_wp,0.2333279544657158e-3_wp,0.2532809539930247e-3_wp, &
                                        0.2692472184211158e-3_wp,0.2819949946811885e-3_wp,0.2920953593973030e-3_wp, &
                                        0.2999889782948352e-3_wp,0.3060292120496902e-3_wp,0.3105109167522192e-3_wp, &
                                        0.3136902387550312e-3_wp,0.3157984652454632e-3_wp,0.3170516518425422e-3_wp, &
                                        0.3176568425633755e-3_wp,0.3177198411207062e-3_wp,0.3175519492394733e-3_wp, &
                                        0.3174654952634756e-3_wp,0.3175676415467654e-3_wp,0.3178923417835410e-3_wp, &
                                        0.3183788287531909e-3_wp,0.3188755151918807e-3_wp,0.3191916889313849e-3_wp], &
                              v5(N2) = [0.1231779611744508e-3_wp,0.1924661373839880e-3_wp,0.2380881867403424e-3_wp, &
                                        0.2693100663037885e-3_wp,0.2908673382834366e-3_wp,0.3053914619381535e-3_wp, &
                                        0.3143916684147777e-3_wp,0.3187042244055363e-3_wp], &
                              v6(N3) = [0.1635219535869790e-3_wp,0.1968109917696070e-3_wp,0.2236754342249974e-3_wp, &
                                        0.2453186687017181e-3_wp,0.2627551791580541e-3_wp,0.2767654860152220e-3_wp, &
                                        0.2879467027765895e-3_wp,0.2967639918918702e-3_wp,0.3035900684660351e-3_wp, &
                                        0.3087338237298308e-3_wp,0.3124608838860167e-3_wp,0.3150084294226743e-3_wp, &
                                        0.3165958398598402e-3_wp,0.3174320440957372e-3_wp,0.2182188909812599e-3_wp, &
                                        0.2399727933921445e-3_wp,0.2579796133514652e-3_wp,0.2727114052623535e-3_wp, &
                                        0.2846327656281355e-3_wp,0.2941491102051334e-3_wp,0.3016049492136107e-3_wp, &
                                        0.3072949726175648e-3_wp,0.3114768142886460e-3_wp,0.3143823673666223e-3_wp, &
                                        0.3162269764661535e-3_wp,0.3172164663759821e-3_wp,0.2554575398967435e-3_wp, &
                                        0.2701704069135677e-3_wp,0.2823693413468940e-3_wp,0.2922898463214289e-3_wp, &
                                        0.3001829062162428e-3_wp,0.3062890864542953e-3_wp,0.3108328279264746e-3_wp, &
                                        0.3140243146201245e-3_wp,0.3160638030977130e-3_wp,0.3171462882206275e-3_wp, &
                                        0.2812388416031796e-3_wp,0.2912137500288045e-3_wp,0.2993241256502206e-3_wp, &
                                        0.3057101738983822e-3_wp,0.3105319326251432e-3_wp,0.3139565514428167e-3_wp, &
                                        0.3161543006806366e-3_wp,0.3172985960613294e-3_wp,0.2989400336901431e-3_wp, &
                                        0.3054555883947677e-3_wp,0.3104764960807702e-3_wp,0.3141015825977616e-3_wp, &
                                        0.3164520621159896e-3_wp,0.3176652305912204e-3_wp,0.3105097161023939e-3_wp, &
                                        0.3143014117890550e-3_wp,0.3168172866287200e-3_wp,0.3181401865570968e-3_wp, &
                                        0.3170663659156037e-3_wp,0.3185447944625510e-3_wp]

  call ld_all(sym,order(50,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3470

subroutine ld3590(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(51,sym)), y(order(51,sym)), z(order(51,sym)), w(order(51,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 24, N2 = 11, N3 = 57
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(51,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3590

subroutine ld3722(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(52,sym)), y(order(52,sym)), z(order(52,sym)), w(order(52,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 24, N2 = 12, N3 = 59
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(52,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3722

subroutine ld3890(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD3890 computes the 3890 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(53,sym)), y(order(53,sym)), z(order(53,sym)), w(order(53,sym))
  integer(kind=iwp), parameter :: N1 = 25, N2 = 8, N3 = 64
  real(kind=wp), parameter :: a4(N1) = [0.1587876419858352e-1_wp,0.4069193593751206e-1_wp,0.7025888115257997e-1_wp, &
                                        0.1027495450028704_wp,0.1371457730893426_wp,0.1727758532671953_wp,0.2091492038929037_wp, &
                                        0.2458813281751915_wp,0.2826545859450066_wp,0.3191957291799622_wp,0.3552621469299578_wp, &
                                        0.3906329503406230_wp,0.4251028614093031_wp,0.4584777520111870_wp,0.4905711358710193_wp, &
                                        0.5212011669847385_wp,0.5501878488737995_wp,0.6025037877479342_wp,0.6254572689549016_wp, &
                                        0.6460107179528248_wp,0.6639541138154251_wp,0.6790688515667495_wp,0.6911338580371512_wp, &
                                        0.6999385956126490_wp,0.7053037748656896_wp], &
                              a5(N2) = [0.4732224387180115e-1_wp,0.1202100529326803_wp,0.2034304820664855_wp, &
                                        0.2912285643573002_wp,0.3802361792726768_wp,0.4680598511056146_wp,0.5528151052155599_wp, &
                                        0.6329386307803041_wp], &
                              a6(N3) = [0.8056516651369069e-1_wp,0.1156476077139389_wp,0.1520473382760421_wp, &
                                        0.1892986699745931_wp,0.2270194446777792_wp,0.2648908185093273_wp,0.3026389259574136_wp, &
                                        0.3400220296151384_wp,0.3768217953335510_wp,0.4128372900921884_wp,0.4478807131815630_wp, &
                                        0.4817742034089257_wp,0.5143472814653344_wp,0.5454346213905650_wp,0.5748739313170252_wp, &
                                        0.1599598738286342_wp,0.1998097412500951_wp,0.2396228952566202_wp,0.2792228341097746_wp, &
                                        0.3184251107546741_wp,0.3570481164426244_wp,0.3949164710492144_wp,0.4318617293970503_wp, &
                                        0.4677221009931678_wp,0.5023417939270955_wp,0.5355701836636128_wp,0.5672608451328771_wp, &
                                        0.5972704202540162_wp,0.2461687022333596_wp,0.2881774566286831_wp,0.3293963604116978_wp, &
                                        0.3697303822241377_wp,0.4090663023135127_wp,0.4472819355411712_wp,0.4842513377231437_wp, &
                                        0.5198477629962928_wp,0.5539453011883145_wp,0.5864196762401251_wp,0.6171484466668390_wp, &
                                        0.3350337830565727_wp,0.3775773224758284_wp,0.4188155229848973_wp,0.4586805892009344_wp, &
                                        0.4970895714224235_wp,0.5339505133960747_wp,0.5691665792531440_wp,0.6026387682680377_wp, &
                                        0.6342676150163307_wp,0.4237951119537067_wp,0.4656918683234929_wp,0.5058857069185980_wp, &
                                        0.5443204666713996_wp,0.5809298813759742_wp,0.6156416039447128_wp,0.6483801351066604_wp, &
                                        0.5103616577251688_wp,0.5506738792580681_wp,0.5889573040995292_wp,0.6251641589516930_wp, &
                                        0.6592414921570178_wp,0.5930314017533384_wp,0.6309812253390175_wp,0.6666296011353230_wp, &
                                        0.6703715271049922_wp], &
                              b6(N3) = [0.2363454684003124e-1_wp,0.5191291632545936e-1_wp,0.8322715736994519e-1_wp, &
                                        0.1165855667993712_wp,0.1513077167409504_wp,0.1868882025807859_wp,0.2229277629776224_wp, &
                                        0.2590951840746235_wp,0.2951047291750847_wp,0.3307019714169930_wp,0.3656544101087634_wp, &
                                        0.3997448951939695_wp,0.4327667110812024_wp,0.4645196123532293_wp,0.4948063555703345_wp, &
                                        0.2792357590048985e-1_wp,0.5877141038139065e-1_wp,0.9164573914691377e-1_wp, &
                                        0.1259049641962687_wp,0.1610594823400863_wp,0.1967151653460898_wp,0.2325404606175168_wp, &
                                        0.2682461141151439_wp,0.3035720116011973_wp,0.3382781859197439_wp,0.3721383065625942_wp, &
                                        0.4049346360466055_wp,0.4364538098633802_wp,0.3070423166833368e-1_wp, &
                                        0.6338034669281885e-1_wp,0.9742862487067941e-1_wp,0.1323799532282290_wp, &
                                        0.1678497018129336_wp,0.2035095105326114_wp,0.2390692566672091_wp,0.2742649818076149_wp, &
                                        0.3088503806580094_wp,0.3425904245906614_wp,0.3752562294789468_wp, &
                                        0.3261589934634747e-1_wp,0.6658438928081572e-1_wp,0.1014565797157954_wp, &
                                        0.1368573320843822_wp,0.1724614851951608_wp,0.2079779381416412_wp,0.2431385788322288_wp, &
                                        0.2776901883049853_wp,0.3113881356386632_wp,0.3394877848664351e-1_wp, &
                                        0.6880219556291447e-1_wp,0.1041946859721635_wp,0.1398039738736393_wp, &
                                        0.1753373381196155_wp,0.2105215793514010_wp,0.2450953312157051_wp, &
                                        0.3485560643800719e-1_wp,0.7026308631512033e-1_wp,0.1059035061296403_wp, &
                                        0.1414823925236026_wp,0.1767207908214530_wp,0.3542189339561672e-1_wp, &
                                        0.7109574040369549e-1_wp,0.1067259792282730_wp,0.3569455268820809e-1_wp], &
                              v1 = 0.1807395252196920e-4_wp, v2 = 0.2848008782238827e-3_wp, v3 = 0.2836065837530581e-3_wp, &
                              v4(N1) = [0.7013149266673816e-4_wp,0.1162798021956766e-3_wp,0.1518728583972105e-3_wp, &
                                        0.1798796108216934e-3_wp,0.2022593385972785e-3_wp,0.2203093105575464e-3_wp, &
                                        0.2349294234299855e-3_wp,0.2467682058747003e-3_wp,0.2563092683572224e-3_wp, &
                                        0.2639253896763318e-3_wp,0.2699137479265108e-3_wp,0.2745196420166739e-3_wp, &
                                        0.2779529197397593e-3_wp,0.2803996086684265e-3_wp,0.2820302356715842e-3_wp, &
                                        0.2830056747491068e-3_wp,0.2834808950776839e-3_wp,0.2835282339078929e-3_wp, &
                                        0.2833819267065800e-3_wp,0.2832858336906784e-3_wp,0.2833268235451244e-3_wp, &
                                        0.2835432677029253e-3_wp,0.2839091722743049e-3_wp,0.2843308178875841e-3_wp, &
                                        0.2846703550533846e-3_wp], &
                              v5(N2) = [0.1051193406971900e-3_wp,0.1657871838796974e-3_wp,0.2064648113714232e-3_wp, &
                                        0.2347942745819741e-3_wp,0.2547775326597726e-3_wp,0.2686876684847025e-3_wp, &
                                        0.2778665755515867e-3_wp,0.2830996616782929e-3_wp], &
                              v6(N3) = [0.1403063340168372e-3_wp,0.1696504125939477e-3_wp,0.1935787242745390e-3_wp, &
                                        0.2130614510521968e-3_wp,0.2289381265931048e-3_wp,0.2418630292816186e-3_wp, &
                                        0.2523400495631193e-3_wp,0.2607623973449605e-3_wp,0.2674441032689209e-3_wp, &
                                        0.2726432360343356e-3_wp,0.2765787685924545e-3_wp,0.2794428690642224e-3_wp, &
                                        0.2814099002062895e-3_wp,0.2826429531578994e-3_wp,0.2832983542550884e-3_wp, &
                                        0.1886695565284976e-3_wp,0.2081867882748234e-3_wp,0.2245148680600796e-3_wp, &
                                        0.2380370491511872e-3_wp,0.2491398041852455e-3_wp,0.2581632405881230e-3_wp, &
                                        0.2653965506227417e-3_wp,0.2710857216747087e-3_wp,0.2754434093903659e-3_wp, &
                                        0.2786579932519380e-3_wp,0.2809011080679474e-3_wp,0.2823336184560987e-3_wp, &
                                        0.2831101175806309e-3_wp,0.2221679970354546e-3_wp,0.2356185734270703e-3_wp, &
                                        0.2469228344805590e-3_wp,0.2562726348642046e-3_wp,0.2638756726753028e-3_wp, &
                                        0.2699311157390862e-3_wp,0.2746233268403837e-3_wp,0.2781225674454771e-3_wp, &
                                        0.2805881254045684e-3_wp,0.2821719877004913e-3_wp,0.2830222502333124e-3_wp, &
                                        0.2457995956744870e-3_wp,0.2551474407503706e-3_wp,0.2629065335195311e-3_wp, &
                                        0.2691900449925075e-3_wp,0.2741275485754276e-3_wp,0.2778530970122595e-3_wp, &
                                        0.2805010567646741e-3_wp,0.2822055834031040e-3_wp,0.2831016901243473e-3_wp, &
                                        0.2624474901131803e-3_wp,0.2688034163039377e-3_wp,0.2738932751287636e-3_wp, &
                                        0.2777944791242523e-3_wp,0.2806011661660987e-3_wp,0.2824181456597460e-3_wp, &
                                        0.2833585216577828e-3_wp,0.2738165236962878e-3_wp,0.2778365208203180e-3_wp, &
                                        0.2807852940418966e-3_wp,0.2827245949674705e-3_wp,0.2837342344829828e-3_wp, &
                                        0.2809233907610981e-3_wp,0.2829930809742694e-3_wp,0.2841097874111479e-3_wp, &
                                        0.2843455206008783e-3_wp]

  call ld_all(sym,order(53,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld3890

subroutine ld4010(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(54,sym)), y(order(54,sym)), z(order(54,sym)), w(order(54,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 25, N2 = 13, N3 = 64
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(54,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4010

subroutine ld4154(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(55,sym)), y(order(55,sym)), z(order(55,sym)), w(order(55,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 26, N2 = 12, N3 = 67
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(55,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4154

subroutine ld4334(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD4334 computes the 4334 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(56,sym)), y(order(56,sym)), z(order(56,sym)), w(order(56,sym))
  integer(kind=iwp), parameter :: N1 = 27, N2 = 9, N3 = 72
  real(kind=wp), parameter :: a4(N1) = [0.1462896151831013e-1_wp,0.3769840812493139e-1_wp,0.6524701904096891e-1_wp, &
                                        0.9560543416134648e-1_wp,0.1278335898929198_wp,0.1613096104466031_wp, &
                                        0.1955806225745371_wp,0.2302935218498028_wp,0.2651584344113027_wp,0.2999276825183209_wp, &
                                        0.3343828669718798_wp,0.3683265013750518_wp,0.4015763206518108_wp,0.4339612026399770_wp, &
                                        0.4653180651114582_wp,0.4954893331080803_wp,0.5243207068924930_wp,0.5516590479041704_wp, &
                                        0.6012371927804176_wp,0.6231574466449819_wp,0.6429416514181271_wp,0.6604124272943595_wp, &
                                        0.6753851470408250_wp,0.6876717970626160_wp,0.6970895061319234_wp,0.7034746912553310_wp, &
                                        0.7067017217542295_wp], &
                              a5(N2) = [0.4382223501131123e-1_wp,0.1117474077400006_wp,0.1897153252911440_wp, &
                                        0.2724023009910331_wp,0.3567163308709902_wp,0.4404784483028087_wp,0.5219833154161411_wp, &
                                        0.5998179868977553_wp,0.6727803154548222_wp], &
                              a6(N3) = [0.7476563943166086e-1_wp,0.1075341482001416_wp,0.1416344885203259_wp, &
                                        0.1766325315388586_wp,0.2121744174481514_wp,0.2479669443408145_wp,0.2837600452294113_wp, &
                                        0.3193344933193984_wp,0.3544935442438745_wp,0.3890571932288154_wp,0.4228581214259090_wp, &
                                        0.4557387211304052_wp,0.4875487950541643_wp,0.5181436529962997_wp,0.5473824095600661_wp, &
                                        0.5751263398976174_wp,0.1489515746840028_wp,0.1863656444351767_wp,0.2238602880356348_wp, &
                                        0.2612723375728160_wp,0.2984332990206190_wp,0.3351786584663333_wp,0.3713505522209120_wp, &
                                        0.4067981098954663_wp,0.4413769993687534_wp,0.4749487182516394_wp,0.5073798105075426_wp, &
                                        0.5385410448878654_wp,0.5683065353670530_wp,0.5965527620663510_wp,0.2299227700856157_wp, &
                                        0.2695752998553267_wp,0.3086178716611389_wp,0.3469649871659077_wp,0.3845153566319655_wp, &
                                        0.4211600033403215_wp,0.4567867834329882_wp,0.4912829319232061_wp,0.5245364793303812_wp, &
                                        0.5564369788915756_wp,0.5868757697775287_wp,0.6157458853519617_wp,0.3138461110672113_wp, &
                                        0.3542495872050569_wp,0.3935751553120181_wp,0.4317634668111147_wp,0.4687413842250821_wp, &
                                        0.5044274237060283_wp,0.5387354077925727_wp,0.5715768898356105_wp,0.6028627200136111_wp, &
                                        0.6325039812653463_wp,0.3981986708423407_wp,0.4382791182133300_wp,0.4769233057218166_wp, &
                                        0.5140823911194238_wp,0.5496977833862983_wp,0.5837047306512727_wp,0.6160349566926879_wp, &
                                        0.6466185353209440_wp,0.4810835158795404_wp,0.5199925041324341_wp,0.5571717692207494_wp, &
                                        0.5925789250836378_wp,0.6261658523859670_wp,0.6578811126669331_wp,0.5609624612998100_wp, &
                                        0.5979959659984670_wp,0.6330523711054002_wp,0.6660960998103972_wp,0.6365384364585819_wp, &
                                        0.6710994302899275_wp], &
                              b6(N3) = [0.2193168509461185e-1_wp,0.4826419281533887e-1_wp,0.7751191883575742e-1_wp, &
                                        0.1087558139247680_wp,0.1413661374253096_wp,0.1748768214258880_wp,0.2089216406612073_wp, &
                                        0.2431987685545972_wp,0.2774497054377770_wp,0.3114460356156915_wp,0.3449806851913012_wp, &
                                        0.3778618641248256_wp,0.4099086391698978_wp,0.4409474925853973_wp,0.4708094517711291_wp, &
                                        0.4993275140354637_wp,0.2599381993267017e-1_wp,0.5479286532462190e-1_wp, &
                                        0.8556763251425254e-1_wp,0.1177257802267011_wp,0.1508168456192700_wp, &
                                        0.1844801892177727_wp,0.2184145236087598_wp,0.2523590641486229_wp,0.2860812976901373_wp, &
                                        0.3193686757808996_wp,0.3520226949547602_wp,0.3838544395667890_wp,0.4146810037640963_wp, &
                                        0.4443224094681121_wp,0.2865757664057584e-1_wp,0.5923421684485993e-1_wp, &
                                        0.9117817776057715e-1_wp,0.1240593814082605_wp,0.1575272058259175_wp, &
                                        0.1912845163525413_wp,0.2250710177858171_wp,0.2586521303440910_wp,0.2918112242865407_wp, &
                                        0.3243439239067890_wp,0.3560536787835351_wp,0.3867480821242581_wp, &
                                        0.3051374637507278e-1_wp,0.6237111233730755e-1_wp,0.9516223952401907e-1_wp, &
                                        0.1285467341508517_wp,0.1622318931656033_wp,0.1959581153836453_wp,0.2294888081183837_wp, &
                                        0.2626031152713945_wp,0.2950904075286713_wp,0.3267458451113286_wp, &
                                        0.3183291458749821e-1_wp,0.6459548193880908e-1_wp,0.9795757037087952e-1_wp, &
                                        0.1316307235126655_wp,0.1653556486358704_wp,0.1988931724126510_wp,0.2320174581438950_wp, &
                                        0.2645106562168662_wp,0.3275917807743992e-1_wp,0.6612546183967181e-1_wp, &
                                        0.9981498331474143e-1_wp,0.1335687001410374_wp,0.1671444402896463_wp, &
                                        0.2003106382156076_wp,0.3337500940231335e-1_wp,0.6708750335901803e-1_wp, &
                                        0.1008792126424850_wp,0.1345050343171794_wp,0.3372799460737052e-1_wp, &
                                        0.6755249309678028e-1_wp], &
                              v1 = 0.1449063022537883e-4_wp, v2 = Zero, v3 = 0.2546377329828424e-3_wp, &
                              v4(N1) = [0.6018432961087496e-4_wp,0.1002286583263673e-3_wp,0.1315222931028093e-3_wp, &
                                        0.1564213746876724e-3_wp,0.1765118841507736e-3_wp,0.1928737099311080e-3_wp, &
                                        0.2062658534263270e-3_wp,0.2172395445953787e-3_wp,0.2262076188876047e-3_wp, &
                                        0.2334885699462397e-3_wp,0.2393355273179203e-3_wp,0.2439559200468863e-3_wp, &
                                        0.2475251866060002e-3_wp,0.2501965558158773e-3_wp,0.2521081407925925e-3_wp, &
                                        0.2533881002388081e-3_wp,0.2541582900848261e-3_wp,0.2545365737525860e-3_wp, &
                                        0.2545726993066799e-3_wp,0.2544456197465555e-3_wp,0.2543481596881064e-3_wp, &
                                        0.2543506451429194e-3_wp,0.2544905675493763e-3_wp,0.2547611407344429e-3_wp, &
                                        0.2551060375448869e-3_wp,0.2554291933816039e-3_wp,0.2556255710686343e-3_wp], &
                              v5(N2) = [0.9041339695118195e-4_wp,0.1438426330079022e-3_wp,0.1802523089820518e-3_wp, &
                                        0.2060052290565496e-3_wp,0.2245002248967466e-3_wp,0.2377059847731150e-3_wp, &
                                        0.2468118955882525e-3_wp,0.2525410872966528e-3_wp,0.2553101409933397e-3_wp], &
                              v6(N3) = [0.1212879733668632e-3_wp,0.1472872881270931e-3_wp,0.1686846601010828e-3_wp, &
                                        0.1862698414660208e-3_wp,0.2007430956991861e-3_wp,0.2126568125394796e-3_wp, &
                                        0.2224394603372113e-3_wp,0.2304264522673135e-3_wp,0.2368854288424087e-3_wp, &
                                        0.2420352089461772e-3_wp,0.2460597113081295e-3_wp,0.2491181912257687e-3_wp, &
                                        0.2513528194205857e-3_wp,0.2528943096693220e-3_wp,0.2538660368488136e-3_wp, &
                                        0.2543868648299022e-3_wp,0.1642595537825183e-3_wp,0.1818246659849308e-3_wp, &
                                        0.1966565649492420e-3_wp,0.2090677905657991e-3_wp,0.2193820409510504e-3_wp, &
                                        0.2278870827661928e-3_wp,0.2348283192282090e-3_wp,0.2404139755581477e-3_wp, &
                                        0.2448227407760734e-3_wp,0.2482110455592573e-3_wp,0.2507192397774103e-3_wp, &
                                        0.2524765968534880e-3_wp,0.2536052388539425e-3_wp,0.2542230588033068e-3_wp, &
                                        0.1944817013047896e-3_wp,0.2067862362746635e-3_wp,0.2172440734649114e-3_wp, &
                                        0.2260125991723423e-3_wp,0.2332655008689523e-3_wp,0.2391699681532458e-3_wp, &
                                        0.2438801528273928e-3_wp,0.2475370504260665e-3_wp,0.2502707235640574e-3_wp, &
                                        0.2522031701054241e-3_wp,0.2534511269978784e-3_wp,0.2541284914955151e-3_wp, &
                                        0.2161509250688394e-3_wp,0.2248778513437852e-3_wp,0.2322388803404617e-3_wp, &
                                        0.2383265471001355e-3_wp,0.2432476675019525e-3_wp,0.2471122223750674e-3_wp, &
                                        0.2500291752486870e-3_wp,0.2521055942764682e-3_wp,0.2534472785575503e-3_wp, &
                                        0.2541599713080121e-3_wp,0.2317380975862936e-3_wp,0.2378550733719775e-3_wp, &
                                        0.2428884456739118e-3_wp,0.2469002655757292e-3_wp,0.2499657574265851e-3_wp, &
                                        0.2521676168486082e-3_wp,0.2535935662645334e-3_wp,0.2543356743363214e-3_wp, &
                                        0.2427353285201535e-3_wp,0.2468258039744386e-3_wp,0.2500060956440310e-3_wp, &
                                        0.2523238365420979e-3_wp,0.2538399260252846e-3_wp,0.2546255927268069e-3_wp, &
                                        0.2500583360048449e-3_wp,0.2524777638260203e-3_wp,0.2540951193860656e-3_wp, &
                                        0.2549524085027472e-3_wp,0.2542569507009158e-3_wp,0.2552114127580376e-3_wp]

  call ld_all(sym,order(56,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4334

subroutine ld4466(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(57,sym)), y(order(57,sym)), z(order(57,sym)), w(order(57,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 27, N2 = 12, N3 = 73
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(57,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4466

subroutine ld4610(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(58,sym)), y(order(58,sym)), z(order(58,sym)), w(order(58,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 27, N2 = 14, N3 = 75
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(58,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4610

subroutine ld4802(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD4802 computes the 4802 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(59,sym)), y(order(59,sym)), z(order(59,sym)), w(order(59,sym))
  integer(kind=iwp), parameter :: N1 = 28, N2 = 9, N3 = 81
  real(kind=wp), parameter :: a4(N1) = [0.2335728608887064e-1_wp,0.4352987836550653e-1_wp,0.6439200521088801e-1_wp, &
                                        0.9003943631993181e-1_wp,0.1196706615548473_wp,0.1511715412838134_wp, &
                                        0.1835982828503801_wp,0.2165081259155405_wp,0.2496208720417563_wp,0.2827200673567900_wp, &
                                        0.3156190823994346_wp,0.3481476793749115_wp,0.3801466086947226_wp,0.4114652119634011_wp, &
                                        0.4419598786519751_wp,0.4714925949329543_wp,0.4999293972879466_wp,0.5271387221431248_wp, &
                                        0.5529896780837761_wp,0.6000856099481712_wp,0.6210562192785175_wp,0.6401165879934240_wp, &
                                        0.6571144029244334_wp,0.6718910821718863_wp,0.6842845591099010_wp,0.6941353476269816_wp, &
                                        0.7012965242212991_wp,0.7056471428242644_wp], &
                              a5(N2) = [0.4595557643585895e-1_wp,0.1049316742435023_wp,0.1773548879549274_wp, &
                                        0.2559071411236127_wp,0.3358156837985898_wp,0.4155835743763893_wp,0.4937894296167472_wp, &
                                        0.5691569694793316_wp,0.6405840854894251_wp], &
                              a6(N3) = [0.7345133894143348e-1_wp,0.1009859834044931_wp,0.1324289619748758_wp, &
                                        0.1654272109607127_wp,0.1990767186776461_wp,0.2330125945523278_wp,0.2670080611108287_wp, &
                                        0.3008753376294316_wp,0.3344475596167860_wp,0.3675709724070786_wp,0.4001000887587812_wp, &
                                        0.4318956350436028_wp,0.4628239056795531_wp,0.4927563229773636_wp,0.5215687136707969_wp, &
                                        0.5491402346984905_wp,0.5753520160126075_wp,0.1388326356417754_wp,0.1743686900537244_wp, &
                                        0.2099737037950268_wp,0.2454492590908548_wp,0.2807219257864278_wp,0.3156842271975842_wp, &
                                        0.3502090945177752_wp,0.3841684849519686_wp,0.4174372367906016_wp,0.4498926465011892_wp, &
                                        0.4814146229807701_wp,0.5118863625734701_wp,0.5411947455119144_wp,0.5692301500357246_wp, &
                                        0.5958857204139576_wp,0.2156270284785766_wp,0.2532385054909710_wp,0.2902564617771537_wp, &
                                        0.3266979823143256_wp,0.3625039627493614_wp,0.3975838937548699_wp,0.4318396099009774_wp, &
                                        0.4651706555732742_wp,0.4974752649620969_wp,0.5286517579627517_wp,0.5586001195731895_wp, &
                                        0.5872229902021319_wp,0.6144258616235123_wp,0.2951676508064861_wp,0.3335085485472725_wp, &
                                        0.3709561760636381_wp,0.4074722861667498_wp,0.4429923648839117_wp,0.4774428052721736_wp, &
                                        0.5107446539535904_wp,0.5428151370542935_wp,0.5735699292556964_wp,0.6029253794562866_wp, &
                                        0.6307998987073145_wp,0.3752652273692719_wp,0.4135383879344028_wp,0.4506113885153907_wp, &
                                        0.4864401554606072_wp,0.5209708076611709_wp,0.5541422135830122_wp,0.5858880915113817_wp, &
                                        0.6161399390603444_wp,0.6448296482255090_wp,0.4544796274917948_wp,0.4919389072146628_wp, &
                                        0.5279313026985183_wp,0.5624169925571135_wp,0.5953484627093287_wp,0.6266730715339185_wp, &
                                        0.6563363204278871_wp,0.5314574716585696_wp,0.5674614932298185_wp,0.6017706004970264_wp, &
                                        0.6343471270264178_wp,0.6651494599127802_wp,0.6050184986005704_wp,0.6390163550880400_wp, &
                                        0.6711199107088448_wp,0.6741354429572275_wp], &
                              b6(N3) = [0.2177844081486067e-1_wp,0.4590362185775188e-1_wp,0.7255063095690877e-1_wp, &
                                        0.1017825451960684_wp,0.1325652320980364_wp,0.1642765374496765_wp,0.1965360374337889_wp, &
                                        0.2290726770542238_wp,0.2616645495370823_wp,0.2941150728843141_wp,0.3262440400919066_wp, &
                                        0.3578835350611916_wp,0.3888751854043678_wp,0.4190678003222840_wp,0.4483151836883852_wp, &
                                        0.4764740676087880_wp,0.5034021310998277_wp,0.2435436510372806e-1_wp, &
                                        0.5118897057342652e-1_wp,0.8014695048539634e-1_wp,0.1105117874155699_wp, &
                                        0.1417950531570966_wp,0.1736604945719597_wp,0.2058466324693981_wp,0.2381284261195919_wp, &
                                        0.2703031270422569_wp,0.3021845683091309_wp,0.3335993355165720_wp,0.3643833735518232_wp, &
                                        0.3943789541958179_wp,0.4234320144403542_wp,0.4513897947419260_wp, &
                                        0.2681225755444491e-1_wp,0.5557495747805614e-1_wp,0.8569368062950249e-1_wp, &
                                        0.1167367450324135_wp,0.1483861994003304_wp,0.1803821503011405_wp,0.2124962965666424_wp, &
                                        0.2445221837805913_wp,0.2762701224322987_wp,0.3075627775211328_wp,0.3382311089826877_wp, &
                                        0.3681108834741399_wp,0.3970397446872839_wp,0.2867499538750441e-1_wp, &
                                        0.5867879341903510e-1_wp,0.8961099205022284e-1_wp,0.1211627927626297_wp, &
                                        0.1530748903554898_wp,0.1851176436721877_wp,0.2170829107658179_wp,0.2487786689026271_wp, &
                                        0.2800239952795016_wp,0.3106445702878119_wp,0.3404689500841194_wp, &
                                        0.2997145098184479e-1_wp,0.6086725898678011e-1_wp,0.9238849548435643e-1_wp, &
                                        0.1242786603851851_wp,0.1563086731483386_wp,0.1882696509388506_wp,0.2199672979126059_wp, &
                                        0.2512165482924867_wp,0.2818368701871888_wp,0.3088970405060312e-1_wp, &
                                        0.6240947677636835e-1_wp,0.9430706144280313e-1_wp,0.1263547818770374_wp, &
                                        0.1583430788822594_wp,0.1900748462555988_wp,0.2213599519592567_wp, &
                                        0.3152508811515374e-1_wp,0.6343865291465561e-1_wp,0.9551503504223951e-1_wp, &
                                        0.1275440099801196_wp,0.1593252037671960_wp,0.3192538338496105e-1_wp, &
                                        0.6402824353962306e-1_wp,0.9609805077002909e-1_wp,0.3211853196273233e-1_wp], &
                              v1 = 0.9687521879420705e-4_wp, v2 = 0.2307897895367918e-3_wp, v3 = 0.2297310852498558e-3_wp, &
                              v4(N1) = [0.7386265944001919e-4_wp,0.8257977698542210e-4_wp,0.9706044762057630e-4_wp, &
                                        0.1302393847117003e-3_wp,0.1541957004600968e-3_wp,0.1704459770092199e-3_wp, &
                                        0.1827374890942906e-3_wp,0.1926360817436107e-3_wp,0.2008010239494833e-3_wp, &
                                        0.2075635983209175e-3_wp,0.2131306638690909e-3_wp,0.2176562329937335e-3_wp, &
                                        0.2212682262991018e-3_wp,0.2240799515668565e-3_wp,0.2261959816187525e-3_wp, &
                                        0.2277156368808855e-3_wp,0.2287351772128336e-3_wp,0.2293490814084085e-3_wp, &
                                        0.2296505312376273e-3_wp,0.2296793832318756e-3_wp,0.2295785443842974e-3_wp, &
                                        0.2295017931529102e-3_wp,0.2295059638184868e-3_wp,0.2296232343237362e-3_wp, &
                                        0.2298530178740771e-3_wp,0.2301579790280501e-3_wp,0.2304690404996513e-3_wp, &
                                        0.2307027995907102e-3_wp], &
                              v5(N2) = [0.9312274696671092e-4_wp,0.1199919385876926e-3_wp,0.1598039138877690e-3_wp, &
                                        0.1822253763574900e-3_wp,0.1988579593655040e-3_wp,0.2112620102533307e-3_wp, &
                                        0.2201594887699007e-3_wp,0.2261622590895036e-3_wp,0.2296458453435705e-3_wp], &
                              v6(N3) = [0.1006006990267000e-3_wp,0.1227676689635876e-3_wp,0.1467864280270117e-3_wp, &
                                        0.1644178912101232e-3_wp,0.1777664890718961e-3_wp,0.1884825664516690e-3_wp, &
                                        0.1973269246453848e-3_wp,0.2046767775855328e-3_wp,0.2107600125918040e-3_wp, &
                                        0.2157416362266829e-3_wp,0.2197557816920721e-3_wp,0.2229192611835437e-3_wp, &
                                        0.2253385110212775e-3_wp,0.2271137107548774e-3_wp,0.2283414092917525e-3_wp, &
                                        0.2291161673130077e-3_wp,0.2295313908576598e-3_wp,0.1438204721359031e-3_wp, &
                                        0.1607738025495257e-3_wp,0.1741483853528379e-3_wp,0.1851918467519151e-3_wp, &
                                        0.1944628638070613e-3_wp,0.2022495446275152e-3_wp,0.2087462382438514e-3_wp, &
                                        0.2141074754818308e-3_wp,0.2184640913748162e-3_wp,0.2219309165220329e-3_wp, &
                                        0.2246123118340624e-3_wp,0.2266062766915125e-3_wp,0.2280072952230796e-3_wp, &
                                        0.2289082025202583e-3_wp,0.2294012695120025e-3_wp,0.1722434488736947e-3_wp, &
                                        0.1830237421455091e-3_wp,0.1923855349997633e-3_wp,0.2004067861936271e-3_wp, &
                                        0.2071817297354263e-3_wp,0.2128250834102103e-3_wp,0.2174513719440102e-3_wp, &
                                        0.2211661839150214e-3_wp,0.2240665257813102e-3_wp,0.2262439516632620e-3_wp, &
                                        0.2277874557231869e-3_wp,0.2287854314454994e-3_wp,0.2293268499615575e-3_wp, &
                                        0.1912628201529828e-3_wp,0.1992499672238701e-3_wp,0.2061275533454027e-3_wp, &
                                        0.2119318215968572e-3_wp,0.2167416581882652e-3_wp,0.2206430730516600e-3_wp, &
                                        0.2237186938699523e-3_wp,0.2260480075032884e-3_wp,0.2277098884558542e-3_wp, &
                                        0.2287845715109671e-3_wp,0.2293547268236294e-3_wp,0.2056073839852528e-3_wp, &
                                        0.2114235865831876e-3_wp,0.2163175629770551e-3_wp,0.2203392158111650e-3_wp, &
                                        0.2235473176847839e-3_wp,0.2260024141501235e-3_wp,0.2277675929329182e-3_wp, &
                                        0.2289102112284834e-3_wp,0.2295027954625118e-3_wp,0.2161281589879992e-3_wp, &
                                        0.2201980477395102e-3_wp,0.2234952066593166e-3_wp,0.2260540098520838e-3_wp, &
                                        0.2279157981899988e-3_wp,0.2291296918565571e-3_wp,0.2297533752536649e-3_wp, &
                                        0.2234927356465995e-3_wp,0.2261288012985219e-3_wp,0.2280818160923688e-3_wp, &
                                        0.2293773295180159e-3_wp,0.2300528767338634e-3_wp,0.2281893855065666e-3_wp, &
                                        0.2295720444840727e-3_wp,0.2303227649026753e-3_wp,0.2304831913227114e-3_wp]

  call ld_all(sym,order(59,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4802

subroutine ld4934(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(60,sym)), y(order(60,sym)), z(order(60,sym)), w(order(60,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 29, N2 = 14, N3 = 81
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(60,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld4934

subroutine ld5090(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(61,sym)), y(order(61,sym)), z(order(61,sym)), w(order(61,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 29, N2 = 14, N3 = 84
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(61,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld5090

subroutine ld5294(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD5294 computes the 5294 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(62,sym)), y(order(62,sym)), z(order(62,sym)), w(order(62,sym))
  integer(kind=iwp), parameter :: N1 = 30, N2 = 10, N3 = 90
  real(kind=wp), parameter :: a4(N1) = [0.2303261686261450e-1_wp,0.3757208620162394e-1_wp,0.5821912033821852e-1_wp, &
                                        0.8403127529194872e-1_wp,0.1122927798060578_wp,0.1420125319192987_wp, &
                                        0.1726396437341978_wp,0.2038170058115696_wp,0.2352849892876508_wp,0.2668363354312461_wp, &
                                        0.2982941279900452_wp,0.3295002922087076_wp,0.3603094918363593_wp,0.3905857895173920_wp, &
                                        0.4202005758160837_wp,0.4490310061597227_wp,0.4769586160311491_wp,0.5038679887049750_wp, &
                                        0.5296454286519961_wp,0.5541776207164850_wp,0.5990467321921213_wp,0.6191467096294587_wp, &
                                        0.6375251212901849_wp,0.6540514381131168_wp,0.6685899064391510_wp,0.6810013009681648_wp, &
                                        0.6911469578730340_wp,0.6988956915141736_wp,0.7041335794868720_wp,0.7067754398018567_wp], &
                              a5(N2) = [0.3840368707853623e-1_wp,0.9835485954117399e-1_wp,0.1665774947612998_wp, &
                                        0.2405702335362910_wp,0.3165270770189046_wp,0.3927386145645443_wp,0.4678825918374656_wp, &
                                        0.5408022024266935_wp,0.6104967445752438_wp,0.6760910702685738_wp], &
                              a6(N3) = [0.6655644120217392e-1_wp,0.9446246161270182e-1_wp,0.1242651925452509_wp, &
                                        0.1553438064846751_wp,0.1871137110542670_wp,0.2192612628836257_wp,0.2515682807206955_wp, &
                                        0.2838535866287290_wp,0.3159578817528521_wp,0.3477370882791392_wp,0.3790576960890540_wp, &
                                        0.4097938317810200_wp,0.4398256572859637_wp,0.4690384114718480_wp,0.4973216048301053_wp, &
                                        0.5245681526132446_wp,0.5506733911803888_wp,0.5755339829522475_wp,0.1305472386056362_wp, &
                                        0.1637327908216477_wp,0.1972734634149637_wp,0.2308694653110130_wp,0.2643899218338160_wp, &
                                        0.2977171599622171_wp,0.3307293903032310_wp,0.3633069198219073_wp,0.3953346955922727_wp, &
                                        0.4267018394184914_wp,0.4573009622571704_wp,0.4870279559856109_wp,0.5157819581450322_wp, &
                                        0.5434651666465393_wp,0.5699823887764627_wp,0.5952403350947741_wp,0.2025152599210369_wp, &
                                        0.2381066653274425_wp,0.2732823383651612_wp,0.3080137692611118_wp,0.3422405614587601_wp, &
                                        0.3758808773890420_wp,0.4088458383438932_wp,0.4410450550841152_wp,0.4723879420561312_wp, &
                                        0.5027843561874343_wp,0.5321453674452458_wp,0.5603839113834030_wp,0.5874150706875146_wp, &
                                        0.6131559381660038_wp,0.2778497016394506_wp,0.3143733562261912_wp,0.3501485810261827_wp, &
                                        0.3851430322303653_wp,0.4193013979470415_wp,0.4525585960458567_wp,0.4848447779622947_wp, &
                                        0.5160871208276894_wp,0.5462112185696926_wp,0.5751425068101757_wp,0.6028073872853596_wp, &
                                        0.6291338275278409_wp,0.3541797528439391_wp,0.3908234972074657_wp,0.4264408450107590_wp, &
                                        0.4609949666553286_wp,0.4944389496536006_wp,0.5267194884346086_wp,0.5577787810220990_wp, &
                                        0.5875563763536670_wp,0.6159910016391269_wp,0.6430219602956268_wp,0.4300647036213646_wp, &
                                        0.4661486308935531_wp,0.5009658555287261_wp,0.5344824270447704_wp,0.5666575997416371_wp, &
                                        0.5974457471404752_wp,0.6267984444116886_wp,0.6546664713575417_wp,0.5042711004437253_wp, &
                                        0.5392127456774380_wp,0.5726819437668618_wp,0.6046469254207278_wp,0.6350716157434952_wp, &
                                        0.6639177679185454_wp,0.5757276040972253_wp,0.6090265823139755_wp,0.6406735344387661_wp, &
                                        0.6706397927793709_wp,0.6435019674426665_wp,0.6747218676375681_wp], &
                              b6(N3) = [0.1936508874588424e-1_wp,0.4252442002115869e-1_wp,0.6806529315354374e-1_wp, &
                                        0.9560957491205369e-1_wp,0.1245931657452888_wp,0.1545385828778978_wp, &
                                        0.1851004249723368_wp,0.2160182608272384_wp,0.2470799012277111_wp,0.2781014208986402_wp, &
                                        0.3089172523515731_wp,0.3393750055472244_wp,0.3693322470987730_wp,0.3986541005609877_wp, &
                                        0.4272112491408562_wp,0.4548781735309936_wp,0.4815315355023251_wp,0.5070486445801855_wp, &
                                        0.2284970375722366e-1_wp,0.4812254338288384e-1_wp,0.7531734457511935e-1_wp, &
                                        0.1039043639882017_wp,0.1334526587117626_wp,0.1636414868936382_wp,0.1942195406166568_wp, &
                                        0.2249752879943753_wp,0.2557218821820032_wp,0.2862897925213193_wp,0.3165224536636518_wp, &
                                        0.3462730221636496_wp,0.3754016870282835_wp,0.4037733784993613_wp,0.4312557784139123_wp, &
                                        0.4577175367122110_wp,0.2520253617719557e-1_wp,0.5223254506119000e-1_wp, &
                                        0.8060669688588620e-1_wp,0.1099335754081255_wp,0.1399120955959857_wp, &
                                        0.1702977801651705_wp,0.2008799256601680_wp,0.2314703052180836_wp,0.2618972111375892_wp, &
                                        0.2920013195600270_wp,0.3216322555190551_wp,0.3506456615934198_wp,0.3789007181306267_wp, &
                                        0.4062580170572782_wp,0.2696271276876226e-1_wp,0.5523469316960465e-1_wp, &
                                        0.8445193201626464e-1_wp,0.1143263119336083_wp,0.1446177898344475_wp, &
                                        0.1751165438438091_wp,0.2056338306745660_wp,0.2359965487229226_wp,0.2660430223139146_wp, &
                                        0.2956193664498032_wp,0.3245763905312779_wp,0.3527670026206972_wp, &
                                        0.2823853479435550e-1_wp,0.5741296374713106e-1_wp,0.8724646633650199e-1_wp, &
                                        0.1175034422915616_wp,0.1479755652628428_wp,0.1784740659484352_wp,0.2088245700431244_wp, &
                                        0.2388628136570763_wp,0.2684308928769185_wp,0.2973740761960252_wp, &
                                        0.2916399920493977e-1_wp,0.5898803024755659e-1_wp,0.8924162698525409e-1_wp, &
                                        0.1197185199637321_wp,0.1502300756161382_wp,0.1806004191913564_wp,0.2106621764786252_wp, &
                                        0.2402526932671914_wp,0.2982529203607657e-1_wp,0.6008728062339922e-1_wp, &
                                        0.9058227674571398e-1_wp,0.1211219235803400_wp,0.1515286404791580_wp, &
                                        0.1816314681255552_wp,0.3026991752575440e-1_wp,0.6078402297870770e-1_wp, &
                                        0.9135459984176636e-1_wp,0.1218024155966590_wp,0.3052608357660639e-1_wp, &
                                        0.6112185773983089e-1_wp], &
                              v1 = 0.9080510764308163e-4_wp, v2 = Zero, v3 = 0.2084824361987793e-3_wp, &
                              v4(N1) = [0.5011105657239616e-4_wp,0.5942520409683854e-4_wp,0.9564394826109721e-4_wp, &
                                        0.1185530657126338e-3_wp,0.1364510114230331e-3_wp,0.1505828825605415e-3_wp, &
                                        0.1619298749867023e-3_wp,0.1712450504267789e-3_wp,0.1789891098164999e-3_wp, &
                                        0.1854474955629795e-3_wp,0.1908148636673661e-3_wp,0.1952377405281833e-3_wp, &
                                        0.1988349254282232e-3_wp,0.2017079807160050e-3_wp,0.2039473082709094e-3_wp, &
                                        0.2056360279288953e-3_wp,0.2068525823066865e-3_wp,0.2076724877534488e-3_wp, &
                                        0.2081694278237885e-3_wp,0.2084157631219326e-3_wp,0.2084381531128593e-3_wp, &
                                        0.2083476277129307e-3_wp,0.2082686194459732e-3_wp,0.2082475686112415e-3_wp, &
                                        0.2083139860289915e-3_wp,0.2084745561831237e-3_wp,0.2087091313375890e-3_wp, &
                                        0.2089718413297697e-3_wp,0.2092003303479793e-3_wp,0.2093336148263241e-3_wp], &
                              v5(N2) = [0.7591708117365267e-4_wp,0.1083383968169186e-3_wp,0.1403019395292510e-3_wp, &
                                        0.1615970179286436e-3_wp,0.1771144187504911e-3_wp,0.1887760022988168e-3_wp, &
                                        0.1973474670768214e-3_wp,0.2033787661234659e-3_wp,0.2072343626517331e-3_wp, &
                                        0.2091177834226918e-3_wp], &
                              v6(N3) = [0.9316684484675566e-4_wp,0.1116193688682976e-3_wp,0.1298623551559414e-3_wp, &
                                        0.1450236832456426e-3_wp,0.1572719958149914e-3_wp,0.1673234785867195e-3_wp, &
                                        0.1756860118725188e-3_wp,0.1826776290439367e-3_wp,0.1885116347992865e-3_wp, &
                                        0.1933457860170574e-3_wp,0.1973060671902064e-3_wp,0.2004987099616311e-3_wp, &
                                        0.2030170909281499e-3_wp,0.2049461460119080e-3_wp,0.2063653565200186e-3_wp, &
                                        0.2073507927381027e-3_wp,0.2079764593256122e-3_wp,0.2083150534968778e-3_wp, &
                                        0.1262715121590664e-3_wp,0.1414386128545972e-3_wp,0.1538740401313898e-3_wp, &
                                        0.1642434942331432e-3_wp,0.1729790609237496e-3_wp,0.1803505190260828e-3_wp, &
                                        0.1865475350079657e-3_wp,0.1917182669679069e-3_wp,0.1959851709034382e-3_wp, &
                                        0.1994529548117882e-3_wp,0.2022138911146548e-3_wp,0.2043518024208592e-3_wp, &
                                        0.2059450313018110e-3_wp,0.2070685715318472e-3_wp,0.2077955310694373e-3_wp, &
                                        0.2081980387824712e-3_wp,0.1521318610377956e-3_wp,0.1622772720185755e-3_wp, &
                                        0.1710498139420709e-3_wp,0.1785911149448736e-3_wp,0.1850125313687736e-3_wp, &
                                        0.1904229703933298e-3_wp,0.1949259956121987e-3_wp,0.1986161545363960e-3_wp, &
                                        0.2015790585641370e-3_wp,0.2038934198707418e-3_wp,0.2056334060538251e-3_wp, &
                                        0.2068705959462289e-3_wp,0.2076753906106002e-3_wp,0.2081179391734803e-3_wp, &
                                        0.1700345216228943e-3_wp,0.1774906779990410e-3_wp,0.1839659377002642e-3_wp, &
                                        0.1894987462975169e-3_wp,0.1941548809452595e-3_wp,0.1980078427252384e-3_wp, &
                                        0.2011296284744488e-3_wp,0.2035888456966776e-3_wp,0.2054516325352142e-3_wp, &
                                        0.2067831033092635e-3_wp,0.2076485320284876e-3_wp,0.2081141439525255e-3_wp, &
                                        0.1834383015469222e-3_wp,0.1889540591777677e-3_wp,0.1936677023597375e-3_wp, &
                                        0.1976176495066504e-3_wp,0.2008536004560983e-3_wp,0.2034280351712291e-3_wp, &
                                        0.2053944466027758e-3_wp,0.2068077642882360e-3_wp,0.2077250949661599e-3_wp, &
                                        0.2082062440705320e-3_wp,0.1934374486546626e-3_wp,0.1974107010484300e-3_wp, &
                                        0.2007129290388658e-3_wp,0.2033736947471293e-3_wp,0.2054287125902493e-3_wp, &
                                        0.2069184936818894e-3_wp,0.2078883689808782e-3_wp,0.2083886366116359e-3_wp, &
                                        0.2006593275470817e-3_wp,0.2033728426135397e-3_wp,0.2055008781377608e-3_wp, &
                                        0.2070651783518502e-3_wp,0.2080953335094320e-3_wp,0.2086284998988521e-3_wp, &
                                        0.2055549387644668e-3_wp,0.2071871850267654e-3_wp,0.2082856600431965e-3_wp, &
                                        0.2088705858819358e-3_wp,0.2083995867536322e-3_wp,0.2090509712889637e-3_wp]

  call ld_all(sym,order(62,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld5294

subroutine ld5438(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(63,sym)), y(order(63,sym)), z(order(63,sym)), w(order(63,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 30, N2 = 14, N3 = 91
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(63,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld5438

subroutine ld5606(sym,x,y,z,w)

  !*********************************************************************
  !
  ! This is a placeholder, the rule does not exist (yet)

  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(64,sym)), y(order(64,sym)), z(order(64,sym)), w(order(64,sym))
  integer(kind=iwp) :: i
  integer(kind=iwp), parameter :: N1 = 30, N2 = 15, N3 = 94
  real(kind=wp), parameter :: a4(N1) = [(Zero,i=1,N1)], &
                              a5(N2) = [(Zero,i=1,N2)], &
                              a6(N3) = [(Zero,i=1,N3)], &
                              b6(N3) = [(Zero,i=1,N3)], &
                              v1 = Zero, v2 = Zero, v3 = Zero, &
                              v4(N1) = [(Zero,i=1,N1)], &
                              v5(N2) = [(Zero,i=1,N2)], &
                              v6(N3) = [(Zero,i=1,N3)]

  write(u6,'(a)') ' '
  write(u6,'(a)') 'LDxxxx - Fatal error!'
  write(u6,'(a)') '  This rule is not implemented.'
  call Abend()

  call ld_all(sym,order(64,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld5606

subroutine ld5810(sym,x,y,z,w)

  !*********************************************************************
  !
  !! LD5810 computes the 5810 point Lebedev angular grid.
  !
  !  Modified:
  !
  !    09 September 2010
  !    26 February 2024
  !
  !  Author:
  !
  !    Dmitri Laikov
  !
  !  Reference:
  !
  !    Vyacheslav Lebedev, Dmitri Laikov,
  !    A quadrature formula for the sphere of the 131st
  !    algebraic order of accuracy,
  !    Russian Academy of Sciences Doklady Mathematics,
  !    Volume 59, Number 3, 1999, pages 477-481.
  !
  !  Parameters:
  !
  !    Output, real ( kind = wp ) X(N), Y(N), Z(N), W(N), the coordinates
  !    and weights of the points.

  integer(kind=iwp), intent(in) :: sym
  real(kind=wp), intent(out) :: x(order(65,sym)), y(order(65,sym)), z(order(65,sym)), w(order(65,sym))
  integer(kind=iwp), parameter :: N1 = 31, N2 = 10, N3 = 100
  real(kind=wp), parameter :: a4(N1) = [0.1182361662400277e-1_wp,0.3062145009138958e-1_wp,0.5329794036834243e-1_wp, &
                                        0.7848165532862220e-1_wp,0.1054038157636201_wp,0.1335577797766211_wp, &
                                        0.1625769955502252_wp,0.1921787193412792_wp,0.2221340534690548_wp,0.2522504912791132_wp, &
                                        0.2823610860679697_wp,0.3123173966267560_wp,0.3419847036953789_wp,0.3712386456999758_wp, &
                                        0.3999627649876828_wp,0.4280466458648093_wp,0.4553844360185711_wp,0.4818736094437834_wp, &
                                        0.5074138709260629_wp,0.5319061304570707_wp,0.5552514978677286_wp,0.5981009025246183_wp, &
                                        0.6173990192228116_wp,0.6351365239411131_wp,0.6512010228227200_wp,0.6654758363948120_wp, &
                                        0.6778410414853370_wp,0.6881760887484110_wp,0.6963645267094598_wp,0.7023010617153579_wp, &
                                        0.7059004636628753_wp], &
                              a5(N2) = [0.3552470312472575e-1_wp,0.9151176620841283e-1_wp,0.1566197930068980_wp, &
                                        0.2265467599271907_wp,0.2988242318581361_wp,0.3717482419703886_wp,0.4440094491758889_wp, &
                                        0.5145337096756642_wp,0.5824053672860230_wp,0.6468283961043370_wp], &
                              a6(N3) = [0.6095964259104373e-1_wp,0.8811962270959388e-1_wp,0.1165936722428831_wp, &
                                        0.1460232857031785_wp,0.1761197110181755_wp,0.2066471190463718_wp,0.2374076026328152_wp, &
                                        0.2682305474337051_wp,0.2989653312142369_wp,0.3294762752772209_wp,0.3596390887276086_wp, &
                                        0.3893383046398812_wp,0.4184653789358347_wp,0.4469172319076166_wp,0.4745950813276976_wp, &
                                        0.5014034601410262_wp,0.5272493404551239_wp,0.5520413051846366_wp,0.5756887237503077_wp, &
                                        0.1225039430588352_wp,0.1539113217321372_wp,0.1856213098637712_wp,0.2174998728035131_wp, &
                                        0.2494128336938330_wp,0.2812321562143480_wp,0.3128372276456111_wp,0.3441145160177973_wp, &
                                        0.3749567714853510_wp,0.4052621732015610_wp,0.4349335453522385_wp,0.4638776641524965_wp, &
                                        0.4920046410462687_wp,0.5192273554861704_wp,0.5454609081136522_wp,0.5706220661424140_wp, &
                                        0.5946286755181518_wp,0.1905370790924295_wp,0.2242518717748009_wp,0.2577190808025936_wp, &
                                        0.2908724534927187_wp,0.3236354020056219_wp,0.3559267359304543_wp,0.3876637123676956_wp, &
                                        0.4187636705218842_wp,0.4491449019883107_wp,0.4787270932425445_wp,0.5074315153055574_wp, &
                                        0.5351810507738336_wp,0.5619001025975381_wp,0.5875144035268046_wp,0.6119507308734495_wp, &
                                        0.2619733870119463_wp,0.2968149743237949_wp,0.3310451504860488_wp,0.3646215567376676_wp, &
                                        0.3974916785279360_wp,0.4295967403772029_wp,0.4608742854473447_wp,0.4912598858949903_wp, &
                                        0.5206882758945558_wp,0.5490940914019819_wp,0.5764123302025542_wp,0.6025786004213506_wp, &
                                        0.6275291964794956_wp,0.3348189479861771_wp,0.3699515545855295_wp,0.4042003071474669_wp, &
                                        0.4375320100182624_wp,0.4699054490335947_wp,0.5012739879431952_wp,0.5315874883754966_wp, &
                                        0.5607937109622117_wp,0.5888393223495521_wp,0.6156705979160163_wp,0.6412338809078123_wp, &
                                        0.4076051259257167_wp,0.4423788125791520_wp,0.4760480917328258_wp,0.5085838725946297_wp, &
                                        0.5399513637391218_wp,0.5701118433636380_wp,0.5990240530606021_wp,0.6266452685139695_wp, &
                                        0.6529320971415942_wp,0.4791583834610126_wp,0.5130373952796940_wp,0.5456252429628476_wp, &
                                        0.5768956329682385_wp,0.6068186944699046_wp,0.6353622248024907_wp,0.6624927035731797_wp, &
                                        0.5484933508028488_wp,0.5810207682142106_wp,0.6120955197181352_wp,0.6416944284294319_wp, &
                                        0.6697926391731260_wp,0.6147594390585488_wp,0.6455390026356783_wp,0.6747258588365477_wp, &
                                        0.6772135750395347_wp], &
                              b6(N3) = [0.1787828275342931e-1_wp,0.3953888740792096e-1_wp,0.6378121797722990e-1_wp, &
                                        0.8985890813745037e-1_wp,0.1172606510576162_wp,0.1456102876970995_wp, &
                                        0.1746153823011775_wp,0.2040383070295584_wp,0.2336788634003698_wp,0.2633632752654219_wp, &
                                        0.2929369098051601_wp,0.3222592785275512_wp,0.3512004791195743_wp,0.3796385677684537_wp, &
                                        0.4074575378263879_wp,0.4345456906027828_wp,0.4607942515205134_wp,0.4860961284181720_wp, &
                                        0.5103447395342790_wp,0.2136455922655793e-1_wp,0.4520926166137188e-1_wp, &
                                        0.7086468177864818e-1_wp,0.9785239488772918e-1_wp,0.1258106396267210_wp, &
                                        0.1544529125047001_wp,0.1835433512202753_wp,0.2128813258619585_wp,0.2422913734880829_wp, &
                                        0.2716163748391453_wp,0.3007127671240280_wp,0.3294470677216479_wp,0.3576932543699155_wp, &
                                        0.3853307059757764_wp,0.4122425044452694_wp,0.4383139587781027_wp,0.4634312536300553_wp, &
                                        0.2371311537781979e-1_wp,0.4917878059254806e-1_wp,0.7595498960495142e-1_wp, &
                                        0.1036991083191100_wp,0.1321348584450234_wp,0.1610316571314789_wp,0.1901912080395707_wp, &
                                        0.2194384950137950_wp,0.2486155334763858_wp,0.2775768931812335_wp,0.3061863786591120_wp, &
                                        0.3343144718152556_wp,0.3618362729028427_wp,0.3886297583620408_wp,0.4145742277792031_wp, &
                                        0.2540047186389353e-1_wp,0.5208107018543989e-1_wp,0.7971828470885599e-1_wp, &
                                        0.1080465999177927_wp,0.1368413849366629_wp,0.1659073184763559_wp,0.1950703730454614_wp, &
                                        0.2241721144376724_wp,0.2530655255406489_wp,0.2816118409731066_wp,0.3096780504593238_wp, &
                                        0.3371348366394987_wp,0.3638547827694396_wp,0.2664841935537443e-1_wp, &
                                        0.5424000066843495e-1_wp,0.8251992715430854e-1_wp,0.1112695182483710_wp, &
                                        0.1402964116467816_wp,0.1694275117584291_wp,0.1985038235312689_wp,0.2273765660020893_wp, &
                                        0.2559041492849764_wp,0.2839497251976899_wp,0.3113791060500690_wp, &
                                        0.2757792290858463e-1_wp,0.5584136834984293e-1_wp,0.8457772087727143e-1_wp, &
                                        0.1135975846359248_wp,0.1427286904765053_wp,0.1718112740057635_wp,0.2006944855985351_wp, &
                                        0.2292335090598907_wp,0.2572871512353714_wp,0.2826094197735932e-1_wp, &
                                        0.5699871359683649e-1_wp,0.8602712528554394e-1_wp,0.1151748137221281_wp, &
                                        0.1442811654136362_wp,0.1731930321657680_wp,0.2017619958756061_wp, &
                                        0.2874219755907391e-1_wp,0.5778312123713695e-1_wp,0.8695262371439526e-1_wp, &
                                        0.1160893767057166_wp,0.1450378826743251_wp,0.2904957622341456e-1_wp, &
                                        0.5823809152617197e-1_wp,0.8740384899884715e-1_wp,0.2919946135808105e-1_wp], &
                              v1 = 0.9735347946175486e-5_wp, v2 = 0.1907581241803167e-3_wp, v3 = 0.1901059546737578e-3_wp, &
                              v4(N1) = [0.3926424538919212e-4_wp,0.6667905467294382e-4_wp,0.8868891315019135e-4_wp, &
                                        0.1066306000958872e-3_wp,0.1214506743336128e-3_wp,0.1338054681640871e-3_wp, &
                                        0.1441677023628504e-3_wp,0.1528880200826557e-3_wp,0.1602330623773609e-3_wp, &
                                        0.1664102653445244e-3_wp,0.1715845854011323e-3_wp,0.1758901000133069e-3_wp, &
                                        0.1794382485256736e-3_wp,0.1823238106757407e-3_wp,0.1846293252959976e-3_wp, &
                                        0.1864284079323098e-3_wp,0.1877882694626914e-3_wp,0.1887716321852025e-3_wp, &
                                        0.1894381638175673e-3_wp,0.1898454899533629e-3_wp,0.1900497929577815e-3_wp, &
                                        0.1900671501924092e-3_wp,0.1899837555533510e-3_wp,0.1899014113156229e-3_wp, &
                                        0.1898581257705106e-3_wp,0.1898804756095753e-3_wp,0.1899793610426402e-3_wp, &
                                        0.1901464554844117e-3_wp,0.1903533246259542e-3_wp,0.1905556158463228e-3_wp, &
                                        0.1907037155663528e-3_wp], &
                              v5(N2) = [0.5992997844249967e-4_wp,0.9749059382456978e-4_wp,0.1241680804599158e-3_wp, &
                                        0.1437626154299360e-3_wp,0.1584200054793902e-3_wp,0.1694436550982744e-3_wp, &
                                        0.1776617014018108e-3_wp,0.1836132434440077e-3_wp,0.1876494727075983e-3_wp, &
                                        0.1899906535336482e-3_wp], &
                              v6(N3) = [0.8143252820767350e-4_wp,0.9998859890887728e-4_wp,0.1156199403068359e-3_wp, &
                                        0.1287632092635513e-3_wp,0.1398378643365139e-3_wp,0.1491876468417391e-3_wp, &
                                        0.1570855679175456e-3_wp,0.1637483948103775e-3_wp,0.1693500566632843e-3_wp, &
                                        0.1740322769393633e-3_wp,0.1779126637278296e-3_wp,0.1810908108835412e-3_wp, &
                                        0.1836529132600190e-3_wp,0.1856752841777379e-3_wp,0.1872270566606832e-3_wp, &
                                        0.1883722645591307e-3_wp,0.1891714324525297e-3_wp,0.1896827480450146e-3_wp, &
                                        0.1899628417059528e-3_wp,0.1123301829001669e-3_wp,0.1253698826711277e-3_wp, &
                                        0.1366266117678531e-3_wp,0.1462736856106918e-3_wp,0.1545076466685412e-3_wp, &
                                        0.1615096280814007e-3_wp,0.1674366639741759e-3_wp,0.1724225002437900e-3_wp, &
                                        0.1765810822987288e-3_wp,0.1800104126010751e-3_wp,0.1827960437331284e-3_wp, &
                                        0.1850140300716308e-3_wp,0.1867333507394938e-3_wp,0.1880178688638289e-3_wp, &
                                        0.1889278925654758e-3_wp,0.1895213832507346e-3_wp,0.1898548277397420e-3_wp, &
                                        0.1349105935937341e-3_wp,0.1444060068369326e-3_wp,0.1526797390930008e-3_wp, &
                                        0.1598208771406474e-3_wp,0.1659354368615331e-3_wp,0.1711279910946440e-3_wp, &
                                        0.1754952725601440e-3_wp,0.1791247850802529e-3_wp,0.1820954300877716e-3_wp, &
                                        0.1844788524548449e-3_wp,0.1863409481706220e-3_wp,0.1877433008795068e-3_wp, &
                                        0.1887444543705232e-3_wp,0.1894009829375006e-3_wp,0.1897683345035198e-3_wp, &
                                        0.1517327037467653e-3_wp,0.1587740557483543e-3_wp,0.1649093382274097e-3_wp, &
                                        0.1701915216193265e-3_wp,0.1746847753144065e-3_wp,0.1784555512007570e-3_wp, &
                                        0.1815687562112174e-3_wp,0.1840864370663302e-3_wp,0.1860676785390006e-3_wp, &
                                        0.1875690583743703e-3_wp,0.1886453236347225e-3_wp,0.1893501123329645e-3_wp, &
                                        0.1897366184519868e-3_wp,0.1643908815152736e-3_wp,0.1696300350907768e-3_wp, &
                                        0.1741553103844483e-3_wp,0.1780015282386092e-3_wp,0.1812116787077125e-3_wp, &
                                        0.1838323158085421e-3_wp,0.1859113119837737e-3_wp,0.1874969220221698e-3_wp, &
                                        0.1886375612681076e-3_wp,0.1893819575809276e-3_wp,0.1897794748256767e-3_wp, &
                                        0.1738963926584846e-3_wp,0.1777442359873466e-3_wp,0.1810010815068719e-3_wp, &
                                        0.1836920318248129e-3_wp,0.1858489473214328e-3_wp,0.1875079342496592e-3_wp, &
                                        0.1887080239102310e-3_wp,0.1894905752176822e-3_wp,0.1898991061200695e-3_wp, &
                                        0.1809065016458791e-3_wp,0.1836297121596799e-3_wp,0.1858426916241869e-3_wp, &
                                        0.1875654101134641e-3_wp,0.1888240751833503e-3_wp,0.1896497383866979e-3_wp, &
                                        0.1900775530219121e-3_wp,0.1858525041478814e-3_wp,0.1876248690077947e-3_wp, &
                                        0.1889404439064607e-3_wp,0.1898168539265290e-3_wp,0.1902779940661772e-3_wp, &
                                        0.1890125641731815e-3_wp,0.1899434637795751e-3_wp,0.1904520856831751e-3_wp, &
                                        0.1905534498734563e-3_wp]

  call ld_all(sym,order(65,sym),N1,N2,N3,a4,a5,a6,b6,v1,v2,v3,v4,v5,v6,x,y,z,w)

  return

end subroutine ld5810

end module Lebedev_quadrature
