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

subroutine dqk61(f,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk61
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  61-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real*8 version
!
!           parameters
!            on entry
!              f      - real*8
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real*8
!                       lower limit of integration
!
!              b      - real*8
!                       upper limit of integration
!
!            on return
!              reslt  - real*8
!                       approximation to the integral i
!                       result is computed by applying the 61-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 30-point gauss rule (resg).
!
!              abserr - real*8
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-reslt)
!
!              resabs - real*8
!                       approximation to the integral j
!
!              resasc - real*8
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk61

use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
interface
  function f(x)
    import :: wp
    real(kind=wp) :: f
    real(kind=wp), intent(in) :: x
  end function f
end interface
real(kind=wp), intent(in) :: a, b
real(kind=wp), intent(out) :: reslt, abserr, resabs, resasc
real(kind=wp) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(30), fv2(30), hlgth, resg, resk, reskh, uflow
integer(kind=iwp) :: j, jtw, jtwm1
real(kind=wp), external :: d1mach

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
! xgk    - abscissae of the 61-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 30-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 30-point gauss rule
!
! wgk    - weights of the 61-point kronrod rule
!
! wg     - weights of the 30-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

real(kind=wp), parameter :: wg(15) = [0.007968192496166605615465883474674_wp, &
                                      0.018466468311090959142302131912047_wp, &
                                      0.028784707883323369349719179611292_wp, &
                                      0.038799192569627049596801936446348_wp, &
                                      0.048402672830594052902938140422808_wp, &
                                      0.057493156217619066481721689402056_wp, &
                                      0.065974229882180495128128515115962_wp, &
                                      0.073755974737705206268243850022191_wp, &
                                      0.080755895229420215354694938460530_wp, &
                                      0.086899787201082979802387530715126_wp, &
                                      0.092122522237786128717632707087619_wp, &
                                      0.096368737174644259639468626351810_wp, &
                                      0.099593420586795267062780282103569_wp, &
                                      0.101762389748405504596428952168554_wp, &
                                      0.102852652893558840341285636705415_wp]
real(kind=wp), parameter :: xgk(31) = [0.999484410050490637571325895705811_wp, &
                                       0.996893484074649540271630050918695_wp, &
                                       0.991630996870404594858628366109486_wp, &
                                       0.983668123279747209970032581605663_wp, &
                                       0.973116322501126268374693868423707_wp, &
                                       0.960021864968307512216871025581798_wp, &
                                       0.944374444748559979415831324037439_wp, &
                                       0.926200047429274325879324277080474_wp, &
                                       0.905573307699907798546522558925958_wp, &
                                       0.882560535792052681543116462530226_wp, &
                                       0.857205233546061098958658510658944_wp, &
                                       0.829565762382768397442898119732502_wp, &
                                       0.799727835821839083013668942322683_wp, &
                                       0.767777432104826194917977340974503_wp, &
                                       0.733790062453226804726171131369528_wp, &
                                       0.697850494793315796932292388026640_wp, &
                                       0.660061064126626961370053668149271_wp, &
                                       0.620526182989242861140477556431189_wp, &
                                       0.579345235826361691756024932172540_wp, &
                                       0.536624148142019899264169793311073_wp, &
                                       0.492480467861778574993693061207709_wp, &
                                       0.447033769538089176780609900322854_wp, &
                                       0.400401254830394392535476211542661_wp, &
                                       0.352704725530878113471037207089374_wp, &
                                       0.304073202273625077372677107199257_wp, &
                                       0.254636926167889846439805129817805_wp, &
                                       0.204525116682309891438957671002025_wp, &
                                       0.153869913608583546963794672743256_wp, &
                                       0.102806937966737030147096751318001_wp, &
                                       0.051471842555317695833025213166723_wp, &
                                       0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(31) = [0.001389013698677007624551591226760_wp, &
                                       0.003890461127099884051267201844516_wp, &
                                       0.006630703915931292173319826369750_wp, &
                                       0.009273279659517763428441146892024_wp, &
                                       0.011823015253496341742232898853251_wp, &
                                       0.014369729507045804812451432443580_wp, &
                                       0.016920889189053272627572289420322_wp, &
                                       0.019414141193942381173408951050128_wp, &
                                       0.021828035821609192297167485738339_wp, &
                                       0.024191162078080601365686370725232_wp, &
                                       0.026509954882333101610601709335075_wp, &
                                       0.028754048765041292843978785354334_wp, &
                                       0.030907257562387762472884252943092_wp, &
                                       0.032981447057483726031814191016854_wp, &
                                       0.034979338028060024137499670731468_wp, &
                                       0.036882364651821229223911065617136_wp, &
                                       0.038678945624727592950348651532281_wp, &
                                       0.040374538951535959111995279752468_wp, &
                                       0.041969810215164246147147541285970_wp, &
                                       0.043452539701356069316831728117073_wp, &
                                       0.044814800133162663192355551616723_wp, &
                                       0.046059238271006988116271735559374_wp, &
                                       0.047185546569299153945261478181099_wp, &
                                       0.048185861757087129140779492298305_wp, &
                                       0.049055434555029778887528165367238_wp, &
                                       0.049795683427074206357811569379942_wp, &
                                       0.050405921402782346840893085653585_wp, &
                                       0.050881795898749606492297473049805_wp, &
                                       0.051221547849258772170656282604944_wp, &
                                       0.051426128537459025933862879215781_wp, &
                                       0.051494729429451567558340433647099_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc  - abscissa
! fval*  - function value
! resg   - result of the 30-point gauss formula
! resk   - result of the 61-point kronrod formula
! reskh  - approximation to the mean value of f over (a,b),
!          i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk61

epmach = d1mach(4)
uflow = d1mach(1)

centr = Half*(a+b)
hlgth = Half*(b-a)
dhlgth = abs(hlgth)

! compute the 61-point kronrod approximation to
! the integral, and estimate the absolute error.

fc = f(centr)
resg = Zero
resk = fc*wgk(31)
resabs = abs(resk)
do j=1,15
  jtw = j*2
  absc = hlgth*xgk(jtw)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1+fval2
  resg = resg+wg(j)*fsum
  resk = resk+wgk(jtw)*fsum
  resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
end do
do j=1,15
  jtwm1 = j*2-1
  absc = hlgth*xgk(jtwm1)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1+fval2
  resk = resk+wgk(jtwm1)*fsum
  resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
end do
reskh = resk*Half
resasc = wgk(31)*abs(fc-reskh)
do j=1,30
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*dhlgth
resasc = resasc*dhlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk61
