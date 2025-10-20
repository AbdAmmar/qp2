
! ---

!complex*16 function crint_1(n, rho)
!
!  implicit none
!  include 'constants.include.F'
!
!  integer,    intent(in) :: n
!  complex*16, intent(in) :: rho
!
!  integer                :: i, mmax
!  double precision       :: rho_mod
!  double precision       :: tmp
!  complex*16             :: rho_inv, rho_exp
!
!  complex*16             :: crint_smallz
!
!  rho_mod = zabs(rho)
!
!  if(rho_mod < 3.5d0) then
!
!    if(rho_mod .lt. 0.35d0) then
!
!      select case(n)
!      case(0)
!        crint_1 = (((((((((1.3122532963802805073d-08 * rho &
!                - 1.450385222315046877d-07) * rho &
!                + 1.458916900093370682d-06) * rho &
!                - 0.132275132275132275d-04) * rho &
!                + 0.106837606837606838d-03) * rho &
!                - 0.757575757575757576d-03) * rho &
!                + 0.462962962962962963d-02) * rho &
!                - 0.238095238095238095d-01) * rho &
!                + 0.10000000000000000000d0) * rho &
!                - 0.33333333333333333333d0) * rho &
!                + 1.0d0
!      case(1)
!        crint_1 = (((((((((1.198144314086343d-08 * rho &
!                - 1.312253296380281d-07) * rho &
!                + 1.305346700083542d-06) * rho &
!                - 1.167133520074696d-05) * rho &
!                + 9.259259259259259d-05) * rho &
!                - 6.410256410256410d-04) * rho &
!                + 3.787878787878788d-03) * rho &
!                - 1.851851851851852d-02) * rho &
!                + 7.142857142857142d-02) * rho &
!                - 2.000000000000000d-01) * rho &
!                + 3.333333333333333d-01
!      case(2)
!        crint_1 = (((((((((1.102292768959436d-08 * rho &
!                - 1.198144314086343d-07) * rho &
!                + 1.181027966742252d-06) * rho &
!                - 1.044277360066834d-05) * rho &
!                + 8.169934640522875d-05) * rho &
!                - 5.555555555555556d-04) * rho &
!                + 3.205128205128205d-03) * rho &
!                - 1.515151515151515d-02) * rho &
!                + 5.555555555555555d-02) * rho &
!                - 1.428571428571428d-01) * rho &
!                + 2.000000000000000d-01 
!      case(3)
!        crint_1 = (((((((((1.020641452740218d-08 * rho &
!                - 1.102292768959436d-07) * rho &
!                + 1.078329882677709d-06) * rho &
!                - 9.448223733938020d-06) * rho &
!                + 7.309941520467836d-05) * rho &
!                - 4.901960784313725d-04) * rho &
!                + 2.777777777777778d-03) * rho &
!                - 1.282051282051282d-02) * rho &
!                + 4.545454545454546d-02) * rho &
!                - 1.111111111111111d-01) * rho &
!                + 1.428571428571428d-01 
!      case default
!        tmp = dble(n + n + 1)
!        crint_1 = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
!                - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
!                + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
!                - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
!                + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
!                - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
!                + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
!                - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
!                + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
!                - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
!                + 1.0d0 / tmp
!      end select 
!
!    else
!
!      crint_1 = crint_smallz(n, rho)
!
!    endif
!
!  else
!
!    rho_exp = 0.5d0 * zexp(-rho)
!    rho_inv = (1.d0, 0.d0) / rho
!
!    call zboysfun00_1(rho, crint_1)
!    do i = 1, n
!      crint_1 = ((dble(i) - 0.5d0) * crint_1 - rho_exp) * rho_inv
!    enddo
!
!  endif
!
!  return
!end
!
!! ---
!
!subroutine crint_1_vec(n_max, rho, vals)
!
!  implicit none
!  include 'constants.include.F'
!
!  integer,    intent(in)  :: n_max
!  complex*16, intent(in)  :: rho
!  complex*16, intent(out) :: vals(0:n_max)
!
!  integer                :: n
!  double precision       :: rho_mod
!  double precision       :: tmp
!  complex*16             :: rho_inv, rho_exp
!
!  complex*16             :: crint_smallz
!
!  rho_mod = zabs(rho)
!
!  if(rho_mod < 3.5d0) then
!
!    if(rho_mod .lt. 0.35d0) then
!
!      vals(0) = (((((((((1.3122532963802805073d-08 * rho &
!              - 1.450385222315046877d-07) * rho &
!              + 1.458916900093370682d-06) * rho &
!              - 0.132275132275132275d-04) * rho &
!              + 0.106837606837606838d-03) * rho &
!              - 0.757575757575757576d-03) * rho &
!              + 0.462962962962962963d-02) * rho &
!              - 0.238095238095238095d-01) * rho &
!              + 0.10000000000000000000d0) * rho &
!              - 0.33333333333333333333d0) * rho &
!              + 1.0d0
!
!      if(n_max > 0) then
!
!        vals(1) = (((((((((1.198144314086343d-08 * rho &
!                - 1.312253296380281d-07) * rho &
!                + 1.305346700083542d-06) * rho &
!                - 1.167133520074696d-05) * rho &
!                + 9.259259259259259d-05) * rho &
!                - 6.410256410256410d-04) * rho &
!                + 3.787878787878788d-03) * rho &
!                - 1.851851851851852d-02) * rho &
!                + 7.142857142857142d-02) * rho &
!                - 2.000000000000000d-01) * rho &
!                + 3.333333333333333d-01
!
!        if(n_max > 1) then
!
!          vals(2) = (((((((((1.102292768959436d-08 * rho &
!                  - 1.198144314086343d-07) * rho &
!                  + 1.181027966742252d-06) * rho &
!                  - 1.044277360066834d-05) * rho &
!                  + 8.169934640522875d-05) * rho &
!                  - 5.555555555555556d-04) * rho &
!                  + 3.205128205128205d-03) * rho &
!                  - 1.515151515151515d-02) * rho &
!                  + 5.555555555555555d-02) * rho &
!                  - 1.428571428571428d-01) * rho &
!                  + 2.000000000000000d-01 
!
!          if(n_max > 2) then
!
!            vals(3) = (((((((((1.020641452740218d-08 * rho &
!                    - 1.102292768959436d-07) * rho &
!                    + 1.078329882677709d-06) * rho &
!                    - 9.448223733938020d-06) * rho &
!                    + 7.309941520467836d-05) * rho &
!                    - 4.901960784313725d-04) * rho &
!                    + 2.777777777777778d-03) * rho &
!                    - 1.282051282051282d-02) * rho &
!                    + 4.545454545454546d-02) * rho &
!                    - 1.111111111111111d-01) * rho &
!                    + 1.428571428571428d-01 
!
!            do n = 4, n_max
!              tmp = dble(n + n + 1)
!              vals(n) = (((((((((2.755731922398589d-07 * rho / (tmp + 20.d0) &
!                      - 2.755731922398589d-06 / (tmp + 18.d0)) * rho &
!                      + 2.480158730158730d-05 / (tmp + 16.d0)) * rho &
!                      - 1.984126984126984d-04 / (tmp + 14.d0)) * rho &
!                      + 1.388888888888889d-03 / (tmp + 12.d0)) * rho &
!                      - 8.333333333333333d-03 / (tmp + 10.d0)) * rho &
!                      + 4.166666666666666d-02 / (tmp +  8.d0)) * rho &
!                      - 1.666666666666667d-01 / (tmp +  6.d0)) * rho &
!                      + 5.000000000000000d-01 / (tmp +  4.d0)) * rho &
!                      - 1.000000000000000d+00 / (tmp +  2.d0)) * rho &
!                      + 1.0d0 / tmp
!            enddo
!
!          endif ! n_max > 2
!        endif ! n_max > 1
!      endif ! n_max > 0
!
!    else
!
!      call crint_smallz_vec(n_max, rho, vals)
!
!    endif
!
!  else
!
!    rho_exp = 0.5d0 * zexp(-rho)
!    rho_inv = (1.d0, 0.d0) / rho
!
!    call zboysfun00_1(rho, vals(0))
!    do n = 1, n_max
!      vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - rho_exp) * rho_inv
!    enddo
!
!  endif
!
!  return
!end

! ---

complex*16 function crint_smallz(n, rho)

  BEGIN_DOC
  ! Standard version of rint
  END_DOC

  implicit none
  integer,    intent(in)      :: n
  complex*16, intent(in)      :: rho

  integer,          parameter :: kmax = 40
  double precision, parameter :: eps = 1.d-10

  integer                     :: k
  double precision            :: delta_mod
  complex*16                  :: rho_k, ct, delta_k

  ct           = 0.5d0 * zexp(-rho) * gamma(dble(n) + 0.5d0)
  crint_smallz = ct / gamma(dble(n) + 1.5d0)

  rho_k = (1.d0, 0.d0)
  do k = 1, kmax

    rho_k        = rho_k * rho
    delta_k      = ct * rho_k / gamma(dble(n+k) + 1.5d0)
    crint_smallz = crint_smallz + delta_k

    delta_mod = dsqrt(real(delta_k)*real(delta_k) + aimag(delta_k)*aimag(delta_k))
    if(delta_mod .lt. eps) return
  enddo

  if(delta_mod > eps) then
    write(*,*) ' pb in crint_smallz !'
    write(*,*) ' n, rho = ', n, rho
    write(*,*) ' value = ', crint_smallz
    write(*,*) ' delta_mod = ', delta_mod
    !stop 1
  endif

end

! ---

complex*16 function crint_2(n, rho)

  implicit none

  integer,    intent(in) :: n
  complex*16, intent(in) :: rho

  logical                :: debug_crint
  double precision       :: tmp, abs_rho
  complex*16             :: val

  complex*16, external   :: crint_smallz

  abs_rho = abs(rho)

  if(abs_rho < 3.5d0) then

    if(abs_rho .lt. 0.35d0) then

      select case(n)
      case(0)

        crint_2 = ((((((((((((((-2.4668270102644571d-14 * rho &
                + 3.9554295164585257d-13) * rho &
                - 5.9477940136376354d-12) * rho &
                + 8.3507027951472397d-11) * rho &
                - 1.0892221037148573d-09) * rho &
                + 1.3122532963802806d-08) * rho &
                - 1.4503852223150468d-07) * rho &
                + 1.4589169000933706d-06) * rho &
                - 1.3227513227513228d-05) * rho &
                + 1.0683760683760684d-04) * rho &
                - 7.5757575757575758d-04) * rho &
                + 4.6296296296296294d-03) * rho &
                - 2.3809523809523808d-02) * rho &
                + 1.0000000000000001d-01) * rho &
                - 3.3333333333333331d-01) * rho &
                + 1.0000000000000000d+00

      case(1)

        crint_2 = ((((((((((((((-2.3173223429757021d-14 * rho &
                + 3.7002405153966856d-13) * rho &
                - 5.5376013230419360d-12) * rho &
                + 7.7321322177289257d-11) * rho &
                - 1.0020843354176688d-09) * rho &
                + 1.1981443140863431d-08) * rho &
                - 1.3122532963802806d-07) * rho &
                + 1.3053467000835422d-06) * rho &
                - 1.1671335200746965d-05) * rho &
                + 9.2592592592592588d-05) * rho &
                - 6.4102564102564103d-04) * rho &
                + 3.7878787878787880d-03) * rho &
                - 1.8518518518518517d-02) * rho &
                + 7.1428571428571425d-02) * rho &
                - 2.0000000000000001d-01) * rho &
                + 3.3333333333333331d-01

      case(2)

        crint_2 = ((((((((((((((-2.1849039233770904d-14 * rho &
                + 3.4759835144635530d-13) * rho &
                - 5.1803367215553597d-12) * rho &
                + 7.1988817199545171d-11) * rho &
                - 9.2785586612747108d-10) * rho &
                + 1.1022927689594356d-08) * rho &
                - 1.1981443140863431d-07) * rho &
                + 1.1810279667422524d-06) * rho &
                - 1.0442773600668338d-05) * rho &
                + 8.1699346405228753d-05) * rho &
                - 5.5555555555555556d-04) * rho &
                + 3.2051282051282050d-03) * rho &
                - 1.5151515151515152d-02) * rho &
                + 5.5555555555555552d-02) * rho &
                - 1.4285714285714285d-01) * rho &
                + 2.0000000000000001d-01

      case(3)

        crint_2 = ((((((((((((((-2.0668010085999505d-14 * rho &
                + 3.2773558850656355d-13) * rho &
                - 4.8663769202489742d-12) * rho &
                + 6.7344377380219679d-11) * rho &
                - 8.6386580639454200d-10) * rho &
                + 1.0206414527402181d-08) * rho &
                - 1.1022927689594356d-07) * rho &
                + 1.0783298826777088d-06) * rho &
                - 9.4482237339380196d-06) * rho &
                + 7.3099415204678364d-05) * rho &
                - 4.9019607843137254d-04) * rho &
                + 2.7777777777777779d-03) * rho &
                - 1.2820512820512820d-02) * rho &
                + 4.5454545454545456d-02) * rho &
                - 1.1111111111111110d-01) * rho &
                + 1.4285714285714285d-01

      case(4)

        crint_2 = ((((((((((((((-1.9608112132871323d-14 * rho &
                + 3.1002015128999257d-13) * rho &
                - 4.5882982390918896d-12) * rho &
                + 6.3262899963236659d-11) * rho &
                - 8.0813252856263604d-10) * rho &
                + 9.5025238703399630d-09) * rho &
                - 1.0206414527402182d-07) * rho &
                + 9.9206349206349202d-07) * rho &
                - 8.6266390614216706d-06) * rho &
                + 6.6137566137566142d-05) * rho &
                - 4.3859649122807018d-04) * rho &
                + 2.4509803921568627d-03) * rho &
                - 1.1111111111111112d-02) * rho &
                + 3.8461538461538464d-02) * rho &
                - 9.0909090909090912d-02) * rho &
                + 1.1111111111111110d-01

      case(5)

        crint_2 = ((((((((((((((-1.8651618858097113d-14 * rho &
                + 2.9412168199306989d-13) * rho &
                - 4.3402821180598959d-12) * rho &
                + 5.9647877108194569d-11) * rho &
                - 7.5915479955883996d-10) * rho &
                + 8.8894578141889977d-09) * rho &
                - 9.5025238703399620d-08) * rho &
                + 9.1857730746619638d-07) * rho &
                - 7.9365079365079362d-06) * rho &
                + 6.0386473429951689d-05) * rho &
                - 3.9682539682539683d-04) * rho &
                + 2.1929824561403508d-03) * rho &
                - 9.8039215686274508d-03) * rho &
                + 3.3333333333333333d-02) * rho &
                - 7.6923076923076927d-02) * rho &
                + 9.0909090909090912d-02

      case(6)

        crint_2 = ((((((((((((((-1.7784101701906551d-14 * rho &
                + 2.7977428287145669d-13) * rho &
                - 4.1177035479029777d-12) * rho &
                + 5.6423667534778645d-11) * rho &
                - 7.1577452529833478d-10) * rho &
                + 8.3507027951472401d-09) * rho &
                - 8.8894578141889964d-08) * rho &
                + 8.5522714833059663d-07) * rho &
                - 7.3486184597295711d-06) * rho &
                + 5.5555555555555558d-05) * rho &
                - 3.6231884057971015d-04) * rho &
                + 1.9841269841269840d-03) * rho &
                - 8.7719298245614030d-03) * rho &
                + 2.9411764705882353d-02) * rho &
                - 6.6666666666666666d-02) * rho &
                + 7.6923076923076927d-02

      case(7)

        crint_2 = ((((((((((((((-1.6993697181821813d-14 * rho &
                + 2.6676152552859824d-13) * rho &
                - 3.9168399602003940d-12) * rho &
                + 5.3530146122738715d-11) * rho &
                - 6.7708401041734372d-10) * rho &
                + 7.8735197782816838d-09) * rho &
                - 8.3507027951472401d-08) * rho &
                + 8.0005120327700973d-07) * rho &
                - 6.8418171866447731d-06) * rho &
                + 5.1440329218106995d-05) * rho &
                - 3.3333333333333332d-04) * rho &
                + 1.8115942028985507d-03) * rho &
                - 7.9365079365079361d-03) * rho &
                + 2.6315789473684209d-02) * rho &
                - 5.8823529411764705d-02) * rho &
                + 6.6666666666666666d-02

      case(8)

        crint_2 = ((((((((((((((-1.6270561131531524d-14 * rho &
                + 2.5490545772732723d-13) * rho &
                - 3.7346613574003751d-12) * rho &
                + 5.0918919482605119d-11) * rho &
                - 6.4236175347286453d-10) * rho &
                + 7.4479241145907816d-09) * rho &
                - 7.8735197782816835d-08) * rho &
                + 7.5156325156325153d-07) * rho &
                - 6.4004096262160778d-06) * rho &
                + 4.7892720306513412d-05) * rho &
                - 3.0864197530864197d-04) * rho &
                + 1.6666666666666668d-03) * rho &
                - 7.2463768115942030d-03) * rho &
                + 2.3809523809523808d-02) * rho &
                - 5.2631578947368418d-02) * rho &
                + 5.8823529411764705d-02

      case(9)

        crint_2 = ((((((((((((((-1.5606456595550644d-14 * rho &
                + 2.4405841697297287d-13) * rho &
                - 3.5686764081825809d-12) * rho &
                + 4.8550597646204880d-11) * rho &
                - 6.1102703379126148d-10) * rho &
                + 7.0659792882015104d-09) * rho &
                - 7.4479241145907807d-08) * rho &
                + 7.0861678004535149d-07) * rho &
                - 6.0125060125060122d-06) * rho &
                + 4.4802867383512545d-05) * rho &
                - 2.8735632183908046d-04) * rho &
                + 1.5432098765432098d-03) * rho &
                - 6.6666666666666671d-03) * rho &
                + 2.1739130434782608d-02) * rho &
                - 4.7619047619047616d-02) * rho &
                + 5.2631578947368418d-02

      case default

        tmp = dble(n + n + 1)
        crint_2 = (((((((((((((((((-2.8114572543455206d-15 / (tmp + 34.d0)) * rho &
                + 4.7794773323873853d-14 / (tmp + 32.d0)) * rho &
                - 7.6471637318198164d-13 / (tmp + 30.d0)) * rho &
                + 1.1470745597729725d-11 / (tmp + 28.d0)) * rho &
                - 1.6059043836821613d-10 / (tmp + 26.d0)) * rho &
                + 2.0876756987868100d-09 / (tmp + 24.d0)) * rho &
                - 2.5052108385441720d-08 / (tmp + 22.d0)) * rho &
                + 2.7557319223985888d-07 / (tmp + 20.d0)) * rho &
                - 2.7557319223985893d-06 / (tmp + 18.d0)) * rho &
                + 2.4801587301587302d-05 / (tmp + 16.d0)) * rho &
                - 1.9841269841269841d-04 / (tmp + 14.d0)) * rho &
                + 1.3888888888888889d-03 / (tmp + 12.d0)) * rho &
                - 8.3333333333333332d-03 / (tmp + 10.d0)) * rho &
                + 4.1666666666666664d-02 / (tmp +  8.d0)) * rho &
                - 1.6666666666666666d-01 / (tmp +  6.d0)) * rho &
                + 5.0000000000000000d-01 / (tmp +  4.d0)) * rho &
                - 1.0000000000000000d+00 / (tmp +  2.d0)) * rho &
                + 1.0000000000000000d+00 / tmp

      end select 

    else

      crint_2 = crint_smallz(n, rho)

    endif

  else

    if(real(rho) .ge. 0.d0) then

      call zboysfun(n, rho, val)
      crint_2 = val

    else

      call zboysfunnrp(n, rho, val)
      crint_2 = val * zexp(-rho)

    endif

  endif

  debug_crint = .true.
  if(debug_crint) then
    write(6666,*) n, real(rho), aimag(rho), real(crint_2), aimag(crint_2)
  endif

  return
end

! ---

subroutine zboysfun_vec(n_max, x, vals)

  BEGIN_DOC
  !
  ! Computes values of the Boys function for n = 0, 1, ..., n_max
  ! for a complex valued argument
  !
  ! Input: x --- argument, complex*16, Re(x) >= 0
  ! Output: vals  --- values of the Boys function, n = 0, 1, ..., n_max
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  complex*16              :: yy, x_inv

  call zboysfun00_2(x, vals(0))

  yy = 0.5d0 * zexp(-x)
  x_inv = (1.d0, 0.d0) / x
  do n = 1, n_max
    vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - yy) * x_inv
  enddo

  return
end

! ---

subroutine zboysfun(n, x, val)

  BEGIN_DOC
  !
  ! Computes values of the Boys function for n
  ! for a complex valued argument
  !
  ! Input: x --- argument, complex*16, Re(x) >= 0
  ! Output: val  --- value of the Boys function
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: val

  integer                 :: i
  complex*16              :: yy, x_inv

  call zboysfun00_2(x, val)

  yy = 0.5d0 * zexp(-x)
  x_inv = (1.d0, 0.d0) / x
  do i = 1, n
    val = ((dble(i) - 0.5d0) * val - yy) * x_inv
  enddo

  return
end

! ---

subroutine zboysfunnrp_vec(n_max, x, vals)

  BEGIN_DOC
  !
  ! Computes values of e^z F(n,z) for n = 0, 1, ..., n_max
  ! (where F(n,z) are the Boys functions)
  ! for a complex valued argument WITH NEGATIVE REAL PART
  !
  ! Input: x  --- argument, complex *16 Re(x)<=0
  ! Output: vals  --- values of e^z F(n,z), n = 0, 1, ..., n_max
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  complex*16              :: x_inv

  call zboysfun00nrp(x, vals(0))

  x_inv = (1.d0, 0.d0) / x
  do n = 1, n_max
    vals(n) = ((dble(n) - 0.5d0) * vals(n-1) - 0.5d0) * x_inv
  enddo

  return
end

! ---

subroutine zboysfunnrp(n, x, val)

  BEGIN_DOC
  !
  ! Computes values of e^z F(n,z) for n
  ! (where F(n,z) are the Boys functions)
  ! for a complex valued argument WITH NEGATIVE REAL PART
  !
  ! Input: x  --- argument, complex *16 Re(x)<=0
  ! Output: val  --- value of e^z F(n,z)
  !
  ! Beylkin & Sharma, J. Chem. Phys. 155, 174117 (2021)
  ! https://doi.org/10.1063/5.0062444
  !
  END_DOC

  implicit none

  integer,    intent(in)  :: n
  complex*16, intent(in)  :: x
  complex*16, intent(out) :: val

  integer                 :: i
  complex*16              :: x_inv

  call zboysfun00nrp(x, val)

  x_inv = (1.d0, 0.d0) / x
  do i = 1, n
    val = ((dble(i) - 0.5d0) * val - 0.5d0) * x_inv
  enddo

  return
end

! ---

complex*16 function crint_sum(n_pt_out, rho, d1)

  implicit none

  integer,    intent(in)  :: n_pt_out
  complex*16, intent(in)  :: rho, d1(0:n_pt_out)
                          
  integer                 :: i
  integer                 :: n_max

  logical                 :: debug_crint

  complex*16, allocatable :: vals(:)

  n_max = shiftr(n_pt_out, 1)
  allocate(vals(0:n_max))

  call crint_2_vec(n_max, rho, vals)

  debug_crint = .true.
  if(debug_crint) then
    do i = 0, n_pt_out, 2
      write(7777,*) shiftr(i, 1), real(rho), aimag(rho), real(vals(shiftr(i, 1))), aimag(vals(shiftr(i, 1)))
    enddo
  endif

  crint_sum = d1(0) * vals(0)
  do i = 2, n_pt_out, 2
    crint_sum += d1(i) * vals(shiftr(i, 1))
  enddo

  deallocate(vals)

  return
end

! ---

subroutine crint_2_vec(n_max, rho, vals)

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n
  double precision        :: tmp, abs_rho
  complex*16              :: erho


  abs_rho = abs(rho)

  if(abs_rho < 3.5d0) then

    if(abs_rho .lt. 0.35d0) then

      vals(0) = ((((((((((((((-2.4668270102644571d-14 * rho &
              + 3.9554295164585257d-13) * rho &
              - 5.9477940136376354d-12) * rho &
              + 8.3507027951472397d-11) * rho &
              - 1.0892221037148573d-09) * rho &
              + 1.3122532963802806d-08) * rho &
              - 1.4503852223150468d-07) * rho &
              + 1.4589169000933706d-06) * rho &
              - 1.3227513227513228d-05) * rho &
              + 1.0683760683760684d-04) * rho &
              - 7.5757575757575758d-04) * rho &
              + 4.6296296296296294d-03) * rho &
              - 2.3809523809523808d-02) * rho &
              + 1.0000000000000001d-01) * rho &
              - 3.3333333333333331d-01) * rho &
              + 1.0000000000000000d+00

      if(n_max > 0) then

        vals(1) = ((((((((((((((-2.3173223429757021d-14 * rho &
                + 3.7002405153966856d-13) * rho &
                - 5.5376013230419360d-12) * rho &
                + 7.7321322177289257d-11) * rho &
                - 1.0020843354176688d-09) * rho &
                + 1.1981443140863431d-08) * rho &
                - 1.3122532963802806d-07) * rho &
                + 1.3053467000835422d-06) * rho &
                - 1.1671335200746965d-05) * rho &
                + 9.2592592592592588d-05) * rho &
                - 6.4102564102564103d-04) * rho &
                + 3.7878787878787880d-03) * rho &
                - 1.8518518518518517d-02) * rho &
                + 7.1428571428571425d-02) * rho &
                - 2.0000000000000001d-01) * rho &
                + 3.3333333333333331d-01

        if(n_max > 1) then

          vals(2) = ((((((((((((((-2.1849039233770904d-14 * rho &
                  + 3.4759835144635530d-13) * rho &
                  - 5.1803367215553597d-12) * rho &
                  + 7.1988817199545171d-11) * rho &
                  - 9.2785586612747108d-10) * rho &
                  + 1.1022927689594356d-08) * rho &
                  - 1.1981443140863431d-07) * rho &
                  + 1.1810279667422524d-06) * rho &
                  - 1.0442773600668338d-05) * rho &
                  + 8.1699346405228753d-05) * rho &
                  - 5.5555555555555556d-04) * rho &
                  + 3.2051282051282050d-03) * rho &
                  - 1.5151515151515152d-02) * rho &
                  + 5.5555555555555552d-02) * rho &
                  - 1.4285714285714285d-01) * rho &
                  + 2.0000000000000001d-01

          if(n_max > 2) then

            vals(3) = ((((((((((((((-2.0668010085999505d-14 * rho &
                    + 3.2773558850656355d-13) * rho &
                    - 4.8663769202489742d-12) * rho &
                    + 6.7344377380219679d-11) * rho &
                    - 8.6386580639454200d-10) * rho &
                    + 1.0206414527402181d-08) * rho &
                    - 1.1022927689594356d-07) * rho &
                    + 1.0783298826777088d-06) * rho &
                    - 9.4482237339380196d-06) * rho &
                    + 7.3099415204678364d-05) * rho &
                    - 4.9019607843137254d-04) * rho &
                    + 2.7777777777777779d-03) * rho &
                    - 1.2820512820512820d-02) * rho &
                    + 4.5454545454545456d-02) * rho &
                    - 1.1111111111111110d-01) * rho &
                    + 1.4285714285714285d-01

            if(n_max > 3) then

              vals(4) = ((((((((((((((-1.9608112132871323d-14 * rho &
                      + 3.1002015128999257d-13) * rho &
                      - 4.5882982390918896d-12) * rho &
                      + 6.3262899963236659d-11) * rho &
                      - 8.0813252856263604d-10) * rho &
                      + 9.5025238703399630d-09) * rho &
                      - 1.0206414527402182d-07) * rho &
                      + 9.9206349206349202d-07) * rho &
                      - 8.6266390614216706d-06) * rho &
                      + 6.6137566137566142d-05) * rho &
                      - 4.3859649122807018d-04) * rho &
                      + 2.4509803921568627d-03) * rho &
                      - 1.1111111111111112d-02) * rho &
                      + 3.8461538461538464d-02) * rho &
                      - 9.0909090909090912d-02) * rho &
                      + 1.1111111111111110d-01

              if(n_max > 4) then

                vals(5) = ((((((((((((((-1.8651618858097113d-14 * rho &
                        + 2.9412168199306989d-13) * rho &
                        - 4.3402821180598959d-12) * rho &
                        + 5.9647877108194569d-11) * rho &
                        - 7.5915479955883996d-10) * rho &
                        + 8.8894578141889977d-09) * rho &
                        - 9.5025238703399620d-08) * rho &
                        + 9.1857730746619638d-07) * rho &
                        - 7.9365079365079362d-06) * rho &
                        + 6.0386473429951689d-05) * rho &
                        - 3.9682539682539683d-04) * rho &
                        + 2.1929824561403508d-03) * rho &
                        - 9.8039215686274508d-03) * rho &
                        + 3.3333333333333333d-02) * rho &
                        - 7.6923076923076927d-02) * rho &
                        + 9.0909090909090912d-02

                if(n_max > 5) then

                  vals(6) = ((((((((((((((-1.7784101701906551d-14 * rho &
                          + 2.7977428287145669d-13) * rho &
                          - 4.1177035479029777d-12) * rho &
                          + 5.6423667534778645d-11) * rho &
                          - 7.1577452529833478d-10) * rho &
                          + 8.3507027951472401d-09) * rho &
                          - 8.8894578141889964d-08) * rho &
                          + 8.5522714833059663d-07) * rho &
                          - 7.3486184597295711d-06) * rho &
                          + 5.5555555555555558d-05) * rho &
                          - 3.6231884057971015d-04) * rho &
                          + 1.9841269841269840d-03) * rho &
                          - 8.7719298245614030d-03) * rho &
                          + 2.9411764705882353d-02) * rho &
                          - 6.6666666666666666d-02) * rho &
                          + 7.6923076923076927d-02

                  if(n_max > 6) then

                    vals(7) = ((((((((((((((-1.6993697181821813d-14 * rho &
                            + 2.6676152552859824d-13) * rho &
                            - 3.9168399602003940d-12) * rho &
                            + 5.3530146122738715d-11) * rho &
                            - 6.7708401041734372d-10) * rho &
                            + 7.8735197782816838d-09) * rho &
                            - 8.3507027951472401d-08) * rho &
                            + 8.0005120327700973d-07) * rho &
                            - 6.8418171866447731d-06) * rho &
                            + 5.1440329218106995d-05) * rho &
                            - 3.3333333333333332d-04) * rho &
                            + 1.8115942028985507d-03) * rho &
                            - 7.9365079365079361d-03) * rho &
                            + 2.6315789473684209d-02) * rho &
                            - 5.8823529411764705d-02) * rho &
                            + 6.6666666666666666d-02

                    if(n_max > 7) then

                      vals(8) = ((((((((((((((-1.6270561131531524d-14 * rho &
                              + 2.5490545772732723d-13) * rho &
                              - 3.7346613574003751d-12) * rho &
                              + 5.0918919482605119d-11) * rho &
                              - 6.4236175347286453d-10) * rho &
                              + 7.4479241145907816d-09) * rho &
                              - 7.8735197782816835d-08) * rho &
                              + 7.5156325156325153d-07) * rho &
                              - 6.4004096262160778d-06) * rho &
                              + 4.7892720306513412d-05) * rho &
                              - 3.0864197530864197d-04) * rho &
                              + 1.6666666666666668d-03) * rho &
                              - 7.2463768115942030d-03) * rho &
                              + 2.3809523809523808d-02) * rho &
                              - 5.2631578947368418d-02) * rho &
                              + 5.8823529411764705d-02

                      if(n_max > 8) then

                        vals(9) = ((((((((((((((-1.5606456595550644d-14 * rho &
                                + 2.4405841697297287d-13) * rho &
                                - 3.5686764081825809d-12) * rho &
                                + 4.8550597646204880d-11) * rho &
                                - 6.1102703379126148d-10) * rho &
                                + 7.0659792882015104d-09) * rho &
                                - 7.4479241145907807d-08) * rho &
                                + 7.0861678004535149d-07) * rho &
                                - 6.0125060125060122d-06) * rho &
                                + 4.4802867383512545d-05) * rho &
                                - 2.8735632183908046d-04) * rho &
                                + 1.5432098765432098d-03) * rho &
                                - 6.6666666666666671d-03) * rho &
                                + 2.1739130434782608d-02) * rho &
                                - 4.7619047619047616d-02) * rho &
                                + 5.2631578947368418d-02

                        do n = 10, n_max
                          tmp = dble(n + n + 1)
                          vals(n) = (((((((((((((((((-2.8114572543455206d-15 / (tmp + 34.d0)) * rho &
                                  + 4.7794773323873853d-14 / (tmp + 32.d0)) * rho &
                                  - 7.6471637318198164d-13 / (tmp + 30.d0)) * rho &
                                  + 1.1470745597729725d-11 / (tmp + 28.d0)) * rho &
                                  - 1.6059043836821613d-10 / (tmp + 26.d0)) * rho &
                                  + 2.0876756987868100d-09 / (tmp + 24.d0)) * rho &
                                  - 2.5052108385441720d-08 / (tmp + 22.d0)) * rho &
                                  + 2.7557319223985888d-07 / (tmp + 20.d0)) * rho &
                                  - 2.7557319223985893d-06 / (tmp + 18.d0)) * rho &
                                  + 2.4801587301587302d-05 / (tmp + 16.d0)) * rho &
                                  - 1.9841269841269841d-04 / (tmp + 14.d0)) * rho &
                                  + 1.3888888888888889d-03 / (tmp + 12.d0)) * rho &
                                  - 8.3333333333333332d-03 / (tmp + 10.d0)) * rho &
                                  + 4.1666666666666664d-02 / (tmp +  8.d0)) * rho &
                                  - 1.6666666666666666d-01 / (tmp +  6.d0)) * rho &
                                  + 5.0000000000000000d-01 / (tmp +  4.d0)) * rho &
                                  - 1.0000000000000000d+00 / (tmp +  2.d0)) * rho &
                                  + 1.0000000000000000d+00 / tmp
                        enddo

                      endif ! n_max > 8
                    endif ! n_max > 7
                  endif ! n_max > 6
                endif ! n_max > 5
              endif ! n_max > 4
            endif ! n_max > 3
          endif ! n_max > 2
        endif ! n_max > 1
      endif ! n_max > 0

    else

      call crint_smallz_vec(n_max, rho, vals)

    endif

  else

    if(real(rho) .ge. 0.d0) then

      call zboysfun_vec(n_max, rho, vals)

    else

      call zboysfunnrp_vec(n_max, rho, vals)
      erho = zexp(-rho)
      do n = 0, n_max
        vals(n) = vals(n) * erho
      enddo

    endif

  endif

  return
end

! ---

subroutine crint_smallz_vec(n_max, rho, vals)

  BEGIN_DOC
  ! Standard version of rint
  END_DOC

  implicit none
  integer,    intent(in)      :: n_max
  complex*16, intent(in)      :: rho
  complex*16, intent(out)     :: vals(0:n_max)

  integer,          parameter :: kmax = 40
  double precision, parameter :: eps = 1.d-10

  integer                     :: k, n
  complex*16                  :: ct, delta_k
  complex*16                  :: rhoe
  complex*16, allocatable     :: rho_k(:)


  allocate(rho_k(0:kmax))

  rho_k(0) = (1.d0, 0.d0)
  do k = 1, kmax
    rho_k(k) = rho_k(k-1) * rho
  enddo

  rhoe = 0.5d0 * zexp(-rho)

  do n = 0, n_max

    ct = rhoe * gamma(dble(n) + 0.5d0)
    vals(n) = ct / gamma(dble(n) + 1.5d0)
  
    do k = 1, kmax
      delta_k = ct * rho_k(k) / gamma(dble(n+k) + 1.5d0)
      vals(n) += delta_k
      if(abs(delta_k) .lt. eps) then
        exit
      endif
    enddo
  
    if(abs(delta_k) > eps) then
      write(*,*) ' pb in crint_smallz_vec !'
      write(*,*) ' n, rho = ', n, rho
      write(*,*) ' value = ', vals(n)
      write(*,*) ' |delta_k| = ', abs(delta_k)
    endif
  enddo

  deallocate(rho_k)

  return
end

! ---

subroutine crint_quad_1(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1)
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1 * tmp1 * tmp1
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp2 * tmp2
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1) * tmp1**n
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv
  end select

end

! ---

subroutine crint_quad_2(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  rhoc = zexp(-rho*tmp_inv)
  rhoe = (1.d0, 0.d0)
  select case (n)

    case (0)
      !do i_quad = 1, n_quad - 1
      !  tmp0 = tmp0 + tmp_inv
      !  rhoe = rhoe * rhoc
      !  tmp1 = (rhoe - 1.d0) / dsqrt(tmp0)
      !  crint_quad = crint_quad + coef(iand(i_quad, 3)) * tmp1
      !enddo
      !crint_quad = 1.d0 + crint_quad * 0.022222222222222223d0 * tmp_inv
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe / dsqrt(tmp0)
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp2
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0**n / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

  end select

end

! ---

subroutine crint_quad_12(n, rho, n_quad, crint_quad)

  implicit none

  integer,    intent(in)  :: n, n_quad
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: crint_quad

  integer                 :: i_quad
  double precision        :: tmp_inv, tmp0, tmp1, tmp2
  double precision        :: coef(0:3) = (/14.d0, 32.d0, 12.d0, 32.d0 /)
  complex*16              :: rhoc, rhoe

  tmp_inv = 1.d0 / dble(n_quad)

  crint_quad = 7.d0 * zexp(-rho)

  tmp0 = 0.d0
  rhoc = zexp(-rho*tmp_inv)
  rhoe = (1.d0, 0.d0)
  select case (n)

    case (0)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * zexp(-rho*tmp1)
      enddo
      crint_quad = crint_quad * 0.044444444444444446d0 * tmp_inv

    case (1)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (2)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (3)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0 * tmp0 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case (4)
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0 * tmp0
        tmp2 = tmp1 * tmp1 / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp2
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

    case default
      do i_quad = 1, n_quad - 1
        tmp0 = tmp0 + tmp_inv
        tmp1 = tmp0**n / dsqrt(tmp0)
        rhoe = rhoe * rhoc
        crint_quad = crint_quad + coef(iand(i_quad, 3)) * rhoe * tmp1
      enddo
      crint_quad = crint_quad * 0.022222222222222223d0 * tmp_inv

  end select

end

! ---

subroutine crint_quad_12_vec(n_max, rho, vals)

  implicit none

  integer,    intent(in)  :: n_max
  complex*16, intent(in)  :: rho
  complex*16, intent(out) :: vals(0:n_max)

  integer                 :: n

  do n = 0, n_max
    call crint_quad_12(n, rho, 10000000, vals(n))
  enddo

  return
end

! ---

