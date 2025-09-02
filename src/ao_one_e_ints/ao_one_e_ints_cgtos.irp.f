! ---

BEGIN_PROVIDER [double precision, ao_overlap_cgtos, (ao_num, ao_num)]

  implicit none

  integer          :: i, j, m, n, l, ii, jj, dim_1, power_A(3), power_B(3)
  double precision :: c, overlap
  double precision :: phiA, phiB
  complex*16       :: alpha, alpha_inv, A_center(3)
  complex*16       :: beta, beta_inv, B_center(3)
  complex*16       :: C1, C2
  complex*16       :: overlap1, overlap_x1, overlap_y1, overlap_z1
  complex*16       :: overlap2, overlap_x2, overlap_y2, overlap_z2

  ao_overlap_cgtos   = 0.d0

  dim_1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                      &
 !$OMP DEFAULT(NONE)                                                     &
 !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, C1, C2,                         &
 !$OMP         alpha, alpha_inv, A_center, power_A, phiA,                &
 !$OMP         beta, beta_inv, B_center, power_B, phiB, overlap,         &
 !$OMP         overlap_x1, overlap_y1, overlap_z1, overlap1,             &
 !$OMP         overlap_x2, overlap_y2, overlap_z2, overlap2)             &
 !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim_1, &
 !$OMP        ao_coef_cgtos_norm_ord_transp, ao_expo_cgtos_ord_transp,   &
 !$OMP        ao_overlap_cgtos, ao_expo_phase_ord_transp)

  do j = 1, ao_num

    jj = ao_nucl(j)
    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)

    do i = 1, ao_num

      ii = ao_nucl(i)
      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      do n = 1, ao_prim_num(j)

        alpha = ao_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        phiA = ao_expo_phase_ord_transp(4,n,j)
        do m = 1, 3
          A_center(m) = dcmplx(nucl_coord(jj,m))
        enddo

        do l = 1, ao_prim_num(i)

          beta = ao_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          phiB = ao_expo_phase_ord_transp(4,l,i)
          do m = 1, 3
            B_center(m) = dcmplx(nucl_coord(ii,m))
          enddo

          c = ao_coef_cgtos_norm_ord_transp(n,j) * ao_coef_cgtos_norm_ord_transp(l,i)

          C1 = zexp((0.d0, 1.d0) * (-phiA - phiB))
          C2 = zexp((0.d0, 1.d0) * ( phiA - phiB))

          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_x1, overlap_y1, overlap_z1, overlap1, dim_1)

          call overlap_cgaussian_xyz(conjg(A_center), B_center, conjg(alpha), beta, power_A, power_B, &
                                     conjg(A_center), B_center, overlap_x2, overlap_y2, overlap_z2, overlap2, dim_1)

          overlap = 0.5d0 * real(C1 * overlap1 + C2 * overlap2)

          ao_overlap_cgtos(i,j) = ao_overlap_cgtos(i,j) + c * overlap

          if(isnan(ao_overlap_cgtos(i,j))) then
            print*,'i, j', i, j
            print*,'l, n', l, n
            print*,'c, overlap', c, overlap
            stop
          endif
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_deriv2_cgtos_x, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_deriv2_cgtos_y, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_deriv2_cgtos_z, (ao_num, ao_num)]

  implicit none
  integer          :: i, j, m, n, l, ii, jj, dim_1, power_A(3), power_B(3)
  double precision :: c, deriv_tmp
  double precision :: phiA, phiB, aa
  complex*16       :: alpha, alpha_inv, A_center(3), C1
  complex*16       :: beta, beta_inv, B_center(3), C2
  complex*16       :: overlap_x, overlap_y, overlap_z, overlap
  complex*16       :: overlap_x0_1, overlap_y0_1, overlap_z0_1
  complex*16       :: overlap_x0_2, overlap_y0_2, overlap_z0_2 
  complex*16       :: overlap_m2_1, overlap_p2_1
  complex*16       :: overlap_m2_2, overlap_p2_2
  complex*16       :: deriv_tmp_1, deriv_tmp_2


  dim_1 = 100

  ! -- Dummy call to provide everything

  A_center(:) = (0.0d0, 0.d0)
  B_center(:) = (1.0d0, 0.d0)
  alpha       = (1.0d0, 0.d0)
  beta        = (0.1d0, 0.d0)
  power_A     = 1
  power_B     = 0
  call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                             A_center, B_center, overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim_1)

  ! ---

  !$OMP PARALLEL                                                          &
  !$OMP DEFAULT(NONE)                                                     &
  !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, aa, C1, C2,                     &
  !$OMP         A_center, power_A, alpha, alpha_inv, phiA,                &
  !$OMP         B_center, power_B, beta, beta_inv, phiB,                  &
  !$OMP         deriv_tmp, deriv_tmp_1, deriv_tmp_2,                      &
  !$OMP         overlap_x, overlap_y, overlap_z, overlap,                 &
  !$OMP         overlap_m2_1, overlap_p2_1, overlap_m2_2, overlap_p2_2,   &
  !$OMP         overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap_x0_2,   &
  !$OMP         overlap_y0_2, overlap_z0_2)                               &
  !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim_1, &
  !$OMP        ao_coef_cgtos_norm_ord_transp, ao_expo_cgtos_ord_transp,   & 
  !$OMP        ao_expo_phase_ord_transp,                                  &
  !$OMP        ao_deriv2_cgtos_x, ao_deriv2_cgtos_y, ao_deriv2_cgtos_z) 
  !$OMP DO SCHEDULE(GUIDED)
  do j = 1, ao_num

    jj = ao_nucl(j)
    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)

    do i = 1, ao_num

      ii = ao_nucl(i)
      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      ao_deriv2_cgtos_x(i,j) = 0.d0
      ao_deriv2_cgtos_y(i,j) = 0.d0
      ao_deriv2_cgtos_z(i,j) = 0.d0

      do n = 1, ao_prim_num(j)

        alpha = ao_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        phiA = ao_expo_phase_ord_transp(4,n,j)
        do m = 1, 3
          A_center(m) = dcmplx(nucl_coord(jj,m))
        enddo

        do l = 1, ao_prim_num(i)

          beta = ao_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          phiB = ao_expo_phase_ord_transp(4,l,i)
          do m = 1, 3
            B_center(m) = dcmplx(nucl_coord(ii,m))
          enddo

          c = ao_coef_cgtos_norm_ord_transp(n,j) * ao_coef_cgtos_norm_ord_transp(l,i)

          C1 = zexp((0.d0, 1.d0) * (-phiA - phiB))
          C2 = zexp((0.d0, 1.d0) * (-phiA + phiB))

          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_x0_1, overlap_y0_1, overlap_z0_1, overlap, dim_1)
          call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                     A_center, conjg(B_center), overlap_x0_2, overlap_y0_2, overlap_z0_2, overlap, dim_1)

          ! ---

          power_A(1) = power_A(1) - 2
          if(power_A(1) > -1) then
            call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                       A_center, B_center, overlap_m2_1, overlap_y, overlap_z, overlap, dim_1)
            call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                       A_center, conjg(B_center), overlap_m2_2, overlap_y, overlap_z, overlap, dim_1)
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(1) = power_A(1) + 2

          power_A(1) = power_A(1) + 2
          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_p2_1, overlap_y, overlap_z, overlap, dim_1)
          call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                     A_center, conjg(B_center), overlap_p2_2, overlap_y, overlap_z, overlap, dim_1)
          power_A(1) = power_A(1) - 2

          aa = dble(power_A(1))

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_x0_1 &
                      + 4.d0 * alpha * alpha * overlap_p2_1
          deriv_tmp_1 = deriv_tmp_1 * overlap_y0_1 * overlap_z0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_x0_2 &
                      + 4.d0 * alpha * alpha * overlap_p2_2
          deriv_tmp_2 = deriv_tmp_2 * overlap_y0_2 * overlap_z0_2

          deriv_tmp = 0.5d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_deriv2_cgtos_x(i,j) += c * deriv_tmp

          ! ---

          power_A(2) = power_A(2) - 2
          if(power_A(2) > -1) then
            call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                       A_center, B_center, overlap_x, overlap_m2_1, overlap_z, overlap, dim_1)
            call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                       A_center, conjg(B_center), overlap_x, overlap_m2_2, overlap_z, overlap, dim_1)
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(2) = power_A(2) + 2

          power_A(2) = power_A(2) + 2
          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_x, overlap_p2_1, overlap_z, overlap, dim_1)
          call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                     A_center, conjg(B_center), overlap_x, overlap_p2_2, overlap_z, overlap, dim_1)
          power_A(2) = power_A(2) - 2

          aa = dble(power_A(2))

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_y0_1 &
                      + 4.d0 * alpha * alpha * overlap_p2_1
          deriv_tmp_1 = deriv_tmp_1 * overlap_x0_1 * overlap_z0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_y0_2 &
                      + 4.d0 * alpha * alpha * overlap_p2_2
          deriv_tmp_2 = deriv_tmp_2 * overlap_x0_2 * overlap_z0_2

          deriv_tmp = 0.5d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_deriv2_cgtos_y(i,j) += c * deriv_tmp

          ! ---

          power_A(3) = power_A(3) - 2
          if(power_A(3) > -1) then
            call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                       A_center, B_center, overlap_x, overlap_y, overlap_m2_1, overlap, dim_1)
            call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                       A_center, conjg(B_center), overlap_x, overlap_y, overlap_m2_2, overlap, dim_1)
          else
            overlap_m2_1 = (0.d0, 0.d0)
            overlap_m2_2 = (0.d0, 0.d0)
          endif
          power_A(3) = power_A(3) + 2

          power_A(3) = power_A(3) + 2
          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_x, overlap_y, overlap_p2_1, overlap, dim_1)
          call overlap_cgaussian_xyz(A_center, conjg(B_center), alpha, conjg(beta), power_A, power_B, &
                                     A_center, conjg(B_center), overlap_x, overlap_y, overlap_p2_2, overlap, dim_1)
          power_A(3) = power_A(3) - 2

          aa = dble(power_A(3))

          deriv_tmp_1 = aa * (aa - 1.d0) * overlap_m2_1            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_z0_1 &
                      + 4.d0 * alpha * alpha * overlap_p2_1
          deriv_tmp_1 = deriv_tmp_1 * overlap_x0_1 * overlap_y0_1

          deriv_tmp_2 = aa * (aa - 1.d0) * overlap_m2_2            &
                      - 4.d0 * alpha * (aa + 0.5d0) * overlap_z0_2 &
                      + 4.d0 * alpha * alpha * overlap_p2_2
          deriv_tmp_2 = deriv_tmp_2 * overlap_x0_2 * overlap_y0_2

          deriv_tmp = 0.5d0 * real(C1 * deriv_tmp_1 + C2 * deriv_tmp_2)

          ao_deriv2_cgtos_z(i,j) += c * deriv_tmp

          ! ---

        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_kinetic_integrals_cgtos, (ao_num, ao_num)]

  BEGIN_DOC
  ! 
  ! Kinetic energy integrals in the cgtos |AO| basis.
  !
  ! $\langle \chi_i |\hat{T}| \chi_j \rangle$
  !
  END_DOC

  implicit none

  integer :: i, j

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP PRIVATE(i, j)             &
  !$OMP SHARED(ao_num, ao_kinetic_integrals_cgtos, ao_deriv2_cgtos_x, ao_deriv2_cgtos_y, ao_deriv2_cgtos_z)
  do j = 1, ao_num
    do i = 1, ao_num
      ao_kinetic_integrals_cgtos(i,j) = -0.5d0 * (ao_deriv2_cgtos_x(i,j) + &
                                                  ao_deriv2_cgtos_y(i,j) + &
                                                  ao_deriv2_cgtos_z(i,j))
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_integrals_n_e_cgtos, (ao_num, ao_num)]

  BEGIN_DOC
  !
  !  Nucleus-electron interaction, in the cgtos |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  END_DOC

  implicit none

  integer              :: power_A(3), power_B(3)
  integer              :: i, j, k, l, m, n, ii, jj
  double precision     :: c, Z, C_center(3)
  double precision     :: phiA, phiB
  complex*16           :: alpha, alpha_inv, A_center(3)
  complex*16           :: beta, beta_inv, B_center(3)
  complex*16           :: C1, C2, I1, I2

  complex*16, external :: NAI_pol_mult_cgtos



  ao_integrals_n_e_cgtos = 0.d0

  !$OMP PARALLEL                                                             &
  !$OMP DEFAULT (NONE)                                                       &
  !$OMP PRIVATE (i, j, k, l, m, n, ii, jj, C_center, Z, c, C1, C2, I1, I2,   &
  !$OMP          alpha, alpha_inv, A_center, phiA, power_A,                  &
  !$OMP          beta, beta_inv, B_center, phiB, power_B)                    &
  !$OMP SHARED (ao_num, ao_prim_num, ao_nucl, nucl_coord,                    &
  !$OMP         ao_power, nucl_num, nucl_charge, n_pt_max_integrals,         &
  !$OMP         ao_expo_cgtos_ord_transp, ao_coef_cgtos_norm_ord_transp,     &
  !$OMP         ao_integrals_n_e_cgtos, ao_expo_phase_ord_transp)
  !$OMP DO SCHEDULE (dynamic)

  do j = 1, ao_num

    jj = ao_nucl(j)
    power_A(1:3) = ao_power(j,1:3)

    do i = 1, ao_num

      ii = ao_nucl(i)
      power_B(1:3) = ao_power(i,1:3)

      do n = 1, ao_prim_num(j)

        alpha = ao_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        phiA = ao_expo_phase_ord_transp(4,n,j)
        do m = 1, 3
          A_center(m) = dcmplx(nucl_coord(jj,m))
        enddo

        do l = 1, ao_prim_num(i)

          beta = ao_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          phiB = ao_expo_phase_ord_transp(4,l,i)
          do m = 1, 3
            B_center(m) = dcmplx(nucl_coord(ii,m))
          enddo

          C1 = zexp((0.d0, 1.d0) * (-phiA - phiB))
          C2 = zexp((0.d0, 1.d0) * ( phiA - phiB))

          c = 0.d0
          do k = 1, nucl_num

            Z = nucl_charge(k)

            C_center(1:3) = nucl_coord(k,1:3)

            I1 = NAI_pol_mult_cgtos(A_center, B_center, power_A, power_B, alpha, beta, &
                                    A_center, B_center, C_center, n_pt_max_integrals)

            I2 = NAI_pol_mult_cgtos(conjg(A_center), B_center, power_A, power_B, conjg(alpha), beta, &
                                    conjg(A_center), B_center, C_center, n_pt_max_integrals)

            c = c - Z * 0.5d0 * real(C1 * I1 + C2 * I2)
          enddo

          ao_integrals_n_e_cgtos(i,j) += c * ao_coef_cgtos_norm_ord_transp(n,j) &
                                           * ao_coef_cgtos_norm_ord_transp(l,i)
        enddo
      enddo
    enddo
  enddo

  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER

! ---

complex*16 function NAI_pol_mult_cgtos(Ae_center, Be_center, power_A, power_B, alpha, beta, &
                                       Ap_center, Bp_center, C_center, n_pt_in)

  BEGIN_DOC
  !
  ! Computes the electron-nucleus attraction with two primitves cgtos.
  !
  ! :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: C_center(3)
  complex*16,       intent(in) :: alpha, Ae_center(3), Ap_center(3)
  complex*16,       intent(in) :: beta, Be_center(3), Bp_center(3)

  integer                      :: i, n_pt, n_pt_out
  double precision             :: dist_AB, dist_AC
  complex*16                   :: p, p_inv, rho, dist, dist_integral, const, const_factor, coeff, factor
  complex*16                   :: P_center(3)
  complex*16                   :: d(0:n_pt_in)

  complex*16, external         :: V_n_e_cgtos
  complex*16, external         :: crint_sum
  complex*16, external         :: crint_1



  dist_AB = 0.d0
  dist_AC = 0.d0
  do i = 1, 3
    dist_AB += abs(Ae_center(i) - Be_center(i))
    dist_AC += abs(Ae_center(i) - C_center(i) * (1.d0, 0.d0))
  enddo

  if((dist_AB .gt. 1d-13) .or. (dist_AC .gt. 1d-13)) then

    continue

  else

    NAI_pol_mult_cgtos = V_n_e_cgtos(power_A(1), power_A(2), power_A(3), &
                                     power_B(1), power_B(2), power_B(3), &
                                     alpha, beta)
    return

  endif

  p     = alpha + beta
  p_inv = (1.d0, 0.d0) / p
  rho   = alpha * beta * p_inv

  dist = (Ae_center(1) - Be_center(1)) * (Ae_center(1) - Be_center(1)) &
       + (Ae_center(2) - Be_center(2)) * (Ae_center(2) - Be_center(2)) &
       + (Ae_center(3) - Be_center(3)) * (Ae_center(3) - Be_center(3))

  const_factor = dist * rho
  if(real(const_factor) > 80.d0) then
    NAI_pol_mult_cgtos = (0.d0, 0.d0)
    return
  endif

  P_center(1) = (alpha * Ae_center(1) + beta * Be_center(1)) * p_inv
  P_center(2) = (alpha * Ae_center(2) + beta * Be_center(2)) * p_inv
  P_center(3) = (alpha * Ae_center(3) + beta * Be_center(3)) * p_inv

  dist_integral = (P_center(1) - C_center(1)) * (P_center(1) - C_center(1)) &
                + (P_center(2) - C_center(2)) * (P_center(2) - C_center(2)) &
                + (P_center(3) - C_center(3)) * (P_center(3) - C_center(3))

  const = p * dist_integral
  factor = zexp(-const_factor)
  coeff = dtwo_pi * factor * p_inv

  n_pt = 2 * ((power_A(1) + power_B(1)) + (power_A(2) + power_B(2)) + (power_A(3) + power_B(3)))
  if(n_pt == 0) then
    NAI_pol_mult_cgtos = coeff * crint_1(0, const)
    return
  endif

  d(0:n_pt_in) = (0.d0, 0.d0)
  call give_cpolynomial_mult_center_one_e(Ap_center, Bp_center, alpha, beta, &
                                          power_A, power_B, C_center, n_pt_in, d, n_pt_out)

  if(n_pt_out < 0) then
    NAI_pol_mult_cgtos = (0.d0, 0.d0)
    return
  endif

  NAI_pol_mult_cgtos = coeff * crint_sum(n_pt_out, const, d)

  return
end

! ---

subroutine give_cpolynomial_mult_center_one_e(A_center, B_center, alpha, beta, &
                                              power_A, power_B, C_center, n_pt_in, d, n_pt_out)

  BEGIN_DOC
  ! Returns the explicit polynomial in terms of the "t" variable of the following
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  END_DOC

  implicit none

  integer,          intent(in) :: n_pt_in, power_A(3), power_B(3)
  double precision, intent(in) :: C_center(3)
  complex*16,       intent(in) :: alpha, beta, A_center(3), B_center(3)
  integer,         intent(out) :: n_pt_out
  complex*16,    intent(inout) :: d(0:n_pt_in)

  integer                      :: a_x, b_x, a_y, b_y, a_z, b_z
  integer                      :: n_pt1, n_pt2, n_pt3, dim, i, n_pt_tmp
  complex*16                   :: p, P_center(3), rho, p_inv, p_inv_2
  complex*16                   :: R1x(0:2), B01(0:2), R1xp(0:2), R2x(0:2)
  complex*16                   :: d1(0:n_pt_in), d2(0:n_pt_in), d3(0:n_pt_in)

  ASSERT (n_pt_in > 1)

  p       = alpha + beta
  p_inv   = (1.d0, 0.d0) / p
  p_inv_2 = 0.5d0 * p_inv

  P_center(1) = (alpha * A_center(1) + beta * B_center(1)) * p_inv
  P_center(2) = (alpha * A_center(2) + beta * B_center(2)) * p_inv
  P_center(3) = (alpha * A_center(3) + beta * B_center(3)) * p_inv

  do i = 0, n_pt_in
    d(i)  = (0.d0, 0.d0)
    d1(i) = (0.d0, 0.d0)
    d2(i) = (0.d0, 0.d0)
    d3(i) = (0.d0, 0.d0)
  enddo

  ! ---

  n_pt1 = n_pt_in

  R1x(0)  =  (P_center(1) - A_center(1))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(1) - C_center(1))

  R1xp(0) =  (P_center(1) - B_center(1))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(1) - C_center(1))

  R2x(0)  =  p_inv_2
  R2x(1)  = (0.d0, 0.d0)
  R2x(2)  = -p_inv_2

  a_x = power_A(1)
  b_x = power_B(1)

  call I_x1_pol_mult_one_e_cgtos(a_x, b_x, R1x, R1xp, R2x, d1, n_pt1, n_pt_in)

  if(n_pt1 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt2 = n_pt_in

  R1x(0)  =  (P_center(2) - A_center(2))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(2) - C_center(2))

  R1xp(0) =  (P_center(2) - B_center(2))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(2) - C_center(2))

  a_y = power_A(2)
  b_y = power_B(2)

  call I_x1_pol_mult_one_e_cgtos(a_y, b_y, R1x, R1xp, R2x, d2, n_pt2, n_pt_in)

  if(n_pt2 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt3 = n_pt_in

  R1x(0)  =  (P_center(3) - A_center(3))
  R1x(1)  =  (0.d0, 0.d0)
  R1x(2)  = -(P_center(3) - C_center(3))

  R1xp(0) =  (P_center(3) - B_center(3))
  R1xp(1) =  (0.d0, 0.d0)
  R1xp(2) = -(P_center(3) - C_center(3))

  a_z = power_A(3)
  b_z = power_B(3)

  call I_x1_pol_mult_one_e_cgtos(a_z, b_z, R1x, R1xp, R2x, d3, n_pt3, n_pt_in)

  if(n_pt3 < 0) then
    n_pt_out = -1
    do i = 0, n_pt_in
      d(i) = (0.d0, 0.d0)
    enddo
    return
  endif

  ! ---

  n_pt_tmp = 0
  call multiply_cpoly(d1, n_pt1, d2, n_pt2, d, n_pt_tmp)

  n_pt_out = 0
  d1(0:n_pt_tmp) = (0.d0, 0.d0)
  call multiply_cpoly(d, n_pt_tmp, d3, n_pt3, d1, n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end

! ---

recursive subroutine I_x1_pol_mult_one_e_cgtos(a, c, R1x, R1xp, R2x, d, nd, n_pt_in)

  BEGIN_DOC
  !  Recursive routine involved in the electron-nucleus potential
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: a, c, n_pt_in
  complex*16, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(inout) :: d(0:n_pt_in)

  integer                   :: nx, ix, dim, iy, ny
  complex*16                :: X(0:max_dim)
  complex*16                :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y

  dim = n_pt_in

  if((a == 0) .and. (c == 0)) then

    nd   = 0
    d(0) = (1.d0, 0.d0)
    return

  elseif((c < 0) .or. (nd < 0)) then

    nd = -1
    return

  elseif((a == 0) .and. (c .ne. 0)) then

    call I_x2_pol_mult_one_e_cgtos(c, R1x, R1xp, R2x, d, nd, n_pt_in)

  elseif(a == 1) then

    nx = nd
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x2_pol_mult_one_e_cgtos(c-1, R1x, R1xp, R2x, X, nx, n_pt_in)

    do ix = 0, nx
      X(ix) *= dble(c)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    call I_x2_pol_mult_one_e_cgtos(c, R1x, R1xp, R2x, Y, ny, n_pt_in)
    call multiply_cpoly(Y, ny, R1x, 2, d, nd)

  else

    nx = 0
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cgtos(a-2, c, R1x, R1xp, R2x, X, nx, n_pt_in)

    do ix = 0, nx
      X(ix) *= dble(a-1)
    enddo
    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    nx = nd
    do ix = 0, n_pt_in
      X(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cgtos(a-1, c-1, R1x, R1xp, R2x, X, nx, n_pt_in)
    do ix = 0, nx
      X(ix) *= dble(c)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    call I_x1_pol_mult_one_e_cgtos(a-1, c, R1x, R1xp, R2x, Y, ny, n_pt_in)
    call multiply_cpoly(Y, ny, R1x, 2, d, nd)

  endif

end

! ---

recursive subroutine I_x2_pol_mult_one_e_cgtos(c, R1x, R1xp, R2x, d, nd, dim)

  BEGIN_DOC
  !  Recursive routine involved in the electron-nucleus potential
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer,    intent(in)    :: dim, c
  complex*16, intent(in)    :: R1x(0:2), R1xp(0:2), R2x(0:2)
  integer,    intent(inout) :: nd
  complex*16, intent(out)   :: d(0:max_dim)

  integer                   :: i, nx, ix, ny
  complex*16                :: X(0:max_dim), Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y

  if(c == 0) then

    nd = 0
    d(0) = (1.d0, 0.d0)
    return

  elseif((nd < 0) .or. (c < 0)) then

    nd = -1
    return

  else

    nx = 0
    do ix = 0, dim
      X(ix) = (0.d0, 0.d0)
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cgtos(0, c-2, R1x, R1xp, R2x, X, nx, dim)

    do ix = 0, nx
      X(ix) *= dble(c-1)
    enddo

    call multiply_cpoly(X, nx, R2x, 2, d, nd)

    ny = 0
    do ix = 0, dim
      Y(ix) = (0.d0, 0.d0)
    enddo

    call I_x1_pol_mult_one_e_cgtos(0, c-1, R1x, R1xp, R2x, Y, ny, dim)

    if(ny .ge. 0) then
      call multiply_cpoly(Y, ny, R1xp, 2, d, nd)
    endif

  endif

end

! ---

complex*16 function V_n_e_cgtos(a_x, a_y, a_z, b_x, b_y, b_z, alpha, beta)

  BEGIN_DOC
  ! Primitve nuclear attraction between the two primitves centered on the same atom.
  !
  ! $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
  !
  ! $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$
  END_DOC

  implicit none

  integer,    intent(in) :: a_x, a_y, a_z, b_x, b_y, b_z
  complex*16, intent(in) :: alpha, beta

  double precision       :: V_phi, V_theta
  complex*16             :: V_r_cgtos

  if( (iand(a_x + b_x, 1) == 1) .or. &
      (iand(a_y + b_y, 1) == 1) .or. &
      (iand(a_z + b_z, 1) == 1) ) then

    V_n_e_cgtos = (0.d0, 0.d0)

  else

    V_n_e_cgtos = V_r_cgtos(a_x + b_x + a_y + b_y + a_z + b_z + 1, alpha + beta) &
                * V_phi(a_x + b_x, a_y + b_y)                                    &
                * V_theta(a_z + b_z, a_x + b_x + a_y + b_y + 1)
  endif

end

! ---

complex*16 function V_r_cgtos(n, alpha)

  BEGIN_DOC
  ! Computes the radial part of the nuclear attraction integral:
  !
  ! $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
  !
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer   , intent(in) :: n
  complex*16, intent(in) :: alpha

  double precision       :: fact

  if(iand(n, 1) .eq. 1) then
    V_r_cgtos = 0.5d0 * fact(shiftr(n, 1)) / (alpha**(shiftr(n, 1) + 1))
  else
    V_r_cgtos = sqpi * fact(n) / fact(shiftr(n, 1)) * (0.5d0/zsqrt(alpha))**(n+1)
  endif

end

! ---

