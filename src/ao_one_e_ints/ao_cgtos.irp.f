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
 !$OMP        ao_expo_pw_ord_transp, ao_expo_phase_ord_transp,           &
 !$OMP        ao_overlap_cgtos)

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
  !$OMP        ao_expo_pw_ord_transp, ao_expo_phase_ord_transp,           & 
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

