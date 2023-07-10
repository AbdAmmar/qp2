
 BEGIN_PROVIDER [real*8, gradvec_detail_right_old, (0:3,nMonoEx)]
&BEGIN_PROVIDER [real*8, gradvec_detail_left_old,  (0:3,nMonoEx)]
  BEGIN_DOC
  ! calculate the orbital gradient <Psi| H E_pq |Psi> by hand, i.e. for
  ! each determinant I we determine the string E_pq |I> (alpha and beta
  ! separately) and generate <Psi|H E_pq |I>
  ! sum_I c_I <Psi|H E_pq |I> is then the pq component of the orbital
  ! gradient
  ! E_pq = a^+_pa_q + a^+_Pa_Q
  END_DOC
  implicit none
  integer                        :: ii,tt,aa,indx,ihole,ipart,istate,ll
  real*8                         :: res_l(0:3), res_r(0:3)
  
 do ii = 1, n_core_inact_orb
  ihole = list_core_inact(ii)
  do aa = 1, n_virt_orb
   ipart = list_virt(aa)
   indx = mat_idx_c_v(ii,aa) 
   call calc_grad_elem_h_tc(ihole,ipart,res_l, res_r)
   do ll = 0, 3
    gradvec_detail_left_old (ll,indx)=res_l(ll)
    gradvec_detail_right_old(ll,indx)=res_r(ll)
   enddo
  enddo
 enddo
!  do indx=1,nMonoEx
!    ihole=excit(1,indx)
!    ipart=excit(2,indx)
!    call calc_grad_elem_h_tc(ihole,ipart,res_l, res_r)
!    do ll = 0, 3
!     gradvec_detail_left_old (ll,indx)=res_l(ll)
!     gradvec_detail_right_old(ll,indx)=res_r(ll)
!    enddo
!  end do
  
  real*8                         :: norm_grad_left, norm_grad_right
  norm_grad_left=0.d0
  norm_grad_right=0.d0
  do indx=1,nMonoEx
    norm_grad_left+=gradvec_detail_left_old(0,indx)*gradvec_detail_left_old(0,indx)
    norm_grad_right+=gradvec_detail_right_old(0,indx)*gradvec_detail_right_old(0,indx)
  end do
  norm_grad_left=sqrt(norm_grad_left)
  norm_grad_right=sqrt(norm_grad_right)
!  if (bavard) then
    write(6,*)
    write(6,*) ' Norm of the LEFT  orbital gradient (via <0|EH|0>) : ', norm_grad_left
    write(6,*) ' Norm of the RIGHT orbital gradient (via <0|HE|0>) : ', norm_grad_right
    write(6,*)
!  endif
  
  
END_PROVIDER

subroutine calc_grad_elem_h_tc(ihole,ipart,res_l, res_r)
  BEGIN_DOC
  ! eq 18 of Siegbahn et al, Physica Scripta 1980
  ! we calculate res_l  = <Phi| H^tc E_pq | Psi>, and res_r = <Phi| E_qp H^tc | Psi>
  ! q=hole, p=particle
  ! res_l(0) =   total matrix element
  ! res_l(1) =   one-electron part 
  ! res_l(2) =   two-electron part 
  ! res_l(3) = three-electron part 
  END_DOC
  implicit none
  integer, intent(in)            :: ihole,ipart
  double precision, intent(out)  :: res_l(0:3), res_r(0:3)
  integer                        :: mu,iii,ispin,ierr,nu,istate,ll
  integer(bit_kind), allocatable :: det_mu(:,:),det_mu_ex(:,:)
  real*8                         :: i_H_chi_array(0:3,N_states),i_H_phi_array(0:3,N_states),phase
  allocate(det_mu(N_int,2))
  allocate(det_mu_ex(N_int,2))
  
  res_l=0.D0
  res_r=0.D0
  
!  print*,'in i_h_psi'
!  print*,ihole,ipart
  do mu=1,n_det
    ! get the string of the determinant
    call det_extract(det_mu,mu,N_int)
    do ispin=1,2
      ! do the monoexcitation on it
      call det_copy(det_mu,det_mu_ex,N_int)
      call do_signed_mono_excitation(det_mu,det_mu_ex,nu             &
          ,ihole,ipart,ispin,phase,ierr)
      if (ierr.eq.1) then

        call i_H_tc_psi_phi(det_mu_ex,psi_det,psi_l_coef_bi_ortho,psi_r_coef_bi_ortho,N_int & 
            ,N_det,N_det,N_states,i_H_chi_array,i_H_phi_array)
!        print*,i_H_chi_array(1,1),i_H_phi_array(1,1)
        do istate=1,N_states
         do ll = 0,3
          res_l(ll)+=i_H_chi_array(ll,istate)*psi_r_coef_bi_ortho(mu,istate)*phase
          res_r(ll)+=i_H_phi_array(ll,istate)*psi_l_coef_bi_ortho(mu,istate)*phase
         enddo
        end do
      end if
    end do
  end do
  
  ! state-averaged gradient
  res_l*=1.d0/dble(N_states)
  res_r*=1.d0/dble(N_states)
  
end 

