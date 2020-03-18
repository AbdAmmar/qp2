


 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_alpha_alpha_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! all_states_act_two_rdm_alpha_alpha_mo(i,j,k,l,istate) = STATE SPECIFIC physicist notation for 2RDM of alpha electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \alpha} a_{l \alpha} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 1 
 all_states_act_two_rdm_alpha_alpha_mo = 0.D0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_alpha_alpha_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_beta_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! all_states_act_two_rdm_beta_beta_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of beta electrons 
! 
! <Psi| a^{\dagger}_{i \beta} a^{\dagger}_{j \beta} a_{l \beta} a_{k \beta} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 2
 all_states_act_two_rdm_beta_beta_mo = 0.d0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_beta_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_alpha_beta_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 double precision, allocatable :: state_weights(:) 
 BEGIN_DOC
! all_states_act_two_rdm_alpha_beta_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM of alpha/beta electrons 
! 
! <Psi| a^{\dagger}_{i \alpha} a^{\dagger}_{j \beta} a_{l \beta} a_{k \alpha} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
 END_DOC 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 print*,''
 print*,''
 print*,''
 print*,'providint all_states_act_two_rdm_alpha_beta_mo '
 ispin = 3 
 print*,'ispin = ',ispin
 all_states_act_two_rdm_alpha_beta_mo = 0.d0
 call orb_range_all_states_two_rdm(all_states_act_two_rdm_alpha_beta_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 


 BEGIN_PROVIDER [double precision, all_states_act_two_rdm_spin_trace_mo, (n_act_orb,n_act_orb,n_act_orb,n_act_orb,N_states)]
 implicit none
 BEGIN_DOC
! all_states_act_two_rdm_spin_trace_mo(i,j,k,l,istate) =  STATE SPECIFIC physicist notation for 2RDM 
! 
! \sum_{\sigma, \sigma'} <Psi| a^{\dagger}_{i \sigma} a^{\dagger}_{j \sigma'} a_{l \sigma'} a_{k \sigma} |Psi>
!
! !!!!! WARNING !!!!! ALL SLATER DETERMINANTS IN PSI_DET MUST BELONG TO AN ACTIVE SPACE DEFINED BY "list_act" 
!
! The active part of the two-electron energy for the state istate can be computed as:
!
!   \sum_{i,j,k,l = 1, n_act_orb} all_states_act_two_rdm_spin_trace_mo(i,j,k,l,istate) * < ii jj | kk ll > 
!
!   with ii = list_act(i), jj = list_act(j), kk = list_act(k), ll = list_act(l)
 END_DOC
 double precision, allocatable :: state_weights(:) 
 allocate(state_weights(N_states))
 state_weights = 1.d0/dble(N_states)
 integer :: ispin
 ! condition for alpha/beta spin
 ispin = 4 
 all_states_act_two_rdm_spin_trace_mo = 0.d0
 integer :: i

 call orb_range_all_states_two_rdm(all_states_act_two_rdm_spin_trace_mo,n_act_orb,n_act_orb,list_act,list_act_reverse,ispin,psi_coef,size(psi_coef,2),size(psi_coef,1))

 END_PROVIDER 
