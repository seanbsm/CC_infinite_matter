program chipot_test
! NNLO_opt chiral interaction
use chiral_constants_nnlo_opt , only : init_chp_constants_nnlo_opt
use chiral_potentials_nnlo_opt

implicit none
integer,parameter :: wp=selected_real_kind(15)
! energy scale
real(wp) :: e_scale
! total number of particle
integer:: npart
! baryon density
real(wp) :: rho
! array of single-particle states occupied in HF
type(t_vector), allocatable :: hf_states(:)
! matrix element
complex*16 :: matel
! output values
real(wp) :: oneb,ehf,temp
! auxiliary variables
integer:: ii,jj


npart =14
rho=0.08_wp
write(6,*) "Density: ",rho
write(6,*) "Number of particles: ",npart

! initialize potential
call init_chp_constants_nnlo_opt(rho,npart,e_scale)
! assign states for HF
allocate(hf_states(npart))
call init_14(hf_states,npart)
! calculate one body energy (per particle)
oneb=oneb_ene(hf_states,npart)*e_scale
! write(6,*) "escale: ", e_scale
! write(6,*) "==> One-body energy per particle     [MeV]: ",oneb/real(npart)
! calculate two-body contribution to HF
ehf=0._wp
temp=0._wp
do ii=1,npart
    do jj=1,npart
        !if ( ii /= jj) cycle
        !ehf = ehf + 0.5*twob_as(hf_states(ii),hf_states(jj),hf_states(ii),hf_states(jj))
        temp = twob_as(hf_states(ii),hf_states(jj),hf_states(ii),hf_states(jj))
        ehf = ehf + 0.5*temp
        ! write(6,*) "element: ", temp, &
        ! ", state : ", ii-1," ",jj-1
        ! write(6,*)
    enddo
enddo
! restore dimensions
write(6,*) "==> One-body energy per particle     [MeV]: ",oneb/real(npart)
write(6,*) "==> Two-body energy per particle     [MeV]: ",ehf*e_scale/real(npart)
!write(6,*) e_scale
ehf = ehf*e_scale
write(6,*) "==> Hartree-Fock energy per particle [MeV]: ",(oneb+ehf)/real(npart)

contains

    ! returns total one-body energy
    function oneb_ene(state_vector,nstates) result(obe)
    implicit none
    type(t_vector), intent(in) :: state_vector(:)
    integer, intent(in) :: nstates
    real(wp) :: obe
    integer :: ii
    
    obe=0._wp
    do ii=1,nstates
        obe=obe+0.5_wp*dot_product(state_vector(ii)%k,state_vector(ii)%k)
        ! write(6,*) "element: ", 0.5_wp*dot_product(state_vector(ii)%k,state_vector(ii)%k), ", state: ", ii
    enddo
    
    end function oneb_ene
    
    ! returns anti-symmetrized matrix element: <ij||ab> = <ij|ab> - <ij|ab>
    function twob_as(stateP,stateQ,stateR,stateS) result(matel)
    implicit none
    type(t_vector), intent(in) :: stateP,stateQ,stateR,stateS
    real(wp) :: matel
    matel=0._wp
    ! check isospin conservation
    if((stateP%t+stateQ%t)/=(stateR%t+stateS%t)) return
    ! check momentum conservation
    if((stateP%k(1)+stateQ%k(1))/=(stateR%k(1)+stateS%k(1))) return
    if((stateP%k(2)+stateQ%k(2))/=(stateR%k(2)+stateS%k(2))) return
    if((stateP%k(3)+stateQ%k(3))/=(stateR%k(3)+stateS%k(3))) return
    ! write(6,*) "first state:  ", stateP
    ! write(6,*) "second state: ", stateQ
    matel=chiral_pot_nnlo_opt_np(stateP,stateQ,stateR,stateS)
    matel = matel - chiral_pot_nnlo_opt_np(stateP,stateQ,stateS,stateR)
    
    ! write(6,*) "element: ", matel, ", state : ", stateP,stateQ,stateR,stateS
    
    end function twob_as
    
    ! brute force way of defining the HF state for 14 neutrons
    ! -- NOTE: neutrons have t=1 and protons t=-1 (sorry for the weird convention)
    subroutine init_14(state_vector,npart)
    implicit none
    integer, intent(in) :: npart
    type(t_vector), intent(inout) :: state_vector(npart)
    
!    if(npart.ne.14) stop "ERROR: only 14 particles allowed!"
    ! spin-up first
    !----k (0,0,0) -----------
    state_vector(1)%k(:)=0._wp
    state_vector(1)%s=1
    state_vector(1)%t=1
    !----k (1,0,0) -----------
    state_vector(2)%k(:)=0._wp
    state_vector(2)%k(1)=1._wp
    state_vector(2)%s=1
    state_vector(2)%t=1
    !----k (-1,0,0) -----------
    state_vector(3)%k(:)=0._wp
    state_vector(3)%k(1)=-1._wp
    state_vector(3)%s=1
    state_vector(3)%t=1
    !----k (0,1,0) -----------
    state_vector(4)%k(:)=0._wp
    state_vector(4)%k(2)=1._wp
    state_vector(4)%s=1
    state_vector(4)%t=1
    !----k (0,-1,0) -----------
    state_vector(5)%k(:)=0._wp
    state_vector(5)%k(2)=-1._wp
    state_vector(5)%s=1
    state_vector(5)%t=1
    !----k (0,0,1) -----------
    state_vector(6)%k(:)=0._wp
    state_vector(6)%k(3)=1._wp
    state_vector(6)%s=1
    state_vector(6)%t=1
    !----k (0,0,-1) -----------
    state_vector(7)%k(:)=0._wp
    state_vector(7)%k(3)=-1._wp
    state_vector(7)%s=1
    state_vector(7)%t=1
    ! spin-down now
    !----k (0,0,0) -----------
    state_vector(8)%k(:)=0._wp
    state_vector(8)%s=-1
    state_vector(8)%t=1
    !----k (1,0,0) -----------
    state_vector(9)%k(:)=0._wp
    state_vector(9)%k(1)=1._wp
    state_vector(9)%s=-1
    state_vector(9)%t=1
    !----k (-1,0,0) -----------
    state_vector(10)%k(:)=0._wp
    state_vector(10)%k(1)=-1._wp
    state_vector(10)%s=-1
    state_vector(10)%t=1
    !----k (0,1,0) -----------
    state_vector(11)%k(:)=0._wp
    state_vector(11)%k(2)=1._wp
    state_vector(11)%s=-1
    state_vector(11)%t=1
    !----k (0,-1,0) -----------
    state_vector(12)%k(:)=0._wp
    state_vector(12)%k(2)=-1._wp
    state_vector(12)%s=-1
    state_vector(12)%t=1
    !----k (0,0,1) -----------
    state_vector(13)%k(:)=0._wp
    state_vector(13)%k(3)=1._wp
    state_vector(13)%s=-1
    state_vector(13)%t=1
    !----k (0,0,-1) -----------
    state_vector(14)%k(:)=0._wp
    state_vector(14)%k(3)=-1._wp
    state_vector(14)%s=-1
    state_vector(14)%t=1
    
    end subroutine init_14

end program
