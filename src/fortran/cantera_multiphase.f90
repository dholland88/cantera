! This file is part of Cantera. See License.txt in the top-level directory or
! at https://cantera.org/license.txt for license and copyright information.

module cantera_multiphase

  use fct
  use cantera_xml
  use cantera_thermo

  type multiphase_t
    integer :: mix_id
    integer :: nsp
    integer :: err
    integer :: nPhases
    type(phase_t) :: phases(10)
  end type multiphase_t

! these definitions are for use with the equilibrate function.
!  integer, parameter :: TV = 100
!  integer, parameter :: HP = 101
!  integer, parameter :: SP = 102
!  integer, parameter :: PV = 103
!  integer, parameter :: TP = 104
!  integer, parameter :: UV = 105
!  integer, parameter :: SV = 107

!  integer, parameter :: VT = -100
!  integer, parameter :: PH = -101
!  integer, parameter :: PS = -102
!  integer, parameter :: VP = -103
!  integer, parameter :: PT = -104
!  integer, parameter :: VU = -105
!  integer, parameter :: VS = -107

contains

    type(multiphase_t) function ctmultiphase_newMultiPhase()
      implicit none
      type(multiphase_t) :: self
      self%mix_id = mphase_addmultiphase()
      self%err = 0
      self%nPhases = 0
      self%nsp = 0
      ctmultiphase_newMultiPhase = self
    end function ctmultiphase_newMultiPhase

    subroutine ctmultiphase_addphase(self, phase, nmole)
      implicit none
      type(multiphase_t), intent(inout) :: self
      type(phase_t) :: phase
      double precision :: nmole
      self%nPhases = self%nPhases + 1
      self%nsp = self%nsp + phase%nsp
      self%phases(self%nPhases) = phase
      self%err = mphase_addphase(self%mix_id, phase%thermo_id, nmole)
    end subroutine ctmultiphase_addphase

    integer function ctmultiphase_elementIndex(self, nm)
      implicit none
      type(multiphase_t), intent(inout) :: self
      character*(*), intent(in) :: nm
      ctmultiphase_elementindex = mphase_elementindex(self%mix_id, nm)
    end function ctmultiphase_elementIndex

    subroutine ctmultiphase_elementName(self, k, nm)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      character*(*), intent(out) :: nm
      self%err = mphase_elementname(self%mix_id, k, nm)
    end subroutine ctmultiphase_elementName

    integer function ctmultiphase_nSpecies(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_nspecies = self%nsp
    end function ctmultiphase_nSpecies

    subroutine ctmultiphase_speciesName(self, k, nm)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      character*(*), intent(out) :: nm
      self%err = mphase_speciesname(self%mix_id, k, nm)
    end subroutine ctmultiphase_speciesname

    double precision function ctmultiphase_nAtoms(self, k, m)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      integer, intent(in) :: m
      ctmultiphase_natoms = mphase_natoms(self%mix_id, k, m)
    end function ctmultiphase_natoms

    subroutine ctmultiphase_getMoleFractions(self, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(out) :: x(self%nsp)
      self%err = mphase_getmolefractions(self%mix_id, x)
    end subroutine ctmultiphase_getMoleFractions

    subroutine ctmultiphase_init(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      self%err = mphase_init(self%mix_id)
    end subroutine ctmultiphase_init

    subroutine ctmultiphase_phaseName(self, p, nm)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in):: p
      character*(*), intent(out) :: nm
      self%err = mphase_phasename(self%mix_id, p, nm)
    end subroutine ctmultiphase_phaseName

    integer function ctmultiphase_phaseIndex(self, name)
      implicit none
      type(multiphase_t), intent(inout) :: self
      character*(*):: name
      ctmultiphase_phaseIndex = mphase_phaseindex(self%mix_id, name)
    end function ctmultiphase_phaseIndex

    double precision function ctmultiphase_phaseMoles(self, p)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: p
      ctmultiphase_phaseMoles = mphase_phasemoles(self%mix_id, p)
    end function ctmultiphase_phaseMoles

    subroutine ctmultiphase_setPhaseMoles(self, n, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: n
      double precision, intent(in) :: x
      self%err = mphase_setphasemoles(self%mix_id, n, x)
    end subroutine ctmultiphase_setPhaseMoles

    type(phase_t) function ctmultiphase_getphase(self, n)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer:: n
      ctmultiphase_getphase = self%phases(n)
    end function ctmultiphase_getphase

    double precision function  ctmultiphase_speciesMoles(self, k)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctmultiphase_speciesMoles = mphase_speciesmoles(self%mix_id, k)
    end function ctmultiphase_speciesMoles

    integer function ctmultiphase_speciesIndex(self, k, p)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k, p
      ctmultiphase_speciesindex = mphase_speciesindex(self%mix_id, k, p)
    end function ctmultiphase_speciesIndex

    integer function ctmultiphase_speciesIndexByName(self, nm, pname)
      implicit none
      type(multiphase_t), intent(inout) :: self
      character*(*), intent(in) :: nm, pname
      ctmultiphase_speciesIndexByName = mphase_speciesindexbyname(self%mix_id, nm, pname)
    end function ctmultiphase_speciesIndexByName

    double precision function ctmultiphase_minTemp(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_minTemp = mphase_mintemp(self%mix_id)
    end function ctmultiphase_minTemp

    double precision function ctmultiphase_maxTemp(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_maxTemp = mphase_maxtemp(self%mix_id)
    end function ctmultiphase_maxTemp

    double precision function ctmultiphase_charge(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_charge = mphase_charge(self%mix_id)
    end function ctmultiphase_charge

    double precision function ctmultiphase_phaseCharge(self, p)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: p
      ctmultiphase_phaseCharge = mphase_phaseCharge(self%mix_id, p)
    end function ctmultiphase_phaseCharge

    double precision function  ctmultiphase_elementMoles(self, n)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: n
      ctmultiphase_elementMoles = mphase_elementmoles(self%mix_id, n)
    end function ctmultiphase_elementMoles

    subroutine ctmultiphase_getChemPotentials(self, mu)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(out) :: mu(self%nsp)
      self%err = mphase_getchempotentials(self%mix_id, mu)
    end subroutine ctmultiphase_getChemPotentials

    double precision function ctmultiphase_temperature(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_temperature = mphase_temperature(self%mix_id)
    end function ctmultiphase_temperature

    subroutine ctmultiphase_equilibrate(self, XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level)
      implicit none
      type(multiphase_t), intent(inout) :: self
      character*(*), intent(in) :: XY
      character*(*), intent(in), optional :: solver
      double precision, intent(in), optional :: rtol
      integer, intent(in), optional :: max_steps
      integer, intent(in), optional :: max_iter
      integer, intent(in), optional :: estimate_equil
      integer, intent(in), optional :: log_level
      character*(50) :: solver_
      double precision :: rtol_
      integer :: max_steps_
      integer :: max_iter_
      integer :: estimate_equil_
      integer :: log_level_
      solver_ = 'auto'
      rtol_ = 1e-9
      max_steps_ = 50000
      max_iter_ = 100
      estimate_equil_ = 0
      log_level_ = 0
      if(present(solver)) then
          solver_ = solver
      endif
      if(present(rtol)) then
          rtol_ = rtol
      endif
      if(present(max_steps)) then
          max_steps_ = max_steps
      endif
      if(present(max_iter)) then
          max_iter_ = max_iter
      endif
      if(present(estimate_equil)) then
          estimate_equil_ = estimate_equil
      endif
      if(present(log_level)) then
          log_level_ = log_level
      endif
      self%err = mphase_equil(self%mix_id, XY, trim(solver_), rtol_, max_steps_, max_iter_, estimate_equil_, log_level_)
    end subroutine ctmultiphase_equilibrate

    subroutine ctmultiphase_setTemperature(self, t)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(in) :: t
      self%err = mphase_settemperature(self%mix_id, t)
    end subroutine ctmultiphase_setTemperature

    subroutine ctmultiphase_setState_TP(self, t, p)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      self%err = mphase_set_tp(self%mix_id, t, p)
    end subroutine ctmultiphase_setState_TP

    subroutine ctmultiphase_setState_TPMoles(self, t, p, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(in) :: t
      double precision, intent(in) :: p
      double precision, intent(in) :: x(self%nsp)
      self%err = mphase_set_tpmoles(self%mix_id, t, p, x)
    end subroutine ctmultiphase_setState_TPMoles

    double precision function ctmultiphase_volume(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_volume = mphase_volume(self%mix_id)
   end function ctmultiphase_volume

    double precision function ctmultiphase_pressure(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_pressure = mphase_pressure(self%mix_id)
    end function ctmultiphase_pressure

    subroutine ctmultiphase_setPressure(self, p)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(in) :: p
      self%err = mphase_setpressure(self%mix_id, p)
    end subroutine ctmultiphase_setpressure

    double precision function ctmultiphase_enthalpy(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_enthalpy = mphase_enthalpy(self%mix_id)
    end function ctmultiphase_enthalpy

    double precision function ctmultiphase_intEnergy(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_intEnergy = mphase_intenergy(self%mix_id)
    end function ctmultiphase_intEnergy

    double precision function ctmultiphase_entropy(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_entropy = mphase_entropy(self%mix_id)
    end function ctmultiphase_entropy

    double precision function ctmultiphase_gibbs(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_gibbs = mphase_gibbs(self%mix_id)
    end function ctmultiphase_gibbs

    double precision function ctmultiphase_cp(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_cp = mphase_cp(self%mix_id)
    end function ctmultiphase_cp

    integer function ctmultiphase_nPhases(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      ctmultiphase_nPhases = self%nPhases
    end function ctmultiphase_nPhases

    integer function ctmultiphase_speciesPhaseIndex(self, k)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctmultiphase_speciesPhaseIndex = mphase_speciesphaseindex(self%mix_id, k)
    end function ctmultiphase_speciesPhaseIndex

    double precision function ctmultiphase_moleFraction(self, k)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      ctmultiphase_molefraction = mphase_molefraction(self%mix_id, k)
    end function ctmultiphase_moleFraction

    subroutine ctmultiphase_setPhaseMoleFractions(self, i, mfs)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: i
      double precision, intent(in) :: mfs(*)
      type(phase_t) :: p
      self%err = mphase_setphasemolefractions(self%mix_id, i, mfs)
    end subroutine ctmultiphase_setPhaseMoleFractions

    subroutine ctmultiphase_setMoles(self, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(in) :: x(*)
      self%err = mphase_setmoles(self%mix_id, x)
    end subroutine ctmultiphase_setMoles

    subroutine ctmultiphase_setMolesByName(self, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      character*(*), intent(in) :: x
      self%err = mphase_setmolesbyname(self%mix_id, x)
    end subroutine ctmultiphase_setMolesByName

    subroutine ctmultiphase_getMoles(self, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(out) :: x(self%nsp)
      self%err = mphase_getmoles(self%mix_id, x)
    end subroutine ctmultiphase_getMoles

    subroutine ctmultiphase_addSpeciesMoles(self, k, x)
      implicit none
      type(multiphase_t), intent(inout) :: self
      integer, intent(in) :: k
      double precision, intent(in) :: x
      self%err = mphase_addSpeciesMoles(self%mix_id, k, x)
    end subroutine ctmultiphase_addSpeciesMoles

    subroutine ctmultiphase_getElemAbundances(self, elemAbundances)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(out) :: elemAbundances(self%nsp)
      self%err = mphase_getelemabundances(self%mix_id, elemAbundances)
    end subroutine ctmultiphase_getElemAbundances

    subroutine ctmultiphase_uploadmolefractions(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      self%err = mphase_uploadmolefractions(self%mix_id)
    end subroutine ctmultiphase_uploadmolefractions

    subroutine ctmultiphase_report(self)
      implicit none
      type(multiphase_t), intent(inout) :: self
      self%err = mphase_report(self%mix_id)
    end subroutine ctmultiphase_report

    subroutine ctmultiphase_getMolecularWeights(self, mwts)
      implicit none
      type(multiphase_t), intent(inout) :: self
      double precision, intent(out) :: mwts(self%nsp)
      type(phase_t) :: p
      integer :: i, pos
      pos = 1
      do i = 1,self%nPhases
         p = ctmultiphase_getphase(self,i)
         call ctthermo_getMolecularWeights(p, mwts(pos:pos+p%nsp-1))
         pos = pos + p%nsp
      end do
    end subroutine ctmultiphase_getMolecularWeights

end module cantera_multiphase
