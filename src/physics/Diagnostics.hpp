// =============================================================================
//  SolarSim - Diagnostics.hpp
//  Conserved quantities of the N-body problem, used to judge integrators:
//      E = sum_i 1/2 m_i v_i^2 - sum_{i<j} G m_i m_j / r_ij
//      L = sum_i m_i r_i x v_i ,   P = sum_i m_i v_i
//  Masses are represented by GM, so the quantities are scaled by G (harmless
//  for relative errors). Only massive bodies contribute.
// =============================================================================
#pragma once

#include "NBodySystem.hpp"

#include "../core/Vec3.hpp"

namespace solsim::physics {

struct Invariants {
    double kinetic = 0.0;
    double potential = 0.0;
    double energy = 0.0;
    Vec3 angularMomentum;
    Vec3 momentum;
    Vec3 centerOfMass;
};

Invariants computeInvariants(const NBodySystem& sys);

struct InvariantErrors {
    double energyRel = 0.0;          // (E - E0) / |E0|
    double angularMomentumRel = 0.0; // |L - L0| / |L0|
    double momentumAbs = 0.0;        // |P - P0| (relative to a typical momentum)
    double comDrift = 0.0;           // AU
    double virialRatio = 0.0;        // 2T / |U|  (~1 for bound orbits on average)
};

InvariantErrors compareInvariants(const Invariants& reference, const Invariants& current);

} // namespace solsim::physics
