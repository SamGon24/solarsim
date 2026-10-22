// =============================================================================
//  SolarSim - Diagnostics.cpp
// =============================================================================
#include "Diagnostics.hpp"

#include <cmath>

namespace solsim::physics {

Invariants computeInvariants(const NBodySystem& sys) {
    Invariants inv;
    const std::size_t nM = sys.massiveCount();
    const auto& x = sys.positions();
    const auto& v = sys.velocities();
    const auto& gm = sys.gm();
    double total = 0.0;
    for (std::size_t i = 0; i < nM; ++i) {
        inv.kinetic += 0.5 * gm[i] * v[i].norm2();
        inv.angularMomentum += x[i].cross(v[i]) * gm[i];
        inv.momentum += v[i] * gm[i];
        inv.centerOfMass += x[i] * gm[i];
        total += gm[i];
        for (std::size_t j = i + 1; j < nM; ++j) inv.potential -= gm[i] * gm[j] / (x[j] - x[i]).norm();
    }
    if (total > 0.0) inv.centerOfMass /= total;
    inv.energy = inv.kinetic + inv.potential;
    return inv;
}

InvariantErrors compareInvariants(const Invariants& ref, const Invariants& cur) {
    InvariantErrors e;
    e.energyRel = ref.energy != 0.0 ? (cur.energy - ref.energy) / std::fabs(ref.energy) : 0.0;
    const double L0 = ref.angularMomentum.norm();
    e.angularMomentumRel = L0 > 0.0 ? (cur.angularMomentum - ref.angularMomentum).norm() / L0 : 0.0;
    e.momentumAbs = (cur.momentum - ref.momentum).norm();
    e.comDrift = (cur.centerOfMass - ref.centerOfMass).norm();
    e.virialRatio = cur.potential != 0.0 ? 2.0 * cur.kinetic / std::fabs(cur.potential) : 0.0;
    return e;
}

} // namespace solsim::physics
