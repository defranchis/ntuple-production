#ifndef ALEPH_UNITS_H
#define ALEPH_UNITS_H
// Unit conventions shared by the ALEPH analyzers: lengths cm, momenta GeV,
// magnetic field tesla. Dependency-free, so standalone builds can include it.
namespace FCCAnalyses {
namespace AlephUnits {
// pT [GeV] = kPtPerTeslaCm * Bz [T] / |omega [1/cm]|  (0.29979 GeV/T/m in cm)
constexpr double kPtPerTeslaCm = 0.0029979;
// LEP1 beam energy [GeV]: the single source for every xE = E/E_beam branch.
constexpr double kEBeam = 45.6;
}  // namespace AlephUnits

// PDG 2024 central values [GeV]: the single source for every analyzer mass.
// The charged-pion constant cannot be called M_PI, that name is a <cmath> macro.
namespace AlephMasses {
constexpr double kPiCh   = 0.13957039;
constexpr double kK      = 0.493677;
constexpr double kProton = 0.93827208;
constexpr double kKs     = 0.497611;
constexpr double kLambda = 1.115683;
constexpr double kPhi    = 1.019461;   // phi(1020); Gamma = 4.25 MeV
constexpr double kD0     = 1.86484;
constexpr double kDstar  = 2.01026;
constexpr double kDeltaMDstarD0 = 0.145426;  // m(D*+) - m(D0)
}  // namespace AlephMasses
}  // namespace FCCAnalyses
#endif
