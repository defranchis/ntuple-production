#ifndef ALEPH_UNITS_H
#define ALEPH_UNITS_H
// Unit conventions shared by the ALEPH analyzers. All lengths are cm, momenta
// GeV, magnetic field tesla. Dependency-free so that standalone builds
// (test/pvnew/bench_compiled.cxx) can include it as well.
namespace FCCAnalyses {
namespace AlephUnits {
// pT [GeV] = kPtPerTeslaCm * Bz [T] / |omega [1/cm]|  (0.29979 GeV/T/m in cm;
// 5-digit value shared with the PV fitter, kept as is)
constexpr double kPtPerTeslaCm = 0.0029979;
}  // namespace AlephUnits
}  // namespace FCCAnalyses
#endif
