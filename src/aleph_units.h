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
// Solenoid field: MC is generated at the nominal value; data runs derive it
// from the magnet current as in ALEPHLIB ALFIEL (nominal 15 kG at 4963750 mA,
// x1.011 when a compensation coil carries less than 800000 mA).
constexpr double kBzNominal = 1.5;   // [T]
constexpr double kGaussPerTesla = 1e4;
constexpr double kFieldNominalkG = kBzNominal * kGaussPerTesla / 1e3;
constexpr double kCurrentNominalmA = 4963750.;
constexpr double kCurrent92OffsetmA = 17700.;   // 1992-93 readout offset
constexpr double kCompCoilMinmA = 800000.;
constexpr double kCompCoilOffCorr = 1.011;
}  // namespace AlephUnits
}  // namespace FCCAnalyses
#endif
