#ifndef JOINTVTX_TABLES_H
#define JOINTVTX_TABLES_H

// GENERATED FILE, do not edit by hand.
// Generated from the parameter files v0_params.json, score_coefs.json and score_coefs_final.json of the
// reference implementation and from the resolved nominal configurations (V0 class on, and fitter without V0
// class; beam-spot settings S1 and LR) by gen_tables.py. Numbers are printed as repr() of the source double.
// The tables were fitted on simulated Z decays (V0 tables and scores on separate simulation files).
// Included: every table and setting read by one of four configurations: V0 class on (the stage1 module) or
// off, each with the beam-spot settings S1 (stage1) and LR.
// The S1 beam-spot widths, the PV track pre-selection window and the PV chi2 cut are taken from
// aleph_reco_config.h; particle masses and the momentum constant from aleph_units.h through jointvtx_common.h.

#include <array>

#include "aleph_reco_config.h"

namespace FCCAnalyses {
namespace AlephJointVtx {

// ---------------------------------------------------------------- types
// log prior odds of one species: in bins of the transverse distance r from the reference PV; below fsRmax
// also in bins of the 3D flight significance (fs[r bin][significance bin])
struct PriorSpecies {
  std::array<double, 15> r;
  std::array<std::array<double, 4>, 15> fs;
};
struct PriorTable {
  std::array<double, 16> rBins;
  std::array<double, 5> fsBins;
  double fsRmax;
  PriorSpecies k0s, lambda, conversion;
};
// V0-against-heavy-flavour term: b0 + ext[min(nExt, 2) - 1] + jet[min(nJet, 3) - 1] + wdl * dl
struct HfSpecies {
  double b0;
  std::array<double, 2> ext;
  std::array<double, 3> jet;
  double wdl;
};
struct HfTable {
  HfSpecies k0s, lambda, conversion;
};
// everything that depends on the beam-spot setting; pointW and pvWscale: K0S, Lambda
struct BeamSpotBundle {
  std::array<double, 3> sigma;
  PriorTable prior;
  HfTable hf;
  std::array<double, 2> pointW;
  std::array<double, 2> pvWscale;
};
// logistic score: ((x - mu) / sd) . coef + b0 over the inputs named in names
struct ScoreCoefs {
  std::array<const char*, 30> names;
  std::array<double, 30> mu, sd, coef;
  double b0;
};

// ---------------------------------------------------------------- tables common to both beam-spot settings
// K0S / Lambda hit term: rate of wrongly attached inner hits per daughter in bins of the daughter momentum [GeV]
inline constexpr std::array<double, 5> JV_HIT_EPS_P_BINS = {{0.0, 1.0, 2.0, 5.0, 1000000000.0}};
inline constexpr std::array<double, 4> JV_HIT_EPS = {{0.0632, 0.0474, 0.0621, 0.1175}};
// PV pointing: tail parameter of the true-V0 pointing significance
constexpr double JV_POINT_NU_SIGNAL = 1.5;
// conversion hit term: log likelihood ratio for 0 / 1 / 2 daughters with inner hits, in bins of the radius [cm];
// the inner edges are the hit-radius bounds JV_HITR_* (checked in jointvtx_lookup.h)
inline constexpr std::array<double, 6> JV_HIT_TABLE_CONV_R_BINS = {{0.0, 7.6, 13.1, 28.9, 45.0, 1000000000.0}};
inline constexpr std::array<std::array<double, 3>, 5> JV_HIT_TABLE_CONV_LOG_LR =
    {{{{0.0, 0.0, 0.0}}, {{0.3576, -1.3204, -1.0631}}, {{0.7866, -0.8019, -0.3267}}, {{0.9779, -1.0793, -1.3695}}, {{0.8262, -2.514, -2.0252}}}};
// proper decay length c * tau [cm] of K0S and Lambda (lifetime term of the heavy-flavour features): the values
// the heavy-flavour tables were fitted with; change them only together with the tables, through the generator
constexpr double JV_CTAU_K0S = 2.6844;
constexpr double JV_CTAU_LAMBDA = 7.89;
// conversion pointing: -w log(1 + max(1 - cos, 0) / tau)
constexpr double JV_CONV_POINT_W = 2.0;
constexpr double JV_CONV_POINT_TAU = 0.003;

// score of the candidate vertices (V0-class and fitter without V0 class alike) and of the final vertices
inline constexpr ScoreCoefs JV_SCORE_CAND = {
  {{"n", "lp", "lc", "dr", "s0", "s1", "sm", "n2v", "n0v", "lf", "ls", "ndf", "lm", "ld", "lpm", "fsig2n", "fneg", "sumpos", "lmc", "lptr", "jc", "jsame", "s0s", "drs", "lfs", "lps", "ldlf", "n2", "n2lf", "n2s1"}},
  {{3.5040530582166545, 1.9182162359957258, -4.557990491955308, 1.1067901513120912, 0.2165660846956447, 3.5145157073557263, 6.556137302613969, 2.227790059772374, 0.3554409236059936, 3.9602067933095437, -3.408524188482586, 1.0553135165661682, 0.5650683713060602, -1.7675423597754036, -0.2761368069403037, 0.4997047846255688, 0.14593281364959182, 3.1654853565172023, 1.1859816031032866, -0.05851965840664883, -4.071768323268696, 0.9666579514129194, 4.931763577521925, 1.9691836029047454, 18.421654476065775, 4.374142326369965, -5.504651645033777, 0.3023008269876361, 1.3335669022137224, 2.146253638733177}},
  {{1.365704119510151, 0.8334199411655129, 2.0523450370080782, 0.8626697884261985, 2.9365731968444253, 4.535385104307441, 5.633087083096763, 1.248835961853372, 0.6722945302218692, 1.6548161922977016, 0.6096110837643469, 0.9943746316130992, 0.6707616467074713, 1.1109835586071637, 0.7454456868751153, 0.32810118633131374, 0.22647303585253165, 2.3880755776992526, 0.7157982391409344, 0.9412949336653514, 2.27470397290817, 0.10121624161082442, 6.114313520154316, 2.8594162507755345, 15.406169201684243, 3.0334240267640435, 3.189997542166829, 0.4592548714931718, 2.210062993542056, 4.306026446122472}},
  {{1.237854048553705, 1.5732725382474937, 0.026338339267798223, -0.6008663230312173, -0.13378755501653086, -0.14836958521387283, -0.54020721781777, 0.003748481195246905, -0.005611088490994372, -0.1950122848933256, 0.22932111666159194, -0.2344918988863268, -0.5933157252226671, 0.5638322725515829, 0.33416265527295996, 0.585625994406432, -0.019389537786195933, 0.8676892850471551, -0.3547400889966168, -0.22094698779211042, -0.2638910367270631, -0.09163384031817504, 0.26191378697080164, 0.622228372628427, 0.7984733835304749, -0.45890535284764267, -0.7224365872454194, -0.4865904330113737, -0.4999817588985901, -0.10397227757205053}},
  -1.2041572738506732
};
inline constexpr ScoreCoefs JV_SCORE_FINAL = {
  {{"n", "lp", "lc", "dr", "s0", "s1", "sm", "n2v", "n0v", "lf", "ls", "ndf", "lm", "ld", "lpm", "fsig2n", "fneg", "sumpos", "lmc", "lptr", "jc", "jsame", "s0s", "drs", "lfs", "lps", "ldlf", "n2", "n2lf", "n2s1"}},
  {{3.5627759646765216, 2.329039803595606, -5.7545115018910975, 0.6889982589701532, 0.9290275017883128, 4.184415574340118, 8.540583114086846, 2.1970627759646764, 0.37329621808408525, 4.631760303220529, -3.1807608093509088, 1.2965862950054616, 0.4681364116668004, -1.1612901818769847, 0.014333868362776295, 0.6569211712329418, 0.1131187779392815, 4.184738126941893, 1.0420574658650996, -0.2683230219680559, -5.213051971984178, 0.9916873497335188, 6.585082735795478, 0.8701864071279414, 23.56360414369257, 5.849910875422863, -4.57628740999979, 0.2390094067959301, 1.1004933245339348, 1.8453545116774308}},
  {{1.2909985821581509, 0.6522917052134245, 1.9638697991589624, 0.6288623110538938, 3.553861302391609, 4.507056159511376, 5.674256956635017, 1.3767405995800703, 0.7325922469262136, 1.4527218031002949, 0.7778981722280877, 1.104603493027261, 0.6267735809136276, 0.9213721173101143, 0.828136758173645, 0.2961661789975795, 0.20042832521000856, 2.4884180074573887, 0.6863832673168051, 0.9511413468260153, 2.1716544150687263, 0.049333316085875906, 6.978235039966297, 1.72696679321361, 15.217169818466507, 2.8691251680347576, 3.5348091705013776, 0.4264784991754827, 2.0777963370831176, 4.0302155001460775}},
  {{1.4680954657025016, 0.9640090461361243, -1.1713377696545744, 0.3882938628876602, 0.05329502035608497, -0.01340149513869401, -0.3110760631116224, -0.14284246195787134, 0.21854822950776487, 0.5562627401460019, -0.5697613784202592, -0.16425235915738529, -1.2127331709465283, 0.3094824469642576, 0.30167146863767735, 0.5249083895941907, -0.004096564251550579, 0.6460412229505741, -0.5028781397997418, 0.6796188233051904, -0.18252766063077888, 0.02217227345262932, 0.3734324716401859, -0.26398436801246056, -0.12027266182291925, -0.2832502655357434, -0.7890809370362307, -0.37979832942007413, -0.33174479471755053, -0.22374794101804676}},
  0.8087910288383107
};

// ---------------------------------------------------------------- beam-spot bundles
inline constexpr BeamSpotBundle JV_BUNDLE_S1 = {
  {{AlephReco::kBeamSigmaX_um / 1e4, AlephReco::kBeamSigmaY_um / 1e4, AlephReco::kBeamSigmaZ_cm}},  // sigma [cm]
  {
    {{0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 13.0, 20.0, 30.0, 45.0, 70.0, 100.0, 1000000000.0}},  // rBins
    {{0.0, 3.0, 5.0, 8.0, 1000000000.0}},  // fsBins
    2.0,  // fsRmax
    {{{-4.2106, -3.5154, -2.7771, -1.8854, -0.6908, 0.4148, 1.2235, 1.9256, 2.3419, 2.6681, 2.0166, 2.1122, 2.0632, 1.4679, 1.0452}}, {{{{-4.1537, -4.5397, -3.5818, -2.5134}}, {{-3.729, -3.1037, -2.2774, -0.0129}}, {{-3.1515, -2.5241, -2.0251, -0.6502}}, {{-2.6687, -2.2072, -1.6113, -0.7437}}, {{-2.6585, -2.0449, -1.3222, 0.6505}}, {{-1.9479, -1.5753, -0.6116, 1.5758}}, {{1.2235, 1.2235, 1.2235, 1.2235}}, {{1.9256, 1.9256, 1.9256, 1.9256}}, {{2.3419, 2.3419, 2.3419, 2.3419}}, {{2.6681, 2.6681, 2.6681, 2.6681}}, {{2.0166, 2.0166, 2.0166, 2.0166}}, {{2.1122, 2.1122, 2.1122, 2.1122}}, {{2.0632, 2.0632, 2.0632, 2.0632}}, {{1.4679, 1.4679, 1.4679, 1.4679}}, {{1.0452, 1.0452, 1.0452, 1.0452}}}}},  // K0S: r, (r, fs)
    {{{-5.7286, -5.3409, -5.1605, -4.5122, -3.5286, -2.7353, -2.2497, -1.6506, -1.0178, -0.2541, -0.9219, -0.9069, -0.3089, -0.3456, -1.5471}}, {{{{-5.5559, -5.6159, -4.6022, -3.4265}}, {{-5.0412, -6.5516, -4.3204, -3.6632}}, {{-4.9843, -5.1295, -6.0685, -4.8689}}, {{-4.6122, -4.5461, -3.9556, -4.728}}, {{-4.1643, -4.3198, -2.9003, -2.8351}}, {{-3.9399, -3.6597, -3.4754, -1.7972}}, {{-2.2497, -2.2497, -2.2497, -2.2497}}, {{-1.6506, -1.6506, -1.6506, -1.6506}}, {{-1.0178, -1.0178, -1.0178, -1.0178}}, {{-0.2541, -0.2541, -0.2541, -0.2541}}, {{-0.9219, -0.9219, -0.9219, -0.9219}}, {{-0.9069, -0.9069, -0.9069, -0.9069}}, {{-0.3089, -0.3089, -0.3089, -0.3089}}, {{-0.3456, -0.3456, -0.3456, -0.3456}}, {{-1.5471, -1.5471, -1.5471, -1.5471}}}}},  // Lambda: r, (r, fs)
    {{{1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026, 1.3026}}, {{{{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}, {{1.3026, 1.3026, 1.3026, 1.3026}}}}},  // conversion: r, (r, fs)
  },
  {
    {1.4940179665528135, {{-0.7402937616912918, -2.79058281438521}}, {{-0.14553681939236265, -0.1348216481090797, -0.21823928939411885}}, 0.07196554284583359},  // K0S
    {1.3816961984929406, {{-0.9361466307814162, -2.1724029874259307}}, {{0.08000994087310911, 0.19683884352489195, 0.025644274124322044}}, 0.18175327639180167},  // Lambda
    {0.18167223399983115, {{-2.333512418596577, -3.015240066721626}}, {{-0.2734437868682242, -0.6409259998835417, -0.7890249725954953}}, 0.0},  // conversion
  },
  {{0.75, 0.375}},  // pointW (K0S, Lambda)
  {{1.75, 1.75}}  // pvWscale (K0S, Lambda)
};

inline constexpr BeamSpotBundle JV_BUNDLE_LR = {
  {{0.0125, 0.0016, 0.7}},  // sigma [cm]
  {
    {{0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 13.0, 20.0, 30.0, 45.0, 70.0, 100.0, 1000000000.0}},  // rBins
    {{0.0, 3.0, 5.0, 8.0, 1000000000.0}},  // fsBins
    2.0,  // fsRmax
    {{{-4.0902, -3.614, -2.7117, -1.8661, -0.6974, 0.4212, 1.2316, 1.9364, 2.347, 2.6702, 2.0147, 2.1112, 2.0641, 1.4651, 1.0479}}, {{{{-4.0984, -4.0186, -3.6392, -1.3104}}, {{-3.7545, -3.4482, -2.2798, -0.354}}, {{-3.0349, -2.5953, -1.8428, -0.8334}}, {{-2.657, -2.2336, -1.6107, -0.7036}}, {{-2.5381, -2.0626, -1.4279, 0.6331}}, {{-1.919, -1.7104, -0.5516, 1.572}}, {{1.2316, 1.2316, 1.2316, 1.2316}}, {{1.9364, 1.9364, 1.9364, 1.9364}}, {{2.347, 2.347, 2.347, 2.347}}, {{2.6702, 2.6702, 2.6702, 2.6702}}, {{2.0147, 2.0147, 2.0147, 2.0147}}, {{2.1112, 2.1112, 2.1112, 2.1112}}, {{2.0641, 2.0641, 2.0641, 2.0641}}, {{1.4651, 1.4651, 1.4651, 1.4651}}, {{1.0479, 1.0479, 1.0479, 1.0479}}}}},  // K0S: r, (r, fs)
    {{{-5.7038, -5.3361, -5.0447, -4.5617, -3.5365, -2.7163, -2.2557, -1.6489, -1.0302, -0.2555, -0.9224, -0.9076, -0.3083, -0.3462, -1.546}}, {{{{-5.5116, -5.7482, -4.5617, -3.427}}, {{-4.9555, -6.7731, -4.3769, -3.8924}}, {{-4.8871, -5.0026, -6.0103, -4.5083}}, {{-4.691, -4.5828, -3.9659, -4.8775}}, {{-4.1144, -4.5554, -2.7771, -2.8629}}, {{-3.9268, -3.5208, -3.5897, -1.7996}}, {{-2.2557, -2.2557, -2.2557, -2.2557}}, {{-1.6489, -1.6489, -1.6489, -1.6489}}, {{-1.0302, -1.0302, -1.0302, -1.0302}}, {{-0.2555, -0.2555, -0.2555, -0.2555}}, {{-0.9224, -0.9224, -0.9224, -0.9224}}, {{-0.9076, -0.9076, -0.9076, -0.9076}}, {{-0.3083, -0.3083, -0.3083, -0.3083}}, {{-0.3462, -0.3462, -0.3462, -0.3462}}, {{-1.546, -1.546, -1.546, -1.546}}}}},  // Lambda: r, (r, fs)
    {{{1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966, 1.2966}}, {{{{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}, {{1.2966, 1.2966, 1.2966, 1.2966}}}}},  // conversion: r, (r, fs)
  },
  {
    {1.4833606061869222, {{-0.7737738166702132, -2.4239277712571514}}, {{-0.1505711703751375, -0.0992784457038006, -0.19411198908668717}}, 0.07228421877023121},  // K0S
    {1.3784526355114348, {{-0.9379505754201185, -2.156054772244486}}, {{0.06759553329562432, 0.2372534445550739, 0.018116028785544826}}, 0.18159359100446928},  // Lambda
    {0.18238204800425245, {{-2.2883719675675334, -3.0147655230275716}}, {{-0.3035465024854113, -0.6184107268257863, -0.7960294786412562}}, 0.0},  // conversion
  },
  {{0.75, 0.375}},  // pointW (K0S, Lambda)
  {{1.75, 1.75}}  // pvWscale (K0S, Lambda)
};

// ---------------------------------------------------------------- settings of the nominal configurations
// JV_<NAME> = value of the configuration key <name> (both beam-spot settings share them); JV_NOV0_<NAME>: the
// value of the fitter without V0 class where it differs. Switches whose other values have no implementation are
// checked by static_assert in jointvtx_algo.h.
constexpr double JV_ACC_COS = 0.5;
constexpr double JV_ACC_DIS_HI = 3.0;
constexpr double JV_ACC_DIS_LO = 0.01;
constexpr double JV_ACC_FSIG2 = 3.0;
constexpr int JV_ACC_MAX_TRK = 8;
constexpr double JV_ACC_SIGL = 0.3;
constexpr bool JV_ACCEPT = true;
constexpr bool JV_ACCEPT_FINAL = true;
constexpr bool JV_ANNEAL_CAND = false;
constexpr double JV_ASSOC_CHI2 = 16.0;
constexpr const char* JV_ASSOC_RULE = "any";
constexpr double JV_CHI2_CUT = 5.0;
constexpr double JV_CHI2_CUT_NEG = 7.0;
constexpr double JV_CHI2_CUT_SV = 16.0;
constexpr const char* JV_CHI2_MODE_PV = "point";
constexpr const char* JV_CHI2_MODE_SV = "vcov";
constexpr double JV_DENS_K = 0.0;
constexpr double JV_EE_WIN = 0.03;
constexpr bool JV_FREEZE_LAST = true;
constexpr int JV_FREEZE_N = 1;
constexpr int JV_GATE_MIN_TPC = 0;
constexpr double JV_GROW_CHI2 = 6.0;
constexpr int JV_GROW_MIN_VDET = 0;
constexpr const char* JV_GROW_RULE = "symmetric";
constexpr bool JV_GROW_SYM_REFIT = false;
constexpr double JV_HARD_MIN = 0.5;
constexpr bool JV_HIT_RADIUS = true;
constexpr double JV_ITER_SIG = 4.0;
constexpr bool JV_ITER_UPDATE_S3 = false;
constexpr bool JV_JOINT = true;
constexpr double JV_KS_WIN = 0.015;
constexpr double JV_LAM_WIN = 0.006;
constexpr bool JV_LIFETIME_PRIOR = true;
constexpr double JV_LOO_EPS = 1.0;
constexpr double JV_MASS_MAX = 6.0;
constexpr bool JV_MERGE_REQUIRE_SHARE = false;
constexpr double JV_MERGE_SHARE = 0.5;
constexpr double JV_MERGE_SIG = 3.0;
constexpr const char* JV_MERGE_WHEN = "each";
constexpr int JV_MIN_SV_TRACKS = 2;
constexpr int JV_N_ITER = 1;
constexpr bool JV_PRIOR_T = false;
constexpr const char* JV_PRIOR_FORM = "lorentz";
constexpr double JV_PRIOR_S = 2.0;
constexpr bool JV_PRIOR_UPDATE = false;
constexpr double JV_PV_TMIN = 4.0;
constexpr double JV_PV_CHI2MAX = AlephReco::kPVChi2Max;
constexpr const char* JV_PV_FINAL = "soft";
constexpr const char* JV_PV_REF = "jetproj_hard";
constexpr double JV_PV_WINDOW_D0 = AlephReco::kPVTrackD0Max;
constexpr double JV_PV_WINDOW_Z0 = AlephReco::kPVTrackZ0MaxBeamspot;
constexpr bool JV_REMOVE_TIGHT_V0 = false;
constexpr bool JV_NOV0_REMOVE_TIGHT_V0 = true;
constexpr const char* JV_S5_MODE = "seeded";
constexpr double JV_SEED_CHI2 = 9.0;
constexpr bool JV_SEED_CROSSINGS = true;
constexpr double JV_SEED_DR = 0.8;
constexpr double JV_SEED_FSIG = 3.0;
constexpr bool JV_SEED_KIN = false;
constexpr const char* JV_SEED_MODE = "signed";
constexpr bool JV_SEED_SAME_JET = false;
constexpr double JV_SEED_SIG = 2.5;
constexpr bool JV_SEED_WINDOW = true;
constexpr double JV_SL_CHI2 = 9.0;
constexpr double JV_SL_FINAL_MIN = 0.55;
constexpr bool JV_SL_FORWARD = true;
constexpr double JV_SL_FSIG = 3.0;
constexpr double JV_SL_SCORE_MIN = -1.0;
constexpr bool JV_SL_SCORE_SEEDED = true;
constexpr bool JV_SL_SYMMETRIC = true;
constexpr double JV_SL_TRIM_CHI2 = 9.0;
constexpr bool JV_SPLIT_CASC = true;
constexpr double JV_SPLIT_DCHI2 = 9.0;
constexpr bool JV_SPLIT_FINAL = true;
constexpr double JV_SPLIT_FSIG = 3.0;
constexpr int JV_SPLIT_NMAX = 7;
constexpr int JV_SPLIT_NMIN = 4;
constexpr double JV_SPLIT_ORDER = 0.98;
constexpr const char* JV_SPLIT_PICK = "valid";
constexpr double JV_SPLIT_SIG = 3.0;
constexpr double JV_STUDENT_NU = 4.0;
constexpr const char* JV_SUPPLY = "vd";
constexpr const char* JV_SV_FINAL = "frozen";
constexpr double JV_SV_FORWARD = 0.8;
constexpr double JV_SV_MIN_P = 0.0;
constexpr int JV_SV_MIN_VDET = 0;
constexpr const char* JV_SV_MODE = "finder";
constexpr double JV_SV_PEN = 0.0;
constexpr double JV_SV_PRIOR_S = 0.0;
constexpr double JV_SV_REFIT_WMIN = 0.0001;
constexpr std::array<double, 3> JV_TEMPS = {{8.0, 2.0, 1.0}};
constexpr double JV_V0_ADMIT = -4.6;
constexpr double JV_V0_CHI2_MAX = 50.0;
constexpr bool JV_V0_CLASS = true;
constexpr bool JV_NOV0_V0_CLASS = false;
constexpr bool JV_V0_CONV = true;
constexpr std::array<double, 2> JV_V0_CONV_HSCALE = {{1.3364, 1.3615}};
constexpr double JV_V0_CONV_RMIN = 4.5;
constexpr double JV_V0_DZ_ABS = 1.0;
constexpr double JV_V0_DZ_SIG = 10.0;
constexpr double JV_V0_HF_CHI2 = 9.0;
constexpr double JV_V0_HF_SIG = 3.0;
constexpr double JV_V0_HIT_SCALE = 3.0;
constexpr bool JV_V0_MASS = true;
constexpr int JV_V0_MAX_COMP = 12;
constexpr double JV_V0_NU = 4.0;
constexpr double JV_V0_OUT_MIN = 0.05;
constexpr bool JV_V0_POINT = true;
constexpr double JV_V0_PRE_KS = 0.15;
constexpr double JV_V0_PRE_LAM = 0.08;
constexpr bool JV_V0_PT_POINT = false;
constexpr std::array<double, 2> JV_V0_QT_MIN = {{0.015, 0.015}};
constexpr double JV_V0_RMAX = 180.0;
constexpr double JV_V0_SM_MIN = 0.0005;
constexpr double JV_V0_TWIN_EXT = 1.5;
constexpr double JV_V0_TWIN_SIG = 3.0;
constexpr bool JV_V0_VETO = false;
constexpr bool JV_NOV0_V0_VETO = true;
constexpr double JV_V0_WIN_KS = 0.06;
constexpr double JV_V0_WIN_LAM = 0.03;
constexpr double JV_VD_R0 = 2.0;
constexpr double JV_VD_FWD_MARGIN = 0.005;
constexpr double JV_VD_GRID_CHI2 = 16.0;
constexpr double JV_VD_GRID_COS = 0.3;
constexpr double JV_VD_GRID_PVCHI2 = 4.0;
constexpr double JV_VD_K0 = 1.0;
constexpr int JV_VD_MAX_STARTS = 15;
constexpr double JV_VD_MEMBER_CHI2 = 9.0;
constexpr int JV_VD_N_CLIMB = 4;
constexpr double JV_VD_TUBE = 1.0;
constexpr double JV_VD_VMIN = 0.3;
constexpr bool JV_VETO_STAGE1_V0 = false;
constexpr double JV_VTX_CHI2NDF = 5.0;
constexpr const char* JV_WEIGHT_FORM = "student";

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
