#pragma once
#include <iostream>

#include "../aqpp/aqpp.h"
#include "../aqpp/assign_budget_for_dimensions.h"
#include "../aqpp/precompute.h"
#include "../aqpp/sampling.h"
#include "../aqpp/sql_interface.h"
#include "../aqpp/tool.h"
#include "../exp/exp_common.h"

namespace exp_comparison {
class ComparisonExperiment {
 private:
  struct Parameters {
    std::string DATASET_NAME;
    std::string DB_NAME;
    std::string TABLE_NAME;
    std::string SAMPLE_NAME;
    std::string SUB_SAMPLE_NAME;
    std::string OUTPUT_PATH;
    double SAMPLE_RATE;       // 0.01
    double SUB_SAMPLE_RATE;   // 0.1
    double RAND_SEED;         // 1
    double CI_INDEX;          // = 1.96;
    double SAMPLE_ROW_NUM;    // -1
    double NF_MAX_ITER;       // 1000z
    bool INIT_DISTINCT_EVEN;  // false
    // TODO revise, include all parameters, including those from PAR and ExpPar
    // bool USE_DB_INDEX;
    // bool DROP_INDEX_BEFORE_CREATE;
    // bool isMTL;
    // bool EXP_UNIFORM;
    // bool QUERY_DB_REAL_VALUE;
    // bool QUERY_DB_SELECTIVELY;
    // bool COVER_DIRECT_RESULT_FILE;  // generally this should be false.
  };

  static const Parameters LoadParameters(void);

 public:
  static const int RunExperiment(SQLHANDLE& sqlconnectionhandle);
};
}  // namespace exp_comparison