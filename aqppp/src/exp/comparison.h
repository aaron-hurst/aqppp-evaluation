#pragma once
#include <iostream>

#include "../aqpp/aqpp.h"
#include "../aqpp/assign_budget_for_dimensions.h"
#include "../aqpp/common_content.h"
#include "../aqpp/precompute.h"
#include "../aqpp/sampling.h"
#include "../aqpp/sql_interface.h"
#include "../aqpp/tool.h"
#include "../exp/exp_common.h"

namespace exp_comparison {
class ComparisonExperiment {
 public:
  ComparisonExperiment(SQLHANDLE& sql_connection_handle);
  const int RunExperiment();

 private:
  const std::vector<std::string> AGGREGATIONS_ = {"COUNT", "SUM", "AVG"};
  const std::string DATASET_NAME_ = "uci-household_power_consumption";
  const std::string DB_NAME_ = "uci_household_power_consumption.dbo";
  const std::string TABLE_NAME_ = "household_power_consumption_100k";
  std::string SAMPLE_TABLE_NAME_;
  std::string SUB_SAMPLE_TABLE_NAME_;
  const std::string QUERIES_FILENAME_ =
      "uci-household_power_consumption-N=100.csv";
  const std::string QUERIES_BASE_PATH_ =
      "C:/Users/au686379/OneDrive - Aarhus Universitet/Documents/04 "
      "Research/queries";
  std::string OUTPUT_PATH_;
  std::string QUERIES_PATH_;
  double N_COLUMNS_ = -1;
  int SAMPLE_ROW_NUM_ = -1;
  const double SAMPLE_RATE_ = 0.01;
  const double SUB_SAMPLE_RATE_ = 0.1;
  const double RAND_SEED_ = 1;
  const double CI_INDEX_ = 1.96;
  const double NF_MAX_ITER_ = 1000;
  const bool INIT_DISTINCT_EVEN_ = false;
  const bool isMTL_ = false;
  const int ALL_MTL_POINTS_ = 50000;
  const int EP_PIECE_NUM_ = 20;
  const int N_RUNS_LOAD_SAMPLES_ = 1;
  SQLHANDLE& SQL_CONNECTION_HANDLE_;

  const int LoadQueries(std::vector<aqppp::Condition>& o_queries) const;
  const double ComputePrefixCube(
      const std::vector<std::vector<double>> sample,
      const std::string aggregate_column_name,
      std::vector<std::string> condition_column_names,
      std::vector<std::vector<aqppp::CA>> o_NF_mtl_points,
      aqppp::MTL_STRU o_NF_mtl_res, FILE* log_file) const;
  void WriteParameters(FILE* fp) const;
  static const double PercentageError(const double estimate,
                                      const double exact_value);
};
}  // namespace exp_comparison