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
  struct Parameters {
    std::string DATASET_NAME = "uci-household_power_consumption";
    std::string DB_NAME = "uci_household_power_consumption.dbo";
    std::string TABLE_NAME = "household_power_consumption";
    std::string SAMPLE_NAME = "household_power_consumption_sample";
    std::string SUB_SAMPLE_NAME = "household_power_consumption_sub_sample";
    std::string QUERIES_FILENAME = "uci-household_power_consumption-N=100.csv";
    std::string QUERIES_BASE_PATH =
        "C:/Users/au686379/OneDrive - Aarhus Universitet/Documents/04 "
        "Research/queries";
    std::string OUTPUT_PATH;
    std::string QUERIES_PATH;
    double N_COLUMNS = -1;
    int SAMPLE_ROW_NUM = -1;
    double SAMPLE_RATE = 0.01;
    double SUB_SAMPLE_RATE = 0.1;
    double RAND_SEED = 1;
    double CI_INDEX = 1.96;
    double NF_MAX_ITER = 1000;
    bool INIT_DISTINCT_EVEN = false;
    bool isMTL = false;
    int ALL_MTL_POINTS = 50000;
    int EP_PIECE_NUM = 20;
  };

  static const int RunExperiment(SQLHANDLE& sqlconnectionhandle);

 private:
  static const std::vector<std::string> AGGREGATIONS;

  static const Parameters LoadParameters(void);
  static const int LoadQueries(const std::string query_filepath,
                               std::vector<aqppp::Query>& o_queries,
                               const int n_columns);
  static const std::pair<double, double> ReadSamples(
      SQLHANDLE& sql_connection_handle, const std::string db_name,
      const std::string sample_table_name,
      const std::string small_sample_table_name,
      std::vector<std::vector<double>>& o_sample,
      std::vector<std::vector<double>>& o_small_sample, const int n_runs = 10);
  static const double ReadDB(SQLHANDLE& sql_connection_handle,
                             std::vector<std::vector<double>>& o_table,
                             const std::string db_name,
                             const std::string table_name);
  static void WriteParameters(FILE* fp, Parameters par);
};
}  // namespace exp_comparison