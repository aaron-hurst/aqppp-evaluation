#pragma once
#include <assert.h>

#include <algorithm>
#include <mutex>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include "common_content.h"
#include "sql_interface.h"
#include "tool.h"

#define BUFFER_COLULMS 10
#define BUFFER_ROWS 100

namespace aqppp {
class Precompute {
 private:
  struct Par {
    std::string DB_NAME;
    std::string TABLE_NAME;
    std::string AGGREGATE_NAME;
    std::vector<std::string> CONDITION_NAMES;
  } PAR;

 public:
  Precompute(std::string db_name, std::string table_name, std::string agg_name,
             std::vector<std::string> condition_names);
  static int MyLowerBound(const std::vector<CA>& cur_col, const CA& key);
  void InitCube(int a, std::vector<int>& cur_indx, MTL_STRU& o_mtl_data,
                const std::vector<std::vector<CA>>& mtl_points);
  double GetPrefixSumCube(const std::vector<std::vector<CA>>& mtl_points,
                          SQLHANDLE& sqlconnectionhandle, MTL_STRU& o_mtl_res,
                          std::string query_type,
                          const DistId& distinct_ids = DistId());
  double GetPrefixSumCube(const std::vector<int>& size_of_dimension,
                          const MTL_STRU& mtl_data, MTL_STRU& o_mtl_res);

 private:
  double ComputeDataCube(const std::vector<std::vector<CA>>& mtl_points,
                         MTL_STRU& o_mtl_data,
                         const std::string query_type, const SQLHANDLE& SQL_CONNECTION_HANDLE_,
                         const DistId& distinct_ids = DistId());
  void ComputeSumCube(int a, int b, const std::vector<int>& size_of_dimension,
                      std::vector<int>& cur_ind, MTL_STRU& o_mtl_res);
  static void ComputeDataCubeOneThread(
      int thread_id, int start_row, int end_row, int CONDITION_DIM,
      const std::vector<std::vector<CA>>& mtl_points,
      double (&table)[BUFFER_COLULMS][BUFFER_ROWS], std::vector<double>& cube,
      std::vector<int>& cube_sizes, std::vector<std::mutex>& cube_locks);
  static void ReadColumnOneThread(SQLHANDLE& sqlstatementhandle);
};
}  // namespace aqppp