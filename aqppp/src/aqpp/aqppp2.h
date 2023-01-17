#pragma once
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "common_content.h"
#include "hill_climbing.h"
#include "sampling.h"
#include "sql_interface.h"
#include "tool.h"

namespace aqppp {
class AQPpp {
 private:
  struct Par {
    int SAMPLE_ROW_NUM;
    double SAMPLE_RATE;
    double CI_INDEX;
  } PAR;

  typedef std::vector<Condition> RGs;

  void GenAllRanges(const std::vector<std::vector<CA>>& mtl_points,
                    const std::vector<int> condition_columns,
                    const std::vector<Condition> conditions,
                    std::vector<RGs>& o_all_rgs);

  double ComputeDifference4Sum(const std::vector<std::vector<double>>& sample,
                               const int aggregate_column,
                               const std::vector<int> condition_columns,
                               const std::vector<Condition>& conditions,
                               const std::vector<Condition>& mtl_choices);

  void FillMaterializeChoices(
      const int col, const std::vector<std::vector<double>>& small_sample,
      Query query, const std::vector<RGs>& all_rgs, double& o_min_ci,
      std::vector<Condition>& mtl_choices,
      std::vector<Condition>& final_mtl_choices);

  void FinalComputeDifference(
      const std::vector<std::vector<double>>& sample,
      const int aggregate_column, const std::vector<int> condition_columns,
      const std::vector<Condition>& conditions,
      const std::vector<Condition>& mtl_choices,
      std::vector<std::pair<double, double>>& final_diff_res);

  void ComputeMquery(const int col, std::vector<int>& inds,
                     const MTL_STRU& mtl_res,
                     const std::vector<Condition>& mtl_choices, double& o_sum,
                     int minus_value);

  Condition GenCondition(const std::vector<CA>& mtl_points, int lbid, int ubid);

 public:
  AQPpp(int sample_row_num, double sample_rate, double ci_index);

  const std::pair<double, double> RunQuery(
      const int query_id, Query query,
      const std::vector<std::vector<double>>& sample,
      const std::vector<std::vector<double>>& small_sample,
      const std::vector<std::vector<CA>>& mtl_points, const MTL_STRU& mtl_res,
      FILE* log_file);
};
}  // namespace aqppp
