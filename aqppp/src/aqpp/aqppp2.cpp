#pragma once
#include "aqppp2.h"

namespace aqppp {

AQPpp::AQPpp(int sample_row_num, double sample_rate, double ci_index) {
  this->PAR.SAMPLE_ROW_NUM = sample_row_num;
  this->PAR.SAMPLE_RATE = sample_rate;
  this->PAR.CI_INDEX = ci_index;
}

const std::pair<double, double> AQPpp::RunQuery(
    const int query_id, aqppp::Query query,
    const std::vector<std::vector<double>>& sample,
    const std::vector<std::vector<double>>& small_sample,
    const std::vector<std::vector<CA>>& mtl_points, const MTL_STRU& mtl_res,
    FILE* log_file) {
  // Setup
  const int DIM = mtl_points.size();
  std::vector<RGs> all_rgs = std::vector<RGs>();
  GenAllRanges(mtl_points, query.condition_columns, query.conditions, all_rgs);

  // In case cannot find mtl range related to user query.
  for (int i = 0; i < all_rgs.size(); i++) {
    if (all_rgs[i].size() < 1) {
      fprintf(log_file,
              "Failed to generate ranges for query %i. Returning sampling "
              "result.\n");
      return Sampling(this->PAR.SAMPLE_RATE, this->PAR.CI_INDEX)
          .SamplingForSumQuery(
              sample, query.conditions);  // TODO this will need re-work after I
                                          // re-work sampling
    }
  }

  // ???
  std::vector<Condition> mtl_choices = std::vector<Condition>(all_rgs.size());
  std::vector<Condition> final_mtl_choices = std::vector<Condition>();
  double min_ci = FLT_MAX;
  FillMaterializeChoices(0, small_sample, query, all_rgs, min_ci, mtl_choices,
                         final_mtl_choices);
  std::vector<std::pair<double, double>> final_diff_res =
      std::vector<std::pair<double, double>>();
  FinalComputeDifference(sample, query.aggregate_column,
                         query.condition_columns, query.conditions,
                         final_mtl_choices, final_diff_res);
  std::pair<double, double> spl_sum_ci = final_diff_res[0];
  std::pair<double, double> hybrid_sum_ci = final_diff_res[1];
  if (spl_sum_ci.second < hybrid_sum_ci.second) return spl_sum_ci;

  std::vector<int> inds = std::vector<int>(DIM);
  if (mtl_res.size() == 0) {
    hybrid_sum_ci.first = -1;
    return hybrid_sum_ci;
  }
  double sum2 = 0;
  ComputeMquery(0, inds, mtl_res, final_mtl_choices, sum2, 0);
  hybrid_sum_ci.first += sum2;
  return hybrid_sum_ci;
}

void AQPpp::GenAllRanges(const std::vector<std::vector<CA>>& mtl_points,
                         const std::vector<int> condition_columns,
                         const std::vector<Condition> conditions,
                         std::vector<RGs>& o_all_rgs) {
  o_all_rgs = std::vector<RGs>(mtl_points.size());
  assert(mtl_points.size() ==
         7);  // TEMPORARY CODE: I assume that mtl_points has a column for each
              // dimension... if it has one less (i.e. excludes the aggregation
              // dimension, then I need to rethink some code)
  for (int i = 0; i < mtl_points.size(); i++) {
    RGs rgs = RGs();
    CA user_lb = CA();
    CA user_ub = CA();
    std::vector<int>::const_iterator itr =
        std::find(condition_columns.begin(), condition_columns.end(), i);
    if (itr != condition_columns.end()) {
      int ci = std::distance(condition_columns.begin(), itr);
      user_lb.condition_value = conditions[ci].lb;
      user_ub.condition_value = conditions[ci].ub;
    } else {
      user_lb.condition_value = -FLT_MAX;
      user_ub.condition_value = FLT_MAX;
    }

    // Get bound IDs
    int lb_lbid = -1;
    int lb_ubid = -1;
    int ub_lbid = -1;
    int ub_ubid = -1;
    std::vector<CA>::const_iterator it;
    it = upper_bound(mtl_points[i].begin(), mtl_points[i].end(), user_lb,
                     CA_compare);  // less & equal
    if (it > mtl_points[i].begin() &&
        it <= mtl_points[i].end())  // note that the position of '='.
    {
      lb_lbid = it - mtl_points[i].begin() - 1;
    } else
      lb_lbid = -1;

    it = lower_bound(mtl_points[i].begin(), mtl_points[i].end(), user_lb,
                     CA_compare);  // great & equal
    if (it >= mtl_points[i].begin() && it < mtl_points[i].end()) {
      lb_ubid = it - mtl_points[i].begin();
    } else
      lb_ubid = -1;

    it = upper_bound(mtl_points[i].begin(), mtl_points[i].end(), user_ub,
                     CA_compare);  // less & equal
    if (it > mtl_points[i].begin() && it <= mtl_points[i].end()) {
      ub_lbid = it - mtl_points[i].begin() - 1;
    } else
      ub_lbid = -1;

    it = lower_bound(mtl_points[i].begin(), mtl_points[i].end(), user_ub,
                     CA_compare);  // great & equal
    if (it >= mtl_points[i].begin() && it < mtl_points[i].end()) {
      ub_ubid = it - mtl_points[i].begin();
    } else
      ub_ubid = -1;

    // Case that only one bound is found and cannot construct a mtl range.
    // Especially when do equal query.
    if (lb_lbid == lb_ubid && lb_lbid != -1) {
      lb_lbid--;  // note that -1 is valid for lbid.
    }

    if (ub_lbid == ub_ubid && ub_lbid != -1) {
      if (ub_ubid + 1 <= mtl_points[i].size() - 1)
        ub_ubid++;
      else
        ub_lbid--;
    }

    // Update ranges
    if (lb_lbid < ub_lbid)
      rgs.push_back(GenCondition(mtl_points[i], lb_lbid, ub_lbid));
    if (lb_lbid < ub_ubid)
      rgs.push_back(GenCondition(mtl_points[i], lb_lbid, ub_ubid));
    if (lb_ubid != -1 && lb_ubid < ub_lbid)
      rgs.push_back(GenCondition(mtl_points[i], lb_ubid, ub_lbid));
    if (lb_ubid != -1 && lb_ubid < ub_ubid)
      rgs.push_back(GenCondition(mtl_points[i], lb_ubid, ub_ubid));
    o_all_rgs[i] = rgs;
  }
  return;
}

// Used for function GenAllRanges
Condition AQPpp::GenCondition(const std::vector<CA>& mtl_points, int lbid,
                              int ubid) {
  // lbid<0 is valid, but ubid<0 is invalid. but when ubid<0 isn't a valid RG,
  // and shouldn't GenCondition. NOTE: We need to store mtl_points.size() as an
  // int, otherwise comparison to lbin fails for lbid<0 due to mtl_points.size()
  // being unsigned.
  const int mtl_points_size = mtl_points.size();
  assert(ubid >= 0 && ubid < mtl_points_size);
  assert(lbid < mtl_points_size);
  Condition d = Condition();
  d.lb_id = lbid;
  d.ub_id = ubid;
  d.lb = lbid < 0 ? -FLT_MAX : mtl_points[lbid].condition_value;
  d.ub = mtl_points[ubid].condition_value;
  return d;
}

// note that this function doesn't consider the case that don't find the mtl
// range. o_min_ci and mtl_choice should be init outside. It choice the mtl
// solution in small sample, and need to further
void AQPpp::FillMaterializeChoices(
    const int col, const std::vector<std::vector<double>>& small_sample,
    Query query, const std::vector<RGs>& all_rgs, double& o_min_ci,
    std::vector<Condition>& mtl_choices,
    std::vector<Condition>& final_mtl_choices) {
  if (col >= all_rgs.size()) {
    double ci = ComputeDifference4Sum(small_sample, query.aggregate_column,
                                      query.condition_columns, query.conditions,
                                      mtl_choices);
    if (ci < o_min_ci) {
      o_min_ci = ci;
      final_mtl_choices = mtl_choices;
    }
    return;
  }

  const int CUR_SIZE = all_rgs[col].size();
  for (int i = 0; i < CUR_SIZE; i++) {
    mtl_choices[col] = all_rgs[col][i];
    FillMaterializeChoices(col + 1, small_sample, query, all_rgs, o_min_ci,
                           mtl_choices, final_mtl_choices);
  }
}

// used for function FillMaterializeChoices.
// should the row num is the real row num or distinct value of row num?
double AQPpp::ComputeDifference4Sum(
    const std::vector<std::vector<double>>& sample, const int aggregate_column,
    const std::vector<int> condition_columns,
    const std::vector<Condition>& conditions,
    const std::vector<Condition>& mtl_choices) {
  double sum = 0;
  double sum2 = 0;
  const int ROW_NUM = sample[0].size();
  const int USER_CONDITION_DIM = conditions.size();
  // Iterate over rows of the sample
  for (int ri = 0; ri < ROW_NUM; ri++) {
    bool umeet = true;
    bool mmeet = true;

    // Check if row satisfies query conditions
    for (int ci = 0; ci < conditions.size(); ci++) {
      double con_data = sample[condition_columns[ci]][ri];
      if (DoubleLess(con_data, conditions[ci].lb) ||
          DoubleGreater(con_data, conditions[ci].ub)) {
        umeet = false;
        break;
      }
    }

    // Check if row ... has something to do with MTL?
    for (int ci = 1; ci < mtl_choices.size() + 1; ci++) {
      double con_data = sample[ci][ri];
      if (DoubleLeq(con_data, mtl_choices[ci - 1].lb) ||
          DoubleGreater(con_data, mtl_choices[ci - 1].ub)) {
        mmeet = false;
        break;
      }
    }

    // Increment sum and sum of squares as appropriate
    double udata = 0;
    double mdata = 0;
    if (umeet) udata = sample[aggregate_column][ri];
    if (mmeet) mdata = sample[aggregate_column][ri];
    double data = udata - mdata;
    sum += data;
    sum2 += data * data;
  }

  // Compute confidence interval
  double variance = 0;
  variance = sum2 / ROW_NUM - (sum / ROW_NUM) * (sum / ROW_NUM);
  double ci = this->PAR.CI_INDEX * sqrt(variance / ROW_NUM) * ROW_NUM /
              this->PAR.SAMPLE_RATE;
  return ci;
}

// when got mtl choice, the final compute the diffrence to chose sampling or
// hybrid
void AQPpp::FinalComputeDifference(
    const std::vector<std::vector<double>>& sample, const int aggregate_column,
    const std::vector<int> condition_columns,
    const std::vector<Condition>& conditions,
    const std::vector<Condition>& mtl_choices,
    std::vector<std::pair<double, double>>& final_diff_res) {
  double hsum = 0;  // h for hybrid
  double hsum2 = 0;
  double ssum = 0;  // s for sampling
  double ssum2 = 0;
  const int ROW_NUM = sample[0].size();
  for (int ri = 0; ri < ROW_NUM; ri++) {
    bool umeet = true;
    bool mmeet = true;
    for (int ci = 0; ci < conditions.size(); ci++) {
      double con_data = sample[condition_columns[ci]][ri];
      if (DoubleLess(con_data, conditions[ci].lb) ||
          DoubleGreater(con_data, conditions[ci].ub)) {
        umeet = false;
        break;
      }
    }
    for (int ci = 1; ci < mtl_choices.size() + 1; ci++) {
      double con_data = sample[ci][ri];
      if (DoubleLeq(con_data, mtl_choices[ci - 1].lb) ||
          DoubleGreater(con_data, mtl_choices[ci - 1].ub)) {
        mmeet = false;
        break;
      }
    }

    double udata = 0;
    double mdata = 0;
    if (umeet) udata = sample[aggregate_column][ri];
    if (mmeet) mdata = sample[aggregate_column][ri];
    double data = udata - mdata;
    hsum += data;  // h for hybrid.
    hsum2 += data * data;

    ssum += udata;  // s for sampling.
    ssum2 += udata * udata;
  }
  double hvariance = hsum2 / ROW_NUM - (hsum / ROW_NUM) * (hsum / ROW_NUM);
  double hci = this->PAR.CI_INDEX * sqrt(hvariance / ROW_NUM) * ROW_NUM /
               this->PAR.SAMPLE_RATE;
  double svariance = ssum2 / ROW_NUM - (ssum / ROW_NUM) * (ssum / ROW_NUM);
  double sci = this->PAR.CI_INDEX * sqrt(svariance / ROW_NUM) * ROW_NUM /
               this->PAR.SAMPLE_RATE;

  final_diff_res = std::vector<std::pair<double, double>>();
  final_diff_res.push_back(
      std::pair<double, double>(ssum / this->PAR.SAMPLE_RATE, sci));
  final_diff_res.push_back(
      std::pair<double, double>(hsum / this->PAR.SAMPLE_RATE, hci));
  return;
}

//*********just use ub-lb rather ub-(lb-1) to get the range sum, because the mtl
// point is sparse. In this case, the result is sum(lb_id,ub_id], note '(]'.
// ans need to be inti 0.
// inds need to be init vector<int>(mtl_res.dim).
// mquery means mtl query.
// cpt given mtl_query in mtl result.
// mtl_choice[i] is the id of i-dimension mtl points in mtl result.
void AQPpp::ComputeMquery(const int col, std::vector<int>& inds,
                          const MTL_STRU& mtl_res,
                          const std::vector<Condition>& mtl_choices,
                          double& o_sum, int minus_value) {
  if (col >= mtl_choices.size())  // dim=mtl_choice.size
  {
    // the following is one term of the final ans.
    int t = 1;
    for (int i = 0; i < inds.size(); i++) {
      if (inds[i] < 0) return;
      if (inds[i] == mtl_choices[i].lb_id - minus_value) t = t * (-1);
    }
    if (mtl_res.count(inds) == 1) {
      o_sum += t * mtl_res.at(inds);
    }
    return;
  }

  for (int i = 0; i < 2; i++) {
    if (i == 0)
      inds[col] = mtl_choices[col].lb_id - minus_value;
    else
      inds[col] = mtl_choices[col].ub_id;
    ComputeMquery(col + 1, inds, mtl_res, mtl_choices, o_sum, minus_value);
  }

  return;
}

}  // namespace aqppp
