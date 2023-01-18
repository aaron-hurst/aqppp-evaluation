// This should contain an entry-point function or class called "RunComparison"
// that is similar to the RunOnce function in comprehensixe_exp, but just
// evaluates the queries that are relevant for my experiments. There should
// probably also be a function that sets the parameters.

#include "comparison.h"

namespace exp_comparison {

ComparisonExperiment::ComparisonExperiment(SQLHANDLE& sql_connection_handle)
    : SQL_CONNECTION_HANDLE_(sql_connection_handle) {
  OUTPUT_PATH_ = "exp_result/comparison/" + DATASET_NAME_ + "/" +
                 "sample_rate_" + std::to_string(SAMPLE_RATE_);
  QUERIES_PATH_ =
      QUERIES_BASE_PATH_ + "/" + DATASET_NAME_ + "/" + QUERIES_FILENAME_;
};

// Load queries
const int ComparisonExperiment::LoadQueries(
    std::vector<aqppp::Condition>& o_queries) const {
  std::ifstream query_file(QUERIES_PATH_);
  std::string line;
  o_queries = std::vector<aqppp::Condition>();  // initialise empty
  if (query_file.is_open()) {
    std::getline(query_file, line);  // skip first line (headers)
    while (std::getline(query_file, line)) {
      // For each row, get the condition column id, low bound and upper bound
      aqppp::Condition condition;
      std::stringstream ss(line);
      std::string str;
      std::getline(ss, str, ',');
      condition.column_id = std::stoi(str);
      std::getline(ss, str, ',');
      condition.lb = std::stod(str);
      std::getline(ss, str, ',');
      condition.ub = std::stod(str);
      o_queries.push_back(condition);
    }
  } else {
    return -1;
  }
  query_file.close();
  return 0;
}

void ComparisonExperiment::WriteParameters(FILE* fp) const {
  fprintf(fp, "---- Databases & files -----\n");
  fprintf(fp, "DATASET_NAME         %s\n", DATASET_NAME_.c_str());
  fprintf(fp, "DB_NAME              %s\n", DB_NAME_.c_str());
  fprintf(fp, "TABLE_NAME           %s\n", TABLE_NAME_.c_str());
  fprintf(fp, "SAMPLE_NAME          %s\n", SAMPLE_TABLE_NAME_.c_str());
  fprintf(fp, "SUB_SAMPLE_NAME      %s\n", SUB_SAMPLE_TABLE_NAME_.c_str());
  fprintf(fp, "QUERIES_FILENAME     %s\n", QUERIES_FILENAME_.c_str());
  fprintf(fp, "OUTPUT_PATH          %s\n", OUTPUT_PATH_.c_str());
  fprintf(fp, "QUERIES_PATH         %s\n", QUERIES_PATH_.c_str());
  fprintf(fp, "----------------------------\n\n");

  fprintf(fp, "---- Sampling & queries ----\n");
  fprintf(fp, "SAMPLE_RATE          %.3f\n", SAMPLE_RATE_);
  fprintf(fp, "SUB_SAMPLE_RATE      %.3f\n", SUB_SAMPLE_RATE_);
  fprintf(fp, "N_RUNS_LOAD_SAMPLES  %i\n", N_RUNS_LOAD_SAMPLES_);
  fprintf(fp, "RAND_SEED            %.3f\n", RAND_SEED_);
  fprintf(fp, "CI_INDEX             %.3f\n", CI_INDEX_);
  fprintf(fp, "SAMPLE_ROW_NUM       %.3f\n", SAMPLE_ROW_NUM_);
  fprintf(fp, "NF_MAX_ITER          %.3f\n", NF_MAX_ITER_);
  fprintf(fp, "INIT_DISTINCT_EVEN   %i\n", INIT_DISTINCT_EVEN_);
  fprintf(fp, "ALL_MTL_POINTS       %i\n", ALL_MTL_POINTS_);
  fprintf(fp, "EP_PIECE_NUM         %i\n", EP_PIECE_NUM_);
  fprintf(fp, "----------------------------\n");
}

// Output files include:
// - parameters (same as loaded from input parameters file)
// - info (timings and high-level statistics)
// - results (query results)
const int ComparisonExperiment::RunExperiment() {
  // Setup
  double t_start = clock();
  FILE *parameters_file, *info_file, *results_file, *log_file;
  aqppp::Tool::MkDirRecursively(OUTPUT_PATH_);  // ensure output path exists
  fopen_s(&parameters_file, (OUTPUT_PATH_ + "/parameters.txt").data(), "w");
  fopen_s(&info_file, (OUTPUT_PATH_ + "/info.txt").data(), "w");
  fopen_s(&results_file, (OUTPUT_PATH_ + "/results.txt").data(), "w");
  fopen_s(&log_file, (OUTPUT_PATH_ + "/log.txt").data(), "w");
  aqppp::SqlInterface::TableColumns table_columns =
      aqppp::SqlInterface::GetTableColumns(SQL_CONNECTION_HANDLE_, TABLE_NAME_);
  std::vector<std::string> column_names = table_columns.all_columns;
  N_COLUMNS_ = column_names.size();
  fprintf(results_file,
          "query_id,aggregation_column,aggregation,estimate_sampling,estimate_"
          "aqppp,ci_sampling,ci_aqppp,"
          "exact_"
          "value,error_sampling,error_"
          "aqppp,time_sampling,time_aqppp,time_exact\n");

  // Export parameters
  ComparisonExperiment::WriteParameters(parameters_file);
  fclose(parameters_file);

  // Create sample database tables
  std::pair<double, double> time_samples_creation =
      aqppp::SqlInterface::CreateDBSamples(
          SQL_CONNECTION_HANDLE_, RAND_SEED_, DB_NAME_, TABLE_NAME_,
          {SAMPLE_RATE_, SUB_SAMPLE_RATE_},
          {SAMPLE_TABLE_NAME_, SUB_SAMPLE_TABLE_NAME_});

  // Load query conditions
  std::vector<aqppp::Condition> queries;
  LoadQueries(queries);

  // Run the queries over each aggregation column
  std::vector<double> sample_load_times = std::vector<double>();
  std::vector<double> sub_sample_load_times = std::vector<double>();
  std::vector<double> compute_prefix_cube_times = std::vector<double>();
  std::vector<double> errors_sampling = std::vector<double>();
  std::vector<double> errors_aqppp = std::vector<double>();
  double total_time_sampling = 0;
  double total_time_aqppp = 0;
  double total_time_exact = 0;
  double n_queries_completed = 0;
  for (int agg_col_id = 0; agg_col_id < N_COLUMNS_; agg_col_id++) {
    std::string aggregate_column_name = table_columns.all_columns[agg_col_id];
    int previous_condition_column_id = -1;
    std::vector<std::vector<double>> sample;
    std::vector<std::vector<double>> small_sample;
    std::vector<std::vector<aqppp::CA>> NF_mtl_points;
    aqppp::MTL_STRU NF_mtl_res;
    // TODO also iterate over each aggregation
    for (int query_id = 0; query_id < queries.size(); query_id++) {
      // Setup
      aqppp::Condition query = queries[query_id];
      std::vector<std::string> condition_column_name = {
          column_names[query.column_id]};

      // If current and previous condition columns are the same, then there is
      // no need to re-load the sample or re-compute the prefix cube
      if (query.column_id != previous_condition_column_id) {
        // Load samples
        // NOTE: This is repeated multiple times to find the median time to load
        // the sample.
        for (int i = 0; i < N_RUNS_LOAD_SAMPLES_; i++) {
          double sample_load_time = aqppp::SqlInterface::ReadDB(
              SQL_CONNECTION_HANDLE_, sample, DB_NAME_, SAMPLE_TABLE_NAME_,
              aggregate_column_name, condition_column_name);
          double sub_sample_load_time = aqppp::SqlInterface::ReadDB(
              SQL_CONNECTION_HANDLE_, small_sample, DB_NAME_,
              SUB_SAMPLE_TABLE_NAME_, aggregate_column_name,
              condition_column_name);
          sample_load_times.push_back(sample_load_time);
          sub_sample_load_times.push_back(sub_sample_load_time);
        }

        // Compute prefix cube
        double compute_prefix_cube_time = ComputePrefixCube(
            sample, aggregate_column_name, condition_column_name, NF_mtl_points,
            NF_mtl_res, log_file);
        compute_prefix_cube_times.push_back(compute_prefix_cube_time);
      }

      // Sampling only
      double t_sampling_start = clock();
      std::pair<double, double> result_sampling = {0, 0};
      // std::pair<double, double> result_sampling =
      //     aqppp::Sampling(PAR.SAMPLE_RATE_, PAR.CI_INDEX_)
      //         .SamplingForSumQuery(sample, queries[query_id]);
      double duration_sampling = (clock() - t_sampling_start) / CLOCKS_PER_SEC;

      // AQP++
      double t_aqppp_start = clock();
      std::pair<double, double> result_aqppp = {0, 0};
      // std::pair<double, double> result_aqppp =
      //     aqppp::Aqpp(PAR.SAMPLE_ROW_NUM, PAR.SAMPLE_RATE_, PAR.CI_INDEX_)
      //         .AqppSumQuery(query_id, info_file, sample, small_sample,
      //                       NF_mtl_points, NF_mtl_res, queries[query_id]);
      double duration_aqppp = (clock() - t_aqppp_start) / CLOCKS_PER_SEC;

      // Exact value... not necessary since I compute this elsewhere already?
      double t_exact_start = clock();
      double exact_value = 0;
      // double exact_value = expDemo::QueryRealValue(
      //     queries[query_id], TABLE_NAME_, SQL_CONNECTION_HANDLE_, PAR,
      //     "sum");
      double duration_exact = (clock() - t_exact_start) / CLOCKS_PER_SEC;

      // Compute errors and store results
      // TODO check if these are actually errors... I think it is just the ratio
      // between the confidence interval width and the exact value
      double error_sampling =
          exact_value > 0 ? result_sampling.second / exact_value : 0;
      double error_aqppp =
          exact_value > 0 ? result_aqppp.second / exact_value : 0;
      errors_sampling.push_back(error_sampling);
      errors_aqppp.push_back(error_aqppp);
      total_time_sampling += duration_sampling;
      total_time_aqppp += duration_aqppp;
      total_time_exact += duration_exact;
      fprintf(results_file, "%d,%d,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", query_id,
        agg_col_id, aggregation, 

              result_sampling.first, result_aqppp.first, result_sampling.second,
              result_aqppp.second, exact_value, error_sampling, error_aqppp, duration_sampling,
              duration_aqppp, duration_exact);

      // Update variables
      n_queries_completed++;
      previous_condition_column_id = query.column_id;
    }
  }
  fclose(results_file);
  fclose(log_file);

  // Export everything
  double time_full_experiment = (clock() - t_start) / CLOCKS_PER_SEC;
  double time_read_sample = aqppp::Tool::get_percentile(sample_load_times, 0.5);
  double time_read_small_sample =
      aqppp::Tool::get_percentile(sub_sample_load_times, 0.5);
  double time_compute_prefix_cube =
      aqppp::Tool::get_percentile(compute_prefix_cube_times, 0.5);
  double time_preprocessing = time_samples_creation.first +
                              time_samples_creation.second +
                              time_compute_prefix_cube;

  double error_sampling_mean = aqppp::Tool::get_avg(errors_sampling);
  double error_aqppp_mean = aqppp::Tool::get_avg(errors_aqppp);
  double error_sampling_median =
      aqppp::Tool::get_percentile(errors_sampling, 0.5);
  double error_aqppp_median = aqppp::Tool::get_percentile(errors_aqppp, 0.5);
  double error_sampling_99p =
      aqppp::Tool::get_percentile(errors_sampling, 0.99);
  double error_aqppp_99p = aqppp::Tool::get_percentile(errors_aqppp, 0.99);

  fprintf(info_file, "Total experiment duration    %.3f s\n",
          time_full_experiment);
  fprintf(info_file, "Sample generation time       %.3f s\n",
          time_samples_creation.first);
  fprintf(info_file, "Subsample generation time    %.3f s\n",
          time_samples_creation.second);
  fprintf(info_file, "\n");

  fprintf(info_file, "Total sample rows            %i\n", SAMPLE_ROW_NUM_);
  fprintf(info_file, "\n");

  fprintf(info_file, "Total pre-processing time    %.3f s\n",
          time_preprocessing);
  fprintf(info_file, "Time to read sample          %.3f s\n", time_read_sample);
  fprintf(info_file, "Time to read sub-sample      %.3f s\n",
          time_read_small_sample);
  fprintf(info_file, "Time to compute BP-cube      %.3f s\n",
          time_compute_prefix_cube);
  fprintf(info_file, "Average query time sampling  %.3f s\n",
          total_time_sampling / n_queries_completed + time_read_sample);
  fprintf(info_file, "Average query time AQP++     %.3f s\n",
          total_time_aqppp / n_queries_completed + time_read_sample +
              time_read_small_sample + time_compute_prefix_cube);
  fprintf(info_file, "Average query time exact     %.3f s\n",
          total_time_exact / n_queries_completed);
  fprintf(info_file, "\n");

  fprintf(info_file, "Average error sampling 	     %.3f %%\n",
          100 * error_sampling_mean);
  fprintf(info_file, "Average error AQP++          %.3f %%\n",
          100 * error_aqppp_mean);
  fprintf(info_file, "Median error sampling 	     %.3f %%\n",
          100 * error_sampling_median);
  fprintf(info_file, "Median error AQP++ 	         %.3f %%\n",
          100 * error_aqppp_median);
  fprintf(info_file, "99th percentile sampling     %.3f %%\n",
          100 * error_sampling_99p);
  fprintf(info_file, "99th percentile AQP++ 	     %.3f %%\n",
          100 * error_aqppp_99p);
  fprintf(info_file, "\n");
  fclose(info_file);
  std::cout << "Comparison experiment completed." << std::endl;
  return 0;
}

const double ComparisonExperiment::ComputePrefixCube(
    const std::vector<std::vector<double>> sample,
    const std::string aggregate_column_name,
    const std::vector<std::string> condition_column_names,
    std::vector<std::vector<aqppp::CA>> o_NF_mtl_points,
    aqppp::MTL_STRU o_NF_mtl_res, FILE* log_file) const {
  double t_start = clock();

  // Compute CAsample
  std::vector<std::vector<aqppp::CA>> CAsample =
      std::vector<std::vector<aqppp::CA>>();
  aqppp::Tool::TransSample(sample, CAsample);

  // Allocate space budget for prefix cube
  double t_prepare_start = clock();
  std::vector<int> mtl_nums;
  aqppp::AssignBudgetForDimensions(SAMPLE_RATE_, ALL_MTL_POINTS_, EP_PIECE_NUM_,
                                   SAMPLE_ROW_NUM_, CI_INDEX_, NF_MAX_ITER_,
                                   INIT_DISTINCT_EVEN_, false)
      .AssignBudget(CAsample, mtl_nums);

  // Select the critical points for the prefix cube
  o_NF_mtl_points = std::vector<std::vector<aqppp::CA>>();
  std::vector<double> max_errs;
  std::vector<int> iter_nums;
  aqppp::HillClimbing(SAMPLE_ROW_NUM_, SAMPLE_RATE_, CI_INDEX_, NF_MAX_ITER_,
                      INIT_DISTINCT_EVEN_)
      .ChoosePoints(CAsample, mtl_nums, o_NF_mtl_points, max_errs, iter_nums);

  // Log max errors
  for (int i = 0; i < CAsample.size(); i++) {
    std::cout << "dim: " << i << " mtl_num: " << mtl_nums[i] << std::endl;
    fprintf(log_file,
            "dim:%d\tdistinct_value:%i\tmtl_num:%d\tmax_err:%f\tclimb_iter_num:"
            "%d\n",
            i, CAsample[i].size(), mtl_nums[i], max_errs[i], iter_nums[i]);
  }

  // Compute prefix cube
  // NOTE: This is not used in the "comprehensive experiment"
  o_NF_mtl_res = aqppp::MTL_STRU();
  if (isMTL_) {
    aqppp::Precompute(DB_NAME_, TABLE_NAME_, aggregate_column_name,
                      condition_column_names)
        .GetPrefixSumCube(o_NF_mtl_points, SQL_CONNECTION_HANDLE_, o_NF_mtl_res,
                          "sum");
  }

  return (clock() - t_start) / CLOCKS_PER_SEC;
}
}  // namespace exp_comparison